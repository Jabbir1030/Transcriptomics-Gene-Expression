"""Core analysis: TLS/Stroma/ECO scores, response AUCs, survival, subtypes. All real data."""
import os, sys, warnings
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.cluster import KMeans
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
sys.path.insert(0, "/home/user/npj_immunotherapy_paper/code")
from signatures import TLS, STROMA, IFNG6, CD8EFF, CHECKPOINTS, AXES

warnings.filterwarnings("ignore")
PROC = "/home/user/npj_immunotherapy_paper/data/processed"
TAB = "/home/user/npj_immunotherapy_paper/results/tables"
os.makedirs(TAB, exist_ok=True)

COHORTS = {
    "TCGA-SKCM": {"cancer": "Melanoma", "type": "TCGA"},
    "TCGA-LUAD": {"cancer": "NSCLC", "type": "TCGA"},
    "TCGA-BLCA": {"cancer": "Bladder", "type": "TCGA"},
    "GSE78220": {"cancer": "Melanoma", "type": "ICI", "drug": "anti-PD-1"},
    "GSE91061": {"cancer": "Melanoma", "type": "ICI", "drug": "nivolumab"},
    "GSE126044": {"cancer": "NSCLC", "type": "ICI", "drug": "anti-PD-1"},
    "GSE135222": {"cancer": "NSCLC", "type": "ICI", "drug": "anti-PD-1/PD-L1"},
    "GSE207422": {"cancer": "NSCLC", "type": "ICI", "drug": "anti-PD-1+chemo (neoadj)"},
    "IMvigor210": {"cancer": "Bladder", "type": "ICI", "drug": "atezolizumab"},
}

def zscore(df):
    m = df.mean(axis=1); s = df.std(axis=1).replace(0, np.nan)
    return df.sub(m, axis=0).div(s, axis=0)

def score(z, genes):
    g = [x for x in genes if x in z.index]
    return z.loc[g].mean(axis=0), len(g)

all_scores, resp_stats, surv_stats = {}, [], []
for cohort, meta in COHORTS.items():
    expr = pd.read_csv(f"{PROC}/{cohort}_expr.csv", index_col=0)
    clin = pd.read_csv(f"{PROC}/{cohort}_clin.csv", dtype={"sample_id": str})
    clin = clin.set_index("sample_id")
    common = [c for c in expr.columns if str(c) in set(clin.index)]
    expr = expr[common]
    clin = clin.loc[[str(c) for c in common]]
    z = zscore(expr)
    tls, n_tls = score(z, TLS)
    stro, n_stro = score(z, STROMA)
    eco = tls - stro
    ifng, _ = score(z, IFNG6)
    cd8, _ = score(z, CD8EFF)
    cd274 = z.loc["CD274"] if "CD274" in z.index else pd.Series(np.nan, index=expr.columns)
    s = pd.DataFrame({"TLS": tls, "STROMA": stro, "ECO": eco, "IFNG6": ifng, "CD8EFF": cd8, "CD274": cd274})
    s["cohort"] = cohort; s["cancer"] = meta["cancer"]
    all_scores[cohort] = (s, clin)
    print(f"[{cohort}] n={len(common)} TLS_genes={n_tls} STROMA_genes={n_stro}", flush=True)
    # response stats
    if "response01" in clin.columns:
        cc = clin["response01"].dropna()
        for name, vals in [("ECO", eco), ("TLS", tls), ("STROMA", stro), ("IFNG6", ifng), ("CD8EFF", cd8), ("CD274", cd274)]:
            v = vals.loc[[c for c in common if str(c) in set(cc.index.map(str))]]
            y = cc.loc[[str(c) for c in v.index]].astype(int)
            if y.nunique() < 2:
                continue
            auc = roc_auc_score(y, v)
            # for STROMA (negative predictor) also report inverted? keep raw; note direction
            try:
                p = mannwhitneyu(v[y == 1], v[y == 0], alternative="two-sided").pvalue
            except Exception:
                p = np.nan
            resp_stats.append({"cohort": cohort, "cancer": meta["cancer"], "n": len(y),
                               "n_resp": int(y.sum()), "metric": name, "AUC": round(float(auc), 3),
                               "MWU_p": float(p),
                               "median_resp": round(float(v[y == 1].median()), 3),
                               "median_nonresp": round(float(v[y == 0].median()), 3)})
    # survival stats (median split on ECO)
    if "time" in clin.columns and "event" in clin.columns:
        t = pd.to_numeric(clin["time"], errors="coerce"); e = pd.to_numeric(clin["event"], errors="coerce")
        ok = t.notna() & e.notna() & (t > 0)
        if ok.sum() > 20:
            med = eco.median()
            grp = (eco >= med).astype(int)
            d = pd.DataFrame({"t": t, "e": e, "grp": grp})[ok]
            lr = logrank_test(d[d.grp == 1]["t"], d[d.grp == 0]["t"], d[d.grp == 1]["e"], d[d.grp == 0]["e"])
            cph = CoxPHFitter().fit(d[["t", "e", "grp"]], "t", "e")
            hr = float(cph.hazard_ratios_["grp"]); ci = cph.confidence_intervals_.loc["grp"].tolist()
            surv_stats.append({"cohort": cohort, "cancer": meta["cancer"], "endpoint": "PFS" if cohort == "GSE135222" else "OS",
                               "n": int(ok.sum()), "events": int(d.e.sum()),
                               "HR_high_vs_low": round(hr, 3), "CI_low": round(ci[0], 3), "CI_high": round(ci[1], 3),
                               "logrank_p": float(lr.p_value)})

resp_df = pd.DataFrame(resp_stats); surv_df = pd.DataFrame(surv_stats)
resp_df.to_csv(f"{TAB}/response_AUCs.csv", index=False)
surv_df.to_csv(f"{TAB}/survival_stats.csv", index=False)
print("\n=== RESPONSE AUCs ==="); print(resp_df.to_string())
print("\n=== SURVIVAL ==="); print(surv_df.to_string())

# pooled ICI analysis (within-cohort standardized ECO)
pool = []
for cohort, meta in COHORTS.items():
    if meta["type"] != "ICI":
        continue
    s, clin = all_scores[cohort]
    if "response01" not in clin.columns:
        continue
    cc = clin["response01"].dropna()
    e = s["ECO"].loc[[c for c in s.index if str(c) in set(cc.map(str).index)]]
    y = cc.loc[[str(c) for c in e.index]].astype(int)
    pool.append(pd.DataFrame({"ECO": (e - e.mean()) / e.std(), "y": y.values, "cohort": cohort}))
pool = pd.concat(pool)
pooled_auc = roc_auc_score(pool.y, pool.ECO)
print(f"\nPOOLED ICI: n={len(pool)} responders={int(pool.y.sum())} pooled_AUC={pooled_auc:.3f}")
pool.to_csv(f"{TAB}/pooled_ICI_ECO.csv", index=False)

# ecosystem subtypes: KMeans on TLS/STROMA pooled TCGA, then map ICI samples to centroids
tcga = pd.concat([all_scores[c][0][["TLS", "STROMA"]] for c in ["TCGA-SKCM", "TCGA-LUAD", "TCGA-BLCA"]])
km = KMeans(n_clusters=3, n_init=20, random_state=42).fit(tcga[["TLS", "STROMA"]].values)
cent = km.cluster_centers_
order = np.lexsort((cent[:, 1], -cent[:, 0]))  # high TLS first
label_map = {order[0]: "TLS-high", order[1]: "Intermediate", order[2]: "Stromal"}
print("centroids (TLS, STROMA):", np.round(cent, 2), "->", [label_map[i] for i in range(3)])
subtype_resp = []
for cohort, meta in COHORTS.items():
    s, clin = all_scores[cohort]
    cl = km.predict(s[["TLS", "STROMA"]].values)
    s["subtype"] = [label_map[i] for i in cl]
    s.to_csv(f"{TAB}/scores_{cohort}.csv")
    if "response01" in clin.columns:
        cc = clin["response01"].dropna()
        for st in ["TLS-high", "Intermediate", "Stromal"]:
            idx = s[s.subtype == st].index
            yy = cc.loc[[str(c) for c in idx if str(c) in set(cc.index.map(str))]]
            if len(yy) > 0:
                subtype_resp.append({"cohort": cohort, "cancer": meta["cancer"], "subtype": st,
                                     "n": len(yy), "resp_rate": round(float(yy.mean()), 3)})
sub_df = pd.DataFrame(subtype_resp)
sub_df.to_csv(f"{TAB}/subtype_response.csv", index=False)
print("\n=== SUBTYPE RESPONSE RATES ==="); print(sub_df.to_string())
# save centroids
pd.DataFrame(cent, columns=["TLS", "STROMA"], index=[label_map[i] for i in range(3)]).to_csv(f"{TAB}/subtype_centroids.csv")
print("\nSaved all tables to", TAB)
