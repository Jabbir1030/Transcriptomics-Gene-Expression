"""Refinements: correct subtype labels, HR CIs, bootstrap AUC CIs, ECO-high/low OR, multivariate Cox, ECO+TMB combo."""
import os, warnings
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from scipy.stats import fisher_exact
from lifelines import CoxPHFitter
warnings.filterwarnings("ignore")

TAB = "/home/user/npj_immunotherapy_paper/results/tables"
PROC = "/home/user/npj_immunotherapy_paper/data/processed"
rng = np.random.default_rng(42)

def boot_auc(y, v, n=2000):
    y = np.asarray(y); v = np.asarray(v)
    aucs = []
    for _ in range(n):
        idx = rng.integers(0, len(y), len(y))
        if len(np.unique(y[idx])) < 2:
            continue
        aucs.append(roc_auc_score(y[idx], v[idx]))
    return float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5))

# ---- 1. fix survival CIs (exp of coef CI) + multivariate Cox ----
try:
    surv_old = pd.read_csv(f"{TAB}/survival_stats.csv")
    keep_cols = [c for c in ["cohort", "logrank_p", "endpoint"] if c in surv_old.columns]
    surv_keep = surv_old[keep_cols].drop_duplicates("cohort") if keep_cols else None
except FileNotFoundError:
    surv_keep = None
print("recomputing Cox properly...", flush=True)
rows = []
for cohort in ["TCGA-SKCM", "TCGA-LUAD", "TCGA-BLCA", "GSE135222", "IMvigor210"]:
    s = pd.read_csv(f"{TAB}/scores_{cohort}.csv", index_col=0)
    s.index = s.index.map(str)
    clin = pd.read_csv(f"{PROC}/{cohort}_clin.csv", dtype=str)
    clin["sample_id"] = clin["sample_id"].astype(str)
    df = s.join(clin.set_index("sample_id"), how="inner")
    df["t"] = pd.to_numeric(df["time"], errors="coerce"); df["e"] = pd.to_numeric(df["event"], errors="coerce")
    df = df[df.t.notna() & df.e.notna() & (df.t > 0)].copy()
    df["ECO_hi"] = (df["ECO"] >= df["ECO"].median()).astype(int)
    cph = CoxPHFitter().fit(df[["t", "e", "ECO_hi"]], "t", "e")
    hr = float(cph.hazard_ratios_["ECO_hi"]); ci = np.exp(cph.confidence_intervals_.loc["ECO_hi"].values)
    p = float(cph.summary.loc["ECO_hi", "p"])
    row = {"cohort": cohort, "n": len(df), "events": int(df.e.sum()), "HR": round(hr, 3),
           "CI_low": round(float(ci[0]), 3), "CI_high": round(float(ci[1]), 3), "p": p}
    # multivariate
    if cohort.startswith("TCGA"):
        df["age"] = pd.to_numeric(df.get("age"), errors="coerce")
        df["sexM"] = (df.get("sex", "").astype(str).str.lower() == "male").astype(int)
        st = df.get("stage", "").astype(str).str.lower()
        df["late"] = st.str.contains("iii|iv").astype(int)
        dd = df.dropna(subset=["age"])[["t", "e", "ECO_hi", "age", "sexM", "late"]]
        if len(dd) > 30:
            c2 = CoxPHFitter().fit(dd, "t", "e")
            row["adj_HR"] = round(float(c2.hazard_ratios_["ECO_hi"]), 3)
            aci = np.exp(c2.confidence_intervals_.loc["ECO_hi"].values)
            row["adj_CI"] = f"{aci[0]:.2f}-{aci[1]:.2f}"
            row["adj_p"] = float(c2.summary.loc["ECO_hi", "p"])
    if cohort == "IMvigor210":
        df["tmb"] = pd.to_numeric(df.get("tmb"), errors="coerce")
        dd = df.dropna(subset=["tmb"])[["t", "e", "ECO_hi", "tmb"]]
        if len(dd) > 30:
            dd["tmb_hi"] = (dd["tmb"] >= dd["tmb"].median()).astype(int)
            c2 = CoxPHFitter().fit(dd[["t", "e", "ECO_hi", "tmb_hi"]], "t", "e")
            row["adj_HR"] = round(float(c2.hazard_ratios_["ECO_hi"]), 3)
            aci = np.exp(c2.confidence_intervals_.loc["ECO_hi"].values)
            row["adj_CI"] = f"{aci[0]:.2f}-{aci[1]:.2f}"
            row["adj_p"] = float(c2.summary.loc["ECO_hi", "p"])
    rows.append(row)
surv2 = pd.DataFrame(rows)
# keep original logrank p / endpoint if available
if surv_keep is not None and "logrank_p" in surv_keep.columns:
    surv2 = surv2.merge(surv_keep, on="cohort", how="left")
else:
    from lifelines.statistics import logrank_test as _lr
    surv2["endpoint"] = surv2["cohort"].map(lambda c: "PFS" if c == "GSE135222" else "OS")
surv2.to_csv(f"{TAB}/survival_stats.csv", index=False)
print(surv2.to_string(), flush=True)

# ---- 2. bootstrap AUC CIs + ECO-high/low OR ----
resp = pd.read_csv(f"{TAB}/response_AUCs.csv")
or_rows = []
for cohort in resp.cohort.unique():
    s = pd.read_csv(f"{TAB}/scores_{cohort}.csv", index_col=0)
    s.index = s.index.map(str)
    clin = pd.read_csv(f"{PROC}/{cohort}_clin.csv", dtype=str)
    cc = clin.set_index("sample_id")["response01"].pipe(lambda x: pd.to_numeric(x, errors="coerce")).dropna().astype(int)
    common = [c for c in s.index if str(c) in set(cc.index)]
    y = cc.loc[[str(c) for c in common]].values
    for metric in ["ECO", "TLS", "STROMA", "IFNG6", "CD8EFF", "CD274"]:
        v = s.loc[common, metric].values
        lo, hi = boot_auc(y, v)
        resp.loc[(resp.cohort == cohort) & (resp.metric == metric), "AUC_CI_low"] = round(lo, 3)
        resp.loc[(resp.cohort == cohort) & (resp.metric == metric), "AUC_CI_high"] = round(hi, 3)
    e = s.loc[common, "ECO"].values
    hi = e >= np.median(e)
    a, b = int(((y == 1) & hi).sum()), int(((y == 0) & hi).sum())
    c, d = int(((y == 1) & ~hi).sum()), int(((y == 0) & ~hi).sum())
    oddsr = (a * d) / max(b * c, 1)
    _, p = fisher_exact([[a, b], [c, d]])
    or_rows.append({"cohort": cohort, "n": len(y), "resp_rate_ECOhi": round(a / max(a + b, 1), 3),
                    "resp_rate_ECOlo": round(c / max(c + d, 1), 3), "OR": round(float(oddsr), 2), "Fisher_p": float(p)})
resp.to_csv(f"{TAB}/response_AUCs.csv", index=False)
or_df = pd.DataFrame(or_rows)
or_df.to_csv(f"{TAB}/ECO_highlow_OR.csv", index=False)
print("\n", resp[resp.metric == "ECO"].to_string(), flush=True)
print("\n", or_df.to_string(), flush=True)

# ---- 3. correct subtype labels by centroid profile ----
cent = pd.read_csv(f"{TAB}/subtype_centroids.csv", index_col=0)
def label_row(t, st):
    if t > 0.3 and st < 0.3:
        return "TLS-high"
    if st > 0.3 and t < 0.35:
        return "Stromal"
    return "Immune-desert"
new_labels = {old: label_row(t, st) for old, (t, st) in cent.iterrows()}
print("\ncentroid relabel:", new_labels, flush=True)
cent.index = [new_labels[i] for i in cent.index]
cent.to_csv(f"{TAB}/subtype_centroids.csv")
sub_rows = []
for f in os.listdir(TAB):
    if f.startswith("scores_"):
        p = f"{TAB}/{f}"
        s = pd.read_csv(p, index_col=0)
        s.index = s.index.map(str)
        s["subtype"] = s["subtype"].map(new_labels)
        s.to_csv(p)
        cohort = f.replace("scores_", "").replace(".csv", "")
        clin = pd.read_csv(f"{PROC}/{cohort}_clin.csv", dtype=str)
        if "response01" in clin.columns:
            cc = clin.set_index("sample_id")["response01"].pipe(lambda x: pd.to_numeric(x, errors="coerce")).dropna().astype(int)
            for st in ["TLS-high", "Stromal", "Immune-desert"]:
                idx = s[s.subtype == st].index
                yy = cc.loc[[str(c) for c in idx if str(c) in set(cc.index)]]
                if len(yy) > 0:
                    sub_rows.append({"cohort": cohort, "subtype": st, "n": len(yy),
                                     "resp_rate": round(float(yy.mean()), 3)})
pd.DataFrame(sub_rows).to_csv(f"{TAB}/subtype_response.csv", index=False)
print(pd.DataFrame(sub_rows).to_string(), flush=True)

# ---- 4. ECO+TMB combo in IMvigor210 ----
s = pd.read_csv(f"{TAB}/scores_IMvigor210.csv", index_col=0)
s.index = s.index.map(str)
clin = pd.read_csv(f"{PROC}/IMvigor210_clin.csv", dtype=str)
cc = clin.set_index("sample_id")
common = [c for c in s.index if str(c) in set(cc.index)]
dd = pd.DataFrame({"ECO": s.loc[common, "ECO"].values,
                   "y": pd.to_numeric(cc.loc[[str(c) for c in common], "response01"], errors="coerce").values,
                   "tmb": pd.to_numeric(cc.loc[[str(c) for c in common], "tmb"], errors="coerce").values})
dd = dd.dropna()
from sklearn.linear_model import LogisticRegression
X = np.column_stack([(dd.ECO - dd.ECO.mean()) / dd.ECO.std(), (dd.tmb - dd.tmb.mean()) / dd.tmb.std()])
for name, col in [("ECO", X[:, [0]]), ("TMB", X[:, [1]]), ("ECO+TMB", X)]:
    lr = LogisticRegression().fit(col, dd.y)
    auc = roc_auc_score(dd.y, lr.predict_proba(col)[:, 1])
    lo, hi = boot_auc(dd.y, lr.predict_proba(col)[:, 1])
    print(f"IMvigor210 {name}: AUC={auc:.3f} [{lo:.3f}-{hi:.3f}] n={len(dd)}", flush=True)
dd.to_csv(f"{TAB}/IMvigor210_ECO_TMB.csv", index=False)
print("DONE refine", flush=True)
