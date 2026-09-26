"""Enhancement analysis + Figures 7-9: pathways, networks, regulation, drugs, TIDE, algorithm."""
import os, warnings
import numpy as np
import pandas as pd
import requests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import networkx as nx
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.linear_model import LogisticRegression
warnings.filterwarnings("ignore")

BASE = "/home/user/npj_immunotherapy_paper"
TAB, PROC, FIG, ENH = f"{BASE}/results/tables", f"{BASE}/data/processed", f"{BASE}/results/figures", f"{BASE}/results/enhancement"
plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold",
                     "figure.dpi": 150, "savefig.dpi": 300, "font.family": "sans-serif"})
TEAL, RUST, GREY = "#1B9E77", "#E15728", "#7570B3"
import sys
sys.path.insert(0, f"{BASE}/code")
from signatures import TLS, STROMA

# ---------- A. PATIENT-level benefit endpoints ----------
BS = "https://www.cbioportal.org/api"
SID = "blca_iatlas_imvigor210_2017"
for a in ["CLINICAL_BENEFIT", "PROGRESSION"]:
    r = requests.get(f"{BS}/studies/{SID}/clinical-data",
                     params={"attributeId": a, "clinicalDataType": "PATIENT", "projection": "SUMMARY", "pageSize": 1000}, timeout=120)
    d = {x.get("patientId"): x["value"] for x in r.json()}
    pd.DataFrame(list(d.items()), columns=["sample_id", a]).to_csv(f"{ENH}/imv_{a}.csv", index=False)
    print(a, "n=", len(d), pd.Series(list(d.values())).value_counts().to_dict(), flush=True)

# ---------- B. TIDE vs ECO ----------
s = pd.read_csv(f"{TAB}/scores_IMvigor210.csv", index_col=0)
s.index = s.index.map(str)
clin = pd.read_csv(f"{PROC}/IMvigor210_clin.csv", dtype=str).set_index("sample_id")
tide = pd.read_csv(f"{ENH}/imv_TIDE_RESPONDER.csv", dtype=str).set_index("sample_id")["TIDE_RESPONDER"]
common = [c for c in s.index if c in set(clin.index) and c in set(tide.index)]
y = pd.to_numeric(clin.loc[common, "response01"], errors="coerce")
tpred = (tide.loc[common].str.upper() == "TRUE").astype(int)
eco = s.loc[common, "ECO"]
dd = pd.DataFrame({"y": y, "tide": tpred, "eco": eco}).dropna()
dd.y = dd.y.astype(int)
tide_auc = (float(((dd.tide == 1) & (dd.y == 1)).sum() / (dd.y == 1).sum()) +
            float(((dd.tide == 0) & (dd.y == 0)).sum() / (dd.y == 0).sum())) / 2
eco_auc = roc_auc_score(dd.y, dd.eco)
lr = LogisticRegression().fit(dd[["eco", "tide"]], dd.y)
combo_auc = roc_auc_score(dd.y, lr.predict_proba(dd[["eco", "tide"]])[:, 1])
print(f"TIDE vs ECO (n={len(dd)}): TIDE_AUC={tide_auc:.3f} ECO_AUC={eco_auc:.3f} COMBO={combo_auc:.3f}", flush=True)
pd.DataFrame([{"n": len(dd), "TIDE_AUC": round(tide_auc, 3), "ECO_AUC": round(eco_auc, 3),
               "COMBO_AUC": round(combo_auc, 3)}]).to_csv(f"{ENH}/tide_comparison.csv", index=False)

# ECO x TMB quadrants
tmb = pd.read_csv(f"{TAB}/IMvigor210_ECO_TMB.csv")
tmb["ECO_hi"] = (tmb.ECO >= tmb.ECO.median()).astype(int)
tmb["TMB_hi"] = (tmb.tmb >= tmb.tmb.median()).astype(int)
quad = tmb.groupby(["ECO_hi", "TMB_hi"]).agg(n=("y", "size"), rate=("y", "mean")).reset_index()
quad.to_csv(f"{ENH}/eco_tmb_quadrants.csv", index=False)
print(quad.to_string(index=False), flush=True)

# ---------- C. co-expression network (TCGA pooled) ----------
frames = []
for c in ["TCGA-SKCM", "TCGA-LUAD", "TCGA-BLCA"]:
    e = pd.read_csv(f"{PROC}/{c}_expr.csv", index_col=0)
    z = e.sub(e.mean(axis=1), axis=0).div(e.std(axis=1).replace(0, np.nan), axis=0)
    frames.append(z)
pool = pd.concat(frames, axis=1)
print("pooled TCGA:", pool.shape, flush=True)
corr = pool.T.corr(method="spearman")
corr.index.name = None
corr.columns.name = None
corr.to_csv(f"{ENH}/coexpression_spearman.csv")
# top edges
mask = np.triu(np.ones(corr.shape, bool), 1)
pairs = corr.where(mask).stack().reset_index()
pairs.columns = ["A", "B", "rho"]
pairs = pairs.reindex(pairs.rho.abs().sort_values(ascending=False).index)
strong = pairs[pairs.rho.abs() >= 0.55]
strong.to_csv(f"{ENH}/coexpression_strong_edges.csv", index=False)
print(f"strong coexpression edges |rho|>=0.55: {len(strong)}", flush=True)
G = nx.Graph()
for _, r in strong.iterrows():
    G.add_edge(r["A"], r["B"], weight=abs(r["rho"]))
deg = pd.Series(dict(G.degree())).sort_values(ascending=False)
deg.to_csv(f"{ENH}/coexpression_hubs.csv", header=["degree"])
print("top coexpression hubs:\n", deg.head(12).to_string(), flush=True)

# STRING hubs
st = pd.read_csv(f"{ENH}/string_edges.csv")
Gs = nx.Graph()
for _, r in st.iterrows():
    Gs.add_edge(r["A"], r["B"], weight=r["score"])
sdeg = pd.Series(dict(Gs.degree())).sort_values(ascending=False)
sdeg.to_csv(f"{ENH}/string_hubs.csv", header=["degree"])
print("top STRING hubs:\n", sdeg.head(12).to_string(), flush=True)

def node_color(g):
    if g in TLS and g in STROMA:
        return "#B7791F"
    if g in TLS:
        return TEAL
    if g in STROMA:
        return RUST
    return "#8DA0B4"

# ---------- FIG 7: pathways + TF + network ----------
fig, axes = plt.subplots(2, 2, figsize=(12, 9))
# A. KEGG
kegg_t = pd.read_csv(f"{ENH}/enrich_TLS_KEGG.csv").head(8)
kegg_s = pd.read_csv(f"{ENH}/enrich_STROMA_KEGG.csv").head(8)
ax = axes[0, 0]
terms = list(kegg_t.term.str.replace("Homo sapiens.*", "", regex=True).str[:42]) + \
        list(kegg_s.term.str.replace("Homo sapiens.*", "", regex=True).str[:42])
vals = list(-np.log10(kegg_t.adj_p + 1e-300)) + list(-np.log10(kegg_s.adj_p + 1e-300))
cols = [TEAL] * len(kegg_t) + [RUST] * len(kegg_s)
y = np.arange(len(terms))
ax.barh(y, vals, color=cols)
ax.set_yticks(y); ax.set_yticklabels(terms, fontsize=7)
ax.set_xlabel("-log10 adj.p (KEGG 2021)"); ax.set_title("A. Pathway programs: TLS (teal) vs Stroma (rust)")
ax.invert_yaxis()
# B. TFs
ax = axes[0, 1]
tf_t = pd.read_csv(f"{ENH}/enrich_TLS_TRRUST.csv").head(6)
tf_s = pd.read_csv(f"{ENH}/enrich_STROMA_TRRUST.csv").head(6)
terms = list(tf_t.term.str[:30]) + list(tf_s.term.str[:30])
vals = list(-np.log10(tf_t.adj_p + 1e-300)) + list(-np.log10(tf_s.adj_p + 1e-300))
cols = [TEAL] * len(tf_t) + [RUST] * len(tf_s)
y = np.arange(len(terms))
ax.barh(y, vals, color=cols)
ax.set_yticks(y); ax.set_yticklabels(terms, fontsize=7)
ax.set_xlabel("-log10 adj.p (TRRUST TFs)"); ax.set_title("B. Upstream regulators differ by program")
ax.invert_yaxis()
# C. STRING network (top edges)
ax = axes[1, 0]
top = st.nlargest(90, "score")
Gn = nx.Graph()
for _, r in top.iterrows():
    Gn.add_edge(r["A"], r["B"])
pos = nx.spring_layout(Gn, seed=7, k=0.6)
degn = dict(Gn.degree())
nx.draw_networkx_nodes(Gn, pos, ax=ax, node_size=[degn[n] * 40 + 30 for n in Gn.nodes()],
                       node_color=[node_color(n) for n in Gn.nodes()], alpha=0.9)
nx.draw_networkx_edges(Gn, pos, ax=ax, alpha=0.25, width=0.7)
nx.draw_networkx_labels(Gn, pos, ax=ax, font_size=6)
ax.set_title("C. Protein interactions (STRING, score≥0.7, top 90 edges)")
ax.axis("off")
# D. coexpression network
ax = axes[1, 1]
topc = strong.nlargest(90, "rho")
Gc = nx.Graph()
for _, r in topc.iterrows():
    Gc.add_edge(r["A"], r["B"])
pos = nx.spring_layout(Gc, seed=7, k=0.6)
degc = dict(Gc.degree())
nx.draw_networkx_nodes(Gc, pos, ax=ax, node_size=[degc[n] * 40 + 30 for n in Gc.nodes()],
                       node_color=[node_color(n) for n in Gc.nodes()], alpha=0.9)
nx.draw_networkx_edges(Gc, pos, ax=ax, alpha=0.25, width=0.7)
nx.draw_networkx_labels(Gc, pos, ax=ax, font_size=6)
ax.set_title("D. Co-expression in TCGA (Spearman, top 90 edges)")
ax.axis("off")
fig.suptitle("Mechanistic wiring of the TLS-versus-stroma balance", fontweight="bold", fontsize=12)
fig.tight_layout(); fig.savefig(f"{FIG}/Fig7_mechanism.png"); fig.savefig(f"{FIG}/Fig7_mechanism.pdf")
print("Fig7 done", flush=True)

# ---------- FIG 8: regulation + drugs ----------
fig, axes = plt.subplots(2, 1, figsize=(12, 10))
ax = axes[0]
ax.set_xlim(0, 12); ax.set_ylim(0, 6); ax.axis("off")
ax.set_title("A. Regulation logic: what turns each program up or down", fontsize=11)
# TLS column
ax.add_patch(patches.FancyBboxPatch((0.3, 0.5), 5.2, 4.6, boxstyle="round,pad=0.05", fc="#D5F0E4", ec=TEAL, lw=2))
ax.text(2.9, 4.7, "TLS PROGRAM → response", ha="center", fontweight="bold", color="#0B6E4F")
ax.text(2.9, 4.1, "UP: NF-κB / RELA · IRF1 · IKBKB\n(interferon sensing, B/Tfh recruitment)", ha="center", fontsize=8,
        bbox=dict(fc="white", ec="0.7"))
ax.text(2.9, 3.1, "Genes: CCL19/CCL21/CCR7 (homing)\nCXCL13/CXCR5 (B follicles) · MS4A1/CD79A (B)\nLTB · CCL2–5/CXCL9–11 (chemotaxis)", ha="center", fontsize=8)
ax.text(2.9, 1.9, "EFFECT: lymphoid organization →\nlocal priming → ICI benefit", ha="center", fontsize=8, style="italic")
ax.text(2.9, 1.0, "ECO = TLS − STROMA", ha="center", fontweight="bold", fontsize=10,
        bbox=dict(fc="#FFF2CC", ec="0.5"))
# Stroma column
ax.add_patch(patches.FancyBboxPatch((6.5, 0.5), 5.2, 4.6, boxstyle="round,pad=0.05", fc="#F9DECF", ec=RUST, lw=2))
ax.text(9.1, 4.7, "STROMAL PROGRAM → resistance", ha="center", fontweight="bold", color="#8C2F0E")
ax.text(9.1, 4.1, "UP: TWIST2 · ETV4 · SP1 · TGF-β/SMAD\nDOWN: miR-29a/b/c (anti-fibrotic)", ha="center", fontsize=8,
        bbox=dict(fc="white", ec="0.7"))
ax.text(9.1, 3.1, "Genes: FAP/ACTA2/PDPN (CAF)\nCOL1/3/5/6 · POSTN/FN1/SPARC (matrix)\nTGFB1–3 · MMP2/11/14 · LOX (remodeling)", ha="center", fontsize=8)
ax.text(9.1, 1.9, "EFFECT: T-cell exclusion →\ncheckpoint resistance", ha="center", fontsize=8, style="italic")
ax.text(9.1, 1.0, "miR-29 loss ⇒ stromal HIGH", ha="center", fontweight="bold", fontsize=9,
        bbox=dict(fc="white", ec=RUST))
ax.annotate("", xy=(6.5, 2.8), xytext=(5.5, 2.8), arrowprops=dict(arrowstyle="<->", lw=1.5))
ax.text(6.0, 3.0, "balance", ha="center", fontsize=8, style="italic")

# B. drug-gene network
ax = axes[1]
dg = pd.read_csv(f"{ENH}/dgidb_drugs.csv")
# pick per gene: approved first (max 3), then investigational (max 2)
pick = []
for g, d in dg.groupby("gene"):
    ap = d[d.approved == True].head(3)
    iv = d[d.approved != True].head(2)
    pick.append(pd.concat([ap, iv]))
pick = pd.concat(pick)
Gb = nx.Graph()
for _, r in pick.iterrows():
    Gb.add_edge(r["gene"], r["drug"].title()[:26], approved=(r["approved"] == True))
pos = nx.spring_layout(Gb, seed=11, k=0.9)
genes = [n for n in Gb.nodes() if n in set(pick.gene)]
drugs = [n for n in Gb.nodes() if n not in set(pick.gene)]
nx.draw_networkx_nodes(Gb, pos, ax=ax, nodelist=genes, node_size=500, node_color=RUST, alpha=0.9)
nx.draw_networkx_nodes(Gb, pos, ax=ax, nodelist=drugs, node_size=120, node_color=TEAL, alpha=0.8)
nx.draw_networkx_edges(Gb, pos, ax=ax, alpha=0.25, width=0.7)
nx.draw_networkx_labels(Gb, pos, labels={n: n for n in genes}, ax=ax, font_size=8, font_weight="bold")
ax.set_title(f"B. Druggable hubs: 12 ECO-network genes ↔ {len(drugs)} drugs (DGIdb; large=rust hubs, small=teal drugs)")
ax.axis("off")
fig.suptitle("Regulation and druggability of the ecosystem programs", fontweight="bold", fontsize=12)
fig.tight_layout(); fig.savefig(f"{FIG}/Fig8_regulation_drugs.png"); fig.savefig(f"{FIG}/Fig8_regulation_drugs.pdf")
print("Fig8 done", flush=True)

# ---------- FIG 9: TIDE + quadrants + algorithm ----------
fig = plt.figure(figsize=(12, 8))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.1])
ax = fig.add_subplot(gs[0, 0])
fpr, tpr, _ = roc_curve(dd.y, dd.eco)
ax.plot(fpr, tpr, color=TEAL, lw=2.5, label=f"ECO {eco_auc:.2f}")
tpr_t = float(((dd.tide == 1) & (dd.y == 1)).sum() / (dd.y == 1).sum())
fpr_t = float(((dd.tide == 1) & (dd.y == 0)).sum() / (dd.y == 0).sum())
ax.scatter([fpr_t], [tpr_t], s=120, color=GREY, zorder=5, label=f"TIDE point (AUC≈{tide_auc:.2f})")
ax.plot([0, 1], [0, 1], "k:", lw=1)
ax.set_xlabel("FPR"); ax.set_ylabel("TPR")
ax.set_title(f"A. ECO vs TIDE (IMvigor210, n={len(dd)})"); ax.legend(fontsize=8)
ax = fig.add_subplot(gs[0, 1])
labels = ["ECO-lo\nTMB-lo", "ECO-lo\nTMB-hi", "ECO-hi\nTMB-lo", "ECO-hi\nTMB-hi"]
vals = [quad[(quad.ECO_hi == i) & (quad.TMB_hi == j)]["rate"].values[0] for i, j in [(0, 0), (0, 1), (1, 0), (1, 1)]]
ns = [quad[(quad.ECO_hi == i) & (quad.TMB_hi == j)]["n"].values[0] for i, j in [(0, 0), (0, 1), (1, 0), (1, 1)]]
bars = ax.bar(labels, vals, color=[GREY, "#9ECAE1", "#FDBB84", TEAL])
ax.set_ylim(0, max(vals) * 1.45); ax.set_ylabel("Response rate")
ax.set_title("B. ECO × TMB quadrants (IMvigor210)")
for b, v, n in zip(bars, vals, ns):
    ax.text(b.get_x() + b.get_width() / 2, v + 0.015, f"{v:.0%}\n(n={n})", ha="center", fontsize=8)
ax = fig.add_subplot(gs[1, :])
ax.set_xlim(0, 12); ax.set_ylim(0, 4); ax.axis("off")
ax.set_title("C. Proposed ECO-guided treatment algorithm (hypothesis for prospective testing)", fontsize=11)
ax.add_patch(patches.FancyBboxPatch((4.2, 2.9), 3.6, 0.7, boxstyle="round,pad=0.05", fc="#DCEAF5", ec="0.4", lw=1.5))
ax.text(6, 3.25, "Bulk RNA → compute ECO", ha="center", fontweight="bold")
# left arm
ax.add_patch(patches.FancyBboxPatch((0.3, 1.5), 5.2, 1.0, boxstyle="round,pad=0.05", fc="#D5F0E4", ec=TEAL, lw=2))
ax.text(2.9, 2.0, "ECO-HIGH / TLS-high → ICI monotherapy\n(expected response ≈28%; inflamed, stroma-low)", ha="center", fontsize=8)
# right arm
ax.add_patch(patches.FancyBboxPatch((6.5, 1.5), 5.2, 1.0, boxstyle="round,pad=0.05", fc="#F9DECF", ec=RUST, lw=2))
ax.text(9.1, 2.0, "ECO-LOW / Stromal → ICI + stromal combo trial\n(TGF-β / FAP / CXCR4 / VEGF / CD73 axes)", ha="center", fontsize=8)
ax.add_patch(patches.FancyBboxPatch((2.2, 0.2), 7.6, 0.7, boxstyle="round,pad=0.05", fc="#FFF2CC", ec="0.5"))
ax.text(6, 0.55, "Re-biopsy on progression → re-compute ECO (adaptive / longitudinal monitoring)", ha="center", fontsize=8, style="italic")
ax.annotate("", xy=(2.9, 2.5), xytext=(5.2, 2.9), arrowprops=dict(arrowstyle="->", lw=1.5))
ax.annotate("", xy=(9.1, 2.5), xytext=(6.8, 2.9), arrowprops=dict(arrowstyle="->", lw=1.5))
ax.annotate("", xy=(6, 0.9), xytext=(6, 1.5), arrowprops=dict(arrowstyle="->", lw=1.2))
fig.suptitle("Clinical translation: benchmarking, stratification and actionability", fontweight="bold", fontsize=12)
fig.tight_layout(); fig.savefig(f"{FIG}/Fig9_translation.png"); fig.savefig(f"{FIG}/Fig9_translation.pdf")
print("Fig9 done", flush=True)
print("ALL ENHANCEMENT DONE", flush=True)
