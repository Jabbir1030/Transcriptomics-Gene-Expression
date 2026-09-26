"""Publication figures (300 DPI PNG + PDF). All from real computed results."""
import os, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.gridspec import GridSpec
from sklearn.metrics import roc_curve, roc_auc_score
from scipy.stats import spearmanr
from lifelines import KaplanMeierFitter
warnings.filterwarnings("ignore")

TAB = "/home/user/npj_immunotherapy_paper/results/tables"
PROC = "/home/user/npj_immunotherapy_paper/data/processed"
FIG = "/home/user/npj_immunotherapy_paper/results/figures"

plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold",
                     "figure.dpi": 150, "savefig.dpi": 300, "font.family": "sans-serif"})
TEAL, ORANGE, GREY, PURPLE, GREEN = "#1B9E77", "#D95F02", "#7570B3", "#984EA3", "#66A61E"
RESP_C = {1: TEAL, 0: "#E15728"}
SUB_C = {"TLS-high": TEAL, "Stromal": "#E15728", "Immune-desert": "#8DA0B4"}

ICI = ["GSE78220", "GSE91061", "GSE126044", "GSE135222", "GSE207422", "IMvigor210"]
CANCER = {"GSE78220": "MEL", "GSE91061": "MEL", "GSE126044": "NSCLC", "GSE135222": "NSCLC",
          "GSE207422": "NSCLC", "IMvigor210": "BLCA"}

def load(cohort):
    s = pd.read_csv(f"{TAB}/scores_{cohort}.csv", index_col=0)
    s.index = s.index.map(str)
    clin = pd.read_csv(f"{PROC}/{cohort}_clin.csv", dtype=str)
    return s, clin.set_index("sample_id")

# ================= FIG 1: design schematic =================
fig, ax = plt.subplots(figsize=(10, 4.5))
ax.set_xlim(0, 10); ax.set_ylim(0, 5); ax.axis("off")
ax.set_title("TLS-versus-Stroma Ecosystem Score (ECO): study design", loc="left", fontsize=12)
boxes = [
    (0.2, 2.6, 1.9, 1.8, "Discovery\nlandscape", "TCGA\nSKCM 460\nLUAD 496\nBLCA 407"),
    (2.5, 2.6, 1.9, 1.8, "ICI validation", "MEL: 25+49\nNSCLC: 16+27+24\nBLCA: 298"),
    (4.8, 2.6, 1.9, 1.8, "Score", "TLS (24 genes)\nSTROMA (30 genes)\nECO = TLS − STROMA"),
    (7.1, 2.6, 1.9, 1.8, "Analyses", "Response AUC\nOS/PFS, HR\nSubtypes, TMB+"),
]
cols = ["#DCEAF5", "#D5F0E4", "#FDE8C8", "#EADCF5"]
for (x, y, w, h, t, b), c in zip(boxes, cols):
    ax.add_patch(patches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05", fc=c, ec="0.4"))
    ax.text(x + w / 2, y + h - 0.3, t, ha="center", va="center", fontweight="bold")
    ax.text(x + w / 2, y + h / 2 - 0.2, b, ha="center", va="center", fontsize=8)
for x0, x1 in [(2.1, 2.5), (4.4, 4.8), (6.7, 7.1)]:
    ax.annotate("", xy=(x1, 3.5), xytext=(x0, 3.5), arrowprops=dict(arrowstyle="->", lw=1.5))
ax.text(5, 0.9, "No wet-lab data generated · cohort-to-cohort validation only · all cohorts public (TCGA/Firehose, GEO, cBioPortal)",
        ha="center", fontsize=8, style="italic")
ax.text(5, 0.4, "TLS: 12-chemokine + B/Tfh imprint (Cabrita 2020; Sautès-Fridman) · STROMA/CAF: TGF-β exclusion + IPRES (Mariathasan 2018; Hugo 2016)",
        ha="center", fontsize=7, color="0.35")
fig.tight_layout(); fig.savefig(f"{FIG}/Fig1_design.png"); fig.savefig(f"{FIG}/Fig1_design.pdf")
print("Fig1 done", flush=True)

# ================= FIG 2: ECO vs response (boxplots + ROCs) =================
fig, axes = plt.subplots(2, 6, figsize=(13, 5.5), gridspec_kw={"hspace": 0.5, "wspace": 0.4})
resp_tab = pd.read_csv(f"{TAB}/response_AUCs.csv")
for j, cohort in enumerate(ICI):
    s, cc = load(cohort)
    y = pd.to_numeric(cc["response01"], errors="coerce").dropna().astype(int)
    common = [c for c in s.index if c in set(y.index)]
    v = s.loc[common, "ECO"]; y = y.loc[common]
    ax = axes[0, j]
    data = [v[y == 0].values, v[y == 1].values]
    bp = ax.boxplot(data, tick_labels=["Non-R", "R"], patch_artist=True, widths=0.55)
    bp["boxes"][0].set_facecolor("#F2B8A0"); bp["boxes"][1].set_facecolor("#A8DDBB")
    auc = resp_tab[(resp_tab.cohort == cohort) & (resp_tab.metric == "ECO")].iloc[0]
    ax.set_title(f"{cohort}\n{CANCER[cohort]} n={len(y)}", fontsize=8)
    ax.text(0.5, 0.94, f"AUC {auc['AUC']:.2f}\np={auc['MWU_p']:.3f}", transform=ax.transAxes,
            ha="center", va="top", fontsize=7, bbox=dict(fc="white", ec="0.7", pad=1.5))
    ax.tick_params(labelsize=7)
    ax2 = axes[1, j]
    fpr, tpr, _ = roc_curve(y, v)
    ax2.plot(fpr, tpr, color=TEAL, lw=2, label=f"ECO {roc_auc_score(y, v):.2f}")
    ax2.plot(fpr, tpr, color="none")
    for m, c in [("IFNG6", GREY), ("CD274", "0.6")]:
        vv = s.loc[common, m].values
        f2, t2, _ = roc_curve(y, vv)
        ax2.plot(f2, t2, color=c, lw=1, ls="--", label=f"{m} {roc_auc_score(y, vv):.2f}")
    ax2.plot([0, 1], [0, 1], "k:", lw=1)
    ax2.set_xlabel("FPR", fontsize=7); ax2.tick_params(labelsize=7)
    if j == 0:
        ax.set_ylabel("ECO score"); ax2.set_ylabel("TPR"); ax2.legend(fontsize=6)
fig.suptitle("ECO predicts checkpoint-inhibitor response across 6 independent cohorts (3 cancers)", fontweight="bold")
fig.savefig(f"{FIG}/Fig2_response.png"); fig.savefig(f"{FIG}/Fig2_response.pdf")
print("Fig2 done", flush=True)

# ================= FIG 3: benchmark forest =================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw={"width_ratios": [1.35, 1]})
pivot = resp_tab.pivot(index="cohort", columns="metric", values="AUC").reindex(ICI)
x = np.arange(len(ICI)); w = 0.12
for i, m in enumerate(["ECO", "TLS", "STROMA", "IFNG6", "CD8EFF", "CD274"]):
    ax1.bar(x + (i - 2.5) * w, pivot[m].values, w, label=m,
            color={ "ECO": TEAL, "TLS": GREEN, "STROMA": "#E15728", "IFNG6": GREY, "CD8EFF": PURPLE, "CD274": "0.6"}[m])
ax1.set_xticks(x); ax1.set_xticklabels([f"{c}\n{CANCER[c]}" for c in ICI], fontsize=7)
ax1.axhline(0.5, color="k", ls=":", lw=1); ax1.set_ylabel("Response AUC"); ax1.set_ylim(0.15, 0.95)
ax1.set_title("A. ECO vs single-axis & published benchmarks"); ax1.legend(fontsize=7, ncol=3)
eco = resp_tab[resp_tab.metric == "ECO"].set_index("cohort").reindex(ICI)
y = np.arange(len(ICI))
ax2.errorbar(eco["AUC"].values, y, xerr=[eco["AUC"].values - eco["AUC_CI_low"].values,
             eco["AUC_CI_high"].values - eco["AUC"].values], fmt="o", color=TEAL, ecolor="0.4", capsize=3)
ax2.axvline(0.5, color="k", ls=":", lw=1); ax2.set_yticks(y); ax2.set_yticklabels(ICI, fontsize=8)
ax2.set_xlabel("ECO AUC (95% bootstrap CI)"); ax2.set_title("B. ECO consistency (all AUC>0.5)")
fig.suptitle("Benchmarking: ECO integrates TLS and stromal signals", fontweight="bold")
fig.tight_layout(); fig.savefig(f"{FIG}/Fig3_benchmark.png"); fig.savefig(f"{FIG}/Fig3_benchmark.pdf")
print("Fig3 done", flush=True)

# ================= FIG 4: KM curves =================
fig, axes = plt.subplots(1, 5, figsize=(14, 3.6))
for ax, cohort in zip(axes, ["TCGA-SKCM", "TCGA-LUAD", "TCGA-BLCA", "IMvigor210", "GSE135222"]):
    s, cc = load(cohort)
    t = pd.to_numeric(cc["time"], errors="coerce"); e = pd.to_numeric(cc["event"], errors="coerce")
    common = [c for c in s.index if c in set(cc.index)]
    d = pd.DataFrame({"t": t.loc[common].values / 30.44, "e": e.loc[common].values,
                      "g": (s.loc[common, "ECO"] >= s["ECO"].median()).astype(int).values})
    d = d[d.t.notna() & d.e.notna() & (d.t > 0)]
    kmf = KaplanMeierFitter()
    for g, c, lab in [(1, TEAL, "ECO-high"), (0, "#E15728", "ECO-low")]:
        kmf.fit(d[d.g == g]["t"], d[d.g == g]["e"], label=f"{lab} (n={(d.g==g).sum()})")
        kmf.plot(ax=ax, color=c, lw=2)
    st = pd.read_csv(f"{TAB}/survival_stats.csv").set_index("cohort").loc[cohort]
    ax.set_title(f"{cohort}\nHR {st['HR']:.2f} [{st['CI_low']:.2f}-{st['CI_high']:.2f}] p={st['p']:.4f}", fontsize=8)
    ax.set_xlabel("Months"); ax.tick_params(labelsize=7)
axes[0].set_ylabel("Survival probability")
fig.suptitle("ECO-high associates with longer OS/PFS in TCGA and ICI cohorts", fontweight="bold")
fig.tight_layout(); fig.savefig(f"{FIG}/Fig4_survival.png"); fig.savefig(f"{FIG}/Fig4_survival.pdf")
print("Fig4 done", flush=True)

# ================= FIG 5: subtypes =================
fig = plt.figure(figsize=(11, 4.2))
gs = GridSpec(1, 2, width_ratios=[1, 1.2])
ax = fig.add_subplot(gs[0])
frames = []
for c in ["TCGA-SKCM", "TCGA-LUAD", "TCGA-BLCA"]:
    s, _ = load(c)
    frames.append(s[["TLS", "STROMA", "subtype"]])
tcga = pd.concat(frames)
for st, c in SUB_C.items():
    d = tcga[tcga.subtype == st]
    ax.scatter(d["TLS"], d["STROMA"], s=6, alpha=0.5, color=c, label=f"{st} (n={len(d)})")
ax.set_xlabel("TLS program (z)"); ax.set_ylabel("STROMA program (z)")
ax.set_title("A. Ecosystem subtypes (TCGA pooled, k=3)"); ax.legend(fontsize=7)
ax2 = fig.add_subplot(gs[1])
sub = pd.read_csv(f"{TAB}/subtype_response.csv")
order = ["TLS-high", "Stromal", "Immune-desert"]
piv = sub.pivot(index="cohort", columns="subtype", values="resp_rate").reindex(ICI)
x = np.arange(len(ICI)); w = 0.22
for i, st in enumerate(order):
    ax2.bar(x + (i - 1) * w, piv[st].values, w, label=st, color=SUB_C[st])
ax2.set_xticks(x); ax2.set_xticklabels([f"{c}\n{CANCER[c]}" for c in ICI], fontsize=7)
ax2.set_ylabel("Response rate"); ax2.set_ylim(0, 1)
ax2.set_title("B. Response rate by subtype (mapped to TCGA centroids)"); ax2.legend(fontsize=7)
fig.suptitle("Tumor ecosystem subtypes and immunotherapy response", fontweight="bold")
fig.tight_layout(); fig.savefig(f"{FIG}/Fig5_subtypes.png"); fig.savefig(f"{FIG}/Fig5_subtypes.pdf")
print("Fig5 done", flush=True)

# ================= FIG 6: TMB + checkpoints =================
fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))
s, cc = load("IMvigor210")
tmb = pd.to_numeric(cc["tmb"], errors="coerce")
common = [c for c in s.index if c in set(cc.index)]
dd = pd.DataFrame({"ECO": s.loc[common, "ECO"].values, "TMB": tmb.loc[common].values,
                   "y": pd.to_numeric(cc.loc[common, "response01"], errors="coerce").values}).dropna()
rho, rhop = spearmanr(dd.ECO, dd.TMB)
axes[0].scatter(dd.TMB, dd.ECO, c=[RESP_C[int(v)] for v in dd.y], s=12, alpha=0.6)
axes[0].set_xlabel("TMB (nonsynonymous)"); axes[0].set_ylabel("ECO")
axes[0].set_title(f"A. ECO vs TMB (ρ={rho:.2f}, p={rhop:.3f})")
med = dd.TMB.median()
res = []
for lab, part in [("TMB-low", dd[dd.TMB < med]), ("TMB-high", dd[dd.TMB >= med])]:
    a = roc_auc_score(part.y, part.ECO)
    res.append((lab, len(part), a))
    print(f"{lab}: n={len(part)} ECO_AUC={a:.3f}", flush=True)
pd.DataFrame(res, columns=["group", "n", "ECO_AUC"]).to_csv(f"{TAB}/TMB_stratified_ECO.csv", index=False)
axes[1].bar([r[0] for r in res], [r[2] for r in res], color=[GREY, TEAL])
axes[1].axhline(0.5, color="k", ls=":"); axes[1].set_ylim(0.3, 0.9); axes[1].set_ylabel("ECO AUC")
axes[1].set_title("B. ECO within TMB strata")
for i, (lab, n, a) in enumerate(res):
    axes[1].text(i, a + 0.02, f"{a:.2f}\nn={n}", ha="center", fontsize=8)
# checkpoints by ECO group
expr = pd.read_csv(f"{PROC}/IMvigor210_expr.csv", index_col=0)
expr.columns = expr.columns.map(str)
ck = [g for g in ["CD274", "PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "PDCD1LG2", "IFNG", "CXCL9", "CXCL13", "FAP", "TGFB1"] if g in expr.index]
z = expr.loc[ck].sub(expr.loc[ck].mean(axis=1), axis=0).div(expr.loc[ck].std(axis=1), axis=0)
hi = [c for c in common if s.loc[c, "ECO"] >= s["ECO"].median()]
lo = [c for c in common if s.loc[c, "ECO"] < s["ECO"].median()]
mat = pd.DataFrame({"ECO-low": z[lo].mean(axis=1), "ECO-high": z[hi].mean(axis=1)})
im = axes[2].imshow(mat.values, aspect="auto", cmap="RdBu_r", vmin=-0.6, vmax=0.6)
axes[2].set_yticks(range(len(ck))); axes[2].set_yticklabels(ck, fontsize=8)
axes[2].set_xticks([0, 1]); axes[2].set_xticklabels(["ECO-low", "ECO-high"])
axes[2].set_title("C. Immuno-stromal genes (IMvigor210)")
fig.colorbar(im, ax=axes[2], shrink=0.8, label="mean z")
fig.suptitle("ECO is orthogonal to TMB and marks an inflamed, stroma-low ecosystem", fontweight="bold")
fig.tight_layout(); fig.savefig(f"{FIG}/Fig6_TMB.png"); fig.savefig(f"{FIG}/Fig6_TMB.pdf")
print("Fig6 done", flush=True)
print("ALL FIGURES DONE")
