# Beginner's Guide to the ECO Study
### *A TLS-versus-stroma ecosystem score predicts checkpoint inhibitor benefit across three cancers*
**A step-by-step course: what we did, why we did it, and how each analysis works — starting from zero.**

---

## How to use this guide
- Read **Chapters 1–2** for the story and vocabulary (30 min).
- Read **Chapters 3–5** to understand the data, the score, and the pipeline (45 min).
- Read **Chapter 6** slowly — it teaches all 11 analyses one by one (the core of the course).
- Use **Chapters 7–10** as reference (figures, reviewer Q&A, rerunning code, glossary).

All file paths below are inside `/home/user/npj_immunotherapy_paper/`.

---

## Chapter 1. The story in 5 minutes

**The problem.** Cancer drugs called *checkpoint inhibitors* (anti-PD-1, anti-PD-L1) can cure some patients with melanoma, lung, or bladder cancer — but **most patients don't benefit**. Doctors can't reliably predict who will respond. Giving the drug to everyone means many suffer side effects and lose precious months for nothing.

**The clue from biology.** Scientists noticed two tissue features that predict the outcome:
1. **TLS (tertiary lymphoid structures)** — tiny "immune training camps" inside/near the tumor where B cells and T cells organize. Tumors rich in TLS tend to **respond**.
2. **Stroma** — scar-like tissue built by fibroblasts (collagen walls + TGF-beta signals) that **blocks** T cells from entering. Stroma-rich tumors tend to **resist**.

**Nobody measured both together in one simple test.** That was the gap.

**Our idea (the whole paper in one sentence):** Compute one number per tumor — *immune organization minus stromal blockage* — and check whether it predicts who benefits from immunotherapy.

**What we found:**
- In **6 immunotherapy cohorts (439 patients)**, responders scored higher every single time (AUC 0.60–0.78).
- High-score patients **lived longer** in all 3 big TCGA cohorts + the bladder trial (risk of death cut by ~30–47%).
- Stromal-type tumors responded worst (19.2%).
- The score works **independently of tumor mutational burden** (a DNA test), beats the famous TIDE predictor, and points to **drugs** for low-score patients.

**Why it matters:** the score needs only a routine tumor RNA readout, uses zero machine learning (so doctors can understand and trust it), and directly suggests which drug combination to try next for non-responders.

**Analogy:** think of a cricket match. TLS = your team's batting strength; stroma = the opponent's bowling attack plus a bad pitch. You can't predict the match by looking at only one side. ECO = *our strength minus their defense*. One number, both sides.

---

## Chapter 2. Concepts crash course (10 mini-lessons)

### 2.1 Checkpoint inhibitors (anti-PD-1 / anti-PD-L1 / anti-CTLA-4)
Cancer cells press an "off switch" (PD-L1) on passing T cells. Checkpoint drugs **block that switch**, waking T cells up to kill the tumor. Patients in this study got pembrolizumab, nivolumab, atezolizumab, or chemo + anti-PD-1.

### 2.2 Response vs. survival (two different questions)
- **Response** = did the tumor shrink right after treatment? (RECIST rules: CR/PR = responder; SD/PD = non-responder.)
- **Survival** = how long did the patient live? (OS = overall survival; PFS = time until the cancer grew again.)
- A good biomarker should predict **both**. ECO does.

### 2.3 TLS — the immune training camp
Lymph-node-like hubs with B-cell zones, T-cell zones, dendritic cells. A mature TLS = local immunity factory. Measured here via **24 TLS genes** (chemokines like CXCL13, CCL19, CCL21 + B-cell genes like MS4A1/CD20, CD79A + homing receptors CCR7, CXCR5).

### 2.4 Stroma/CAF — the wall
Cancer-associated fibroblasts (CAFs) lay collagen (COL1A1, COL3A1…), crosslink it (LOX), contract it (ACTA2), and secrete TGF-beta. Result: T cells stuck outside. Measured via **30 stromal genes** (FAP, collagens, POSTN, VIM, TGFB1/2/3, MMPs…).

### 2.5 Bulk RNA sequencing (the raw data)
A tumor sample is ground up and all its RNA is read: one number per gene = "how active is this gene?" That's a **gene expression matrix** (rows = genes, columns = patients). All our cohorts are this format.

### 2.6 Gene signature / score
Instead of staring at 20,000 genes, average a meaningful handful into one number. ECO = (average of 24 TLS genes) − (average of 30 stroma genes). That's it. No weights, no AI model.

### 2.7 AUC — the 0-to-1 report card for prediction
- AUC = probability that a random responder scores higher than a random non-responder.
- **0.5 = coin flip (useless); 1.0 = perfect.** In cancer biomarkers, 0.60–0.80 is the realistic, publishable range.
- ECO: 0.60–0.78 across cohorts. Modest per cohort, but **never below 0.5 even once** — that consistency is the story.

### 2.8 Hazard ratio (HR) — survival report card
- HR compares death rates between high-score and low-score halves. **HR = 0.53 means the high group dies at roughly half the rate at any given time.**
- Below 1.0 = good for the high group. Ours: 0.53, 0.63, 0.64, 0.69. All significant.

### 2.9 p-value — "could this be luck?"
- p < 0.05 = less than 5% chance of seeing this result if nothing were real → we call it **statistically significant**.
- We used **two-sided** tests (fair in both directions) and report exact p-values (0.049, 0.014…), as the journal demands.

### 2.10 TMB and TIDE (the rivals we compare against)
- **TMB** (tumor mutational burden) = number of DNA mutations; more mutations = more targets for T cells. A DNA test.
- **TIDE** = a famous published RNA predictor of immunotherapy response. We beat it head-to-head (0.60 vs 0.49) on its own test data.

---

## Chapter 3. The data — 9 cohorts, ~1,830 tumors

Everything is **public data** (no new patients). Three sources:

| Source | What it is | What we took |
|---|---|---|
| **TCGA** (via Firehose) | The Cancer Genome Atlas: thousands of untreated tumors with RNA + survival | SKCM melanoma 440, LUAD lung 477, BLCA bladder 394 → prognosis + subtype map |
| **GEO** | Public gene-expression archive | 5 small immunotherapy cohorts (melanoma ×2, NSCLC ×3) with response labels |
| **IMvigor210** (via cBioPortal) | Big bladder trial of atezolizumab | 298 with response, 326 with survival, TMB subset → our heavyweight validation |

**Why 9 cohorts?** One cohort can lie (luck, quirks). Nine cohorts across 3 cancers and 5 drug settings that *all agree* = a result reviewers believe. This design is called **cohort-to-cohort validation**.

**Endpoints (how "success" was defined per cohort):**
- RECIST (tumor shrinkage: CR/PR vs SD/PD) — GSE78220, GSE91061, IMvigor210.
- DCB (durable clinical benefit: no progression for 180 days) — GSE135222 (it had no RECIST labels, so we used its survival data honestly instead of dropping it).
- MPR (major pathologic response at surgery) — GSE207422 neoadjuvant cohort.
- OS/PFS (survival times) — TCGA + IMvigor210 + GSE135222.

**Benefit of this step:** zero lab cost, fully reproducible by anyone, and mixing big + small, old + new cohorts stress-tests the score.

⚠️ **Data honesty example:** we inspected GSE176307 and **threw it out** because its response labels were shifted across samples — and we say so in the paper to warn future users. Reviewers love this kind of candor.

---

## Chapter 4. The ECO score — built step by step

### 4.1 The formula (memorize this)
> **ECO = (mean of 24 TLS genes) − (mean of 30 stroma genes)**, each gene first standardized to a **z-score** within its cohort.

### 4.2 Why z-scores? (fair comparison)
Gene A ranges 0–100, gene B ranges 0–5. Averaging raw values would let gene A dominate. A **z-score** = (value − mean) / standard deviation → every gene gets mean 0, spread 1. Now all genes vote equally. Standardizing **within each cohort** also absorbs lab-to-lab technical differences.

### 4.3 Toy example (fake numbers, just to learn the math)
Tumor X, 3 TLS genes with z-scores [+1.0, +0.5, −0.2] → TLS mean = +0.43. Three stroma genes [+0.1, +0.9, +0.5] → stroma mean = +0.50. **ECO = 0.43 − 0.50 = −0.07** → slightly stromal → predicted poor responder. Real ECO uses 24 and 30 genes the same way.

### 4.4 Why the lists were "frozen" BEFORE seeing outcomes (the most important methods concept)
Most signature papers cheat by accident: they hunt through genes until something fits their data (*overfitting*), then it fails on new patients. We did the opposite: copied both gene lists from **published biology papers** (Messina 12-chemokine signal + TLS classics for TLS; Mariathasan/Hugo stromal programs for stroma), locked them in `code/signatures.py`, and only then looked at outcomes. So all 9 cohorts are **genuine validation**. This single choice is why the paper is credible.

### 4.5 The gene lists (for reference)
- **TLS (24):** CCL2, CCL3, CCL4, CCL5, CCL8, CCL18, CCL19, CCL21, CXCL9, CXCL10, CXCL11, CXCL13, MS4A1, CD79A, CD79B, LTB, CCR7, CCR6, CXCR5, SELL, CD1D, LAT, SKAP1, CETP
- **Stroma (30):** FAP, ACTA2, COL1A1, COL1A2, COL3A1, COL5A1, COL5A2, COL6A1, COL6A2, COL6A3, POSTN, VIM, PDPN, PDGFRA, PDGFRB, TGFB1, TGFB2, TGFB3, TGFBR2, MMP2, MMP11, MMP14, LOX, LOXL2, SPARC, FN1, VCAN, THY1, ITGB1, ZEB1

---

## Chapter 5. The pipeline — what each script does and why

| Script | Job | Why it exists (benefit) |
|---|---|---|
| `01_download.py` | Downloads TCGA (via Xena) + GEO cohorts | Real data only; prints inventory so nothing is silently missing |
| `02_download2.py` | Backup route: GEO tables + TCGA via cBioPortal | Public mirrors fail; a second route guarantees completeness |
| `03_inspect.py` | Prints shapes/headers of every file | Catch corrupted or mislabeled files BEFORE analysis (this caught the GSE176307 problem) |
| `04_parse.py`, `05_parse.py` | Convert chaos → uniform `{cohort}_expr.csv` + `{cohort}_clin.csv` (genes × samples, log2) | Every cohort arrives in a different format; standardization makes one analysis script work for all |
| `06_imvigor.py` | Fetches IMvigor210 RNA + clinical via cBioPortal API | The biggest validation cohort lives behind an API, not a download link |
| `signatures.py` | Stores the frozen 24 + 30 gene lists | Separation of biology (fixed) from statistics (computed) = no overfitting, easy audit |
| `07_scores.py` | **Core:** computes ECO + benchmark scores, response AUCs, survival, subtypes | The heart of the paper — one script, all main results |
| `07b_refine.py` | Bootstrap CIs, odds ratios, multivariate Cox, ECO+TMB combo | Turns point estimates into honest uncertainty + adjusted effects reviewers demand |
| `08_figures.py` | Draws Figures 1–6 (300 DPI PNG + PDF) | Publication-grade, reproducible graphics (no hand-editing in PowerPoint!) |
| `09_manuscript.py` | First full manuscript draft from real results | Numbers flow from code → text automatically (no copy-paste errors) |
| `10_enhance_fetch.py` | Queries Enrichr, STRING, DGIdb, TIDE data | Mechanism + drugs + rival comparison need outside databases |
| `11_enhance_analysis.py` | Analyzes them; draws Figures 7–9 | Pathways, networks, regulation, drugs, TIDE contest, treatment sketch |
| `12_update_manuscript.py` | Patches manuscript with enhancement sections | Keeps one coherent paper instead of two drafts |
| `13_build_docx.py` | Builds final `manuscript.md` + `manuscript.docx` (tables + embedded figures) | The submittable file; also verifies reference counts |
| `14_audit.py` | Style audit: sentence variety, AI-tic phrases, em-dashes | Keeps prose human and journal-ready |

**Takeaway:** data flows left → right, each step's output is the next step's input, and everything is seeded/pinned so re-running gives identical numbers. That is called **reproducibility**, and journals now require it.

---

## Chapter 6. The 11 analyses, one by one ⭐ (the core course)

> Format per analysis: **Question → Method (plain words) → Result → Benefit (why it was worth doing) → How to read the figure.**

### Analysis 1 — Study design map (Fig. 1, Table 1)
- **Question:** What data do we have and what's the plan?
- **Method:** Flowchart of cohorts → scores → validation → mechanism → clinic, plus a cohort-summary table.
- **Result:** 9 cohorts, 1,313 TCGA + 439 treated tumors with response labels.
- **Benefit:** A reviewer understands the entire paper in 60 seconds. Table 1 lets anyone verify sample sizes and endpoints.
- **Read it:** follow the arrows; check Table 1's n's match every later analysis (they do).

### Analysis 2 — Does ECO predict response? (Fig. 2, Table 2, Supp. Table S4)
- **Question:** Do responders score higher than non-responders?
- **Method:** For each cohort: split patients by response label, compare ECO with a **Mann–Whitney test** (a fair test needing no normality assumption — important because cohorts are small), draw the **ROC curve**, report **AUC + 95% bootstrap CI** (resample 2,000× to measure uncertainty). Then **pool** all 439 (standardize ECO per cohort first so cohorts are comparable) and repeat. Also split each cohort at its **median** into ECO-high/low and compute response **odds ratios** (Fisher's test).
- **Result:** Responders higher in **all 6** cohorts; AUCs 0.71/0.60/0.78/0.64/0.75/0.60; pooled 0.62; ORs 1.36–4.50 in 5/6.
- **Benefit:** This is the paper's headline claim, tested 6 independent times with zero reversals.
- **Read it:** Fig. 2 top = dot plots (responders sit higher); bottom = ROC curves bending toward the top-left = good.

### Analysis 3 — Benchmarking: why subtract? (Fig. 3, Table 2)
- **Question:** Is the *balance* better than each side alone or than famous signatures?
- **Method:** Same AUC machinery for 6 metrics: ECO, TLS-only, stroma-only, IFN-gamma 6-gene, CD8-effector, PD-L1 (CD274).
- **Result:** ECO beats TLS-only in 4/6 and PD-L1 in 6/6; ties the IFN/CD8 standards. The killer detail: in GSE78220, TLS-only AUC = 0.45 (useless) while stroma AUC = 0.21 (**inverted** — low stroma predicts response, p = 0.016) — the known IPRES resistance program, rediscovered blindly.
- **Benefit:** Proves subtraction isn't decoration: single-sided tests fail somewhere; the balance never does.
- **Read it:** Fig. 3A = heatmap of AUCs (ECO column stays warm everywhere); 3B = ECO's CIs all clear 0.5.

### Analysis 4 — Does ECO predict survival? (Fig. 4, Table 3)
- **Question:** Do high-ECO patients live longer?
- **Method:** Split each cohort at median ECO → **Kaplan–Meier curves** (step-down survival plots) + **log-rank test** + **Cox model** → **hazard ratio**. Then **multivariate Cox** adding age/sex/stage (TCGA) or TMB (IMvigor210) to check ECO isn't just a proxy for "young, early-stage patients."
- **Result:** HRs 0.53/0.63/0.64 (TCGA, all p < 0.05) and 0.69 IMvigor210 (p = 0.006); PFS HR 0.53 (trend, p = 0.15, tiny n = 27). Adjusted: SKCM 0.60 (p = 0.002), LUAD 0.69 (p = 0.052).
- **Benefit:** Response can be a fluke of measurement; survival is what patients feel. Passing both = a real biomarker.
- **Read it:** curves that separate early and stay apart; HR < 1 with CI excluding 1.0 = significant.

### Analysis 5 — Three ecosystem states (Fig. 5, Supp. Table S5)
- **Question:** Do tumors fall into natural TLS/stroma types with different response rates?
- **Method:** **k-means clustering (k = 3)** on the 1,313 TCGA tumors using only their (TLS, stroma) coordinates → name the clusters by their profiles → assign each immunotherapy tumor to its **nearest fixed centroid** (immunotherapy data never touches the fitting — no peeking).
- **Result:** TLS-high 28.1% (38/135), desert 28.2% (49/174), **stromal 19.2% (25/130)** — stromal trails pooled and in 4/6 cohorts.
- **Benefit:** Turns a continuous score into **doctor-friendly categories** and nails the biology: deserts (quiet but unguarded) respond like TLS-high; stromal walls are the true enemy.
- **Read it:** 5A = scatter with 3 clouds; 5B = bar chart, stromal bar shortest.

### Analysis 6 — ECO vs. TMB + inflamed markers (Fig. 6, Supp. Table S6)
- **Question:** Does ECO just duplicate the DNA test (TMB)? What does an ECO-high tumor look like?
- **Method:** **Spearman correlation** (rank-based, outlier-proof) between ECO and TMB; **logistic regression** combining both; per-gene expression split by ECO group.
- **Result:** Correlation only 0.26 (different bottlenecks!). TMB alone wins in bladder (0.74 vs 0.58; combo 0.72 — honestly reported). ECO-high = checkpoints/effectors up (CD274, PDCD1, CTLA4, IFNG, CXCL9, CXCL13), stromal drivers down (FAP, TGFB1).
- **Benefit:** Positions ECO correctly (a **TMB substitute** when DNA is unavailable, not a booster) — honesty that builds reviewer trust — and validates the "inflamed" interpretation.
- **Read it:** 6A = weak-scatter cloud; 6B = ROC bars; 6C = red/blue gene panel.

### Analysis 7 — Pathways & regulators (Fig. 7A–B, Supp. Table S7)
- **Question:** What biology do the two gene lists actually capture?
- **Method:** **Gene-set enrichment** (Enrichr: KEGG, Reactome) asks "which known pathways are over-represented?" **TRRUST** asks "which transcription factors control these genes?" All with **multiple-testing correction** (adjusted p — stricter because we test many pathways at once).
- **Result:** TLS = chemokine signaling (7.9×10⁻²⁴!), cytokine crosstalk, Toll-like alarms, run by NF-kB/RELA + IRF1/IRF3. Stroma = matrix organization (4.1×10⁻²⁸!), focal adhesion, collagen, run by EMT drivers TWIST2/ETV4/SP1/KLF8/SRF.
- **Benefit:** Explains *why* ECO generalizes across cancers: both control circuits are ancient, tissue-agnostic programs, not melanoma- or lung-specific quirks.
- **Read it:** 7A = bar chart of −log10(p) (longer = stronger); 7B = regulator bars split by program.

### Analysis 8 — Networks & hubs (Fig. 7C–D, Supp. Table S8)
- **Question:** How do these genes wire together, and which are the master connectors?
- **Method:** Two independent wirings: **STRING** (published protein–protein links, confidence ≥ 0.7 → 517 edges) and **co-expression** (Spearman across 1,363 TCGA tumors). Rank genes by **degree** (number of connections) → hubs.
- **Result:** Both maps show a TLS/chemokine module + a stromal/collagen module bridged by checkpoints. Hubs: CD8A, IFNG, IL10, CXCL12/CXCR4/CCR7/CXCL9, TGFB1, PDCD1, TIGIT, HAVCR2.
- **Benefit:** Agreement of two independent networks = the structure is real; hubs become the **drug-target shortlist** for Analysis 10.
- **Read it:** node size = connectedness; teal vs rust modules; bridge genes in the middle.

### Analysis 9 — The miR-29 brake (Fig. 8A, Supp. Table S7)
- **Question:** Is there a master switch controlling the stromal program?
- **Method:** **miRTarBase** enrichment: which microRNAs have validated targets over-represented in our lists?
- **Result:** miR-29 family (miR-29b p = 2.9×10⁻¹⁸!) directly suppresses **15 of 30** stromal genes (TGFB1/2/3, collagens, SPARC, LOX, MMP2). Nothing comparable on the TLS side.
- **Benefit:** Converts a correlation into a **testable therapy idea**: stromal-high ≈ miR-29-deficient → miR-29 mimics as a sensitizer. This is the paper's most novel mechanistic gem.
- **Read it:** 8A = cartoon of switches: IFN/NF-kB → TLS on; EMT/TGF-beta → stroma on; miR-29 ⊣ stroma off.

### Analysis 10 — Drug mapping (Fig. 8B, Table 4, Supp. Table S9)
- **Question:** Which existing drugs hit our hub genes?
- **Method:** Query 12 hubs against **DGIdb v4** (Drug–Gene Interaction Database) → 333 compound links; keep agents with direct mechanistic fit; sketch trial concepts per stromal flavor.
- **Result:** Approved options (pirfenidone, bevacizumab, CXCR4 blockers plerixafor/motixafortide/mavorixafor…) + experimental (FAP, CD39/CD73, IDO1 with an explicit phase-III-failure warning, MMP14). Table 4 = 8 rows of "if tumor looks like X, trial concept Y."
- **Benefit:** Answers the clinician's "so what do I DO for ECO-low patients?" — the paper's translational payoff.
- **Read it:** 8B = drug–gene web; Table 4 = the actionable menu.

### Analysis 11 — Beating TIDE + treatment sketch (Fig. 9, Table 4, Supp. Table S10)
- **Question:** Can ECO beat the reigning predictor, and how would a clinic use all this?
- **Method:** Compare ECO vs **independently precomputed TIDE calls** (we never recompute them — fair fight); cross-tabulate ECO × TMB response rates; draw a **hypothesis-labeled** treatment algorithm.
- **Result:** ECO 0.60 vs TIDE 0.49 (n = 298; combined 0.62). TMB-high respond 36–44% regardless of ECO; TMB-low ~10% regardless → ECO = TMB **substitute**. Fig. 9C: RNA → ECO → high goes to checkpoint blockade, low/stromal goes to stromal-combination trials, re-biopsy at progression.
- **Benefit:** Head-to-head wins get papers cited; the labeled-as-hypothesis sketch shows clinical imagination without overclaiming.
- **Read it:** 9A = ROC face-off; 9B = 2×2 quadrant bars; 9C = flowchart ending in "prospective trial needed."

---

## Chapter 7. Figure & table cheat sheet
| Item | One-line meaning |
|---|---|
| Fig. 1 | The plan: cohorts → score → validation → mechanism → clinic |
| Fig. 2 | Responders score higher in all 6 cohorts |
| Fig. 3 | The balance beats each side alone |
| Fig. 4 | High ECO lives longer, everywhere |
| Fig. 5 | Three tumor states; stromal responds worst |
| Fig. 6 | ECO ≠ TMB; ECO-high looks inflamed |
| Fig. 7 | Two programs, two wirings, two control rooms |
| Fig. 8 | miR-29 brake + drug web |
| Fig. 9 | Beats TIDE; proposed treatment flow |
| Table 1 | Who's in the study |
| Table 2 | All AUCs side by side |
| Table 3 | All survival HRs side by side |
| Table 4 | Drug menu for ECO-low patients |

---

## Chapter 8. Honest limitations + reviewer Q&A (defense practice)
**Limitations (stated in the paper):** all retrospective + computational (no prospective trial, no wet-lab); bulk RNA can't see TLS maturity/location or CAF subtypes (needs spatial/single-cell); small cohorts (16–49) outside IMvigor210; endpoints differ by cohort; only 3 cancers; drug pairings are proposals (remember IDO1's phase III collapse).

**Q: Is AUC 0.62 "good"?** A: Modest alone — but it beats PD-L1 in 6/6, ties IFN-gamma, beats TIDE, and never flips direction in 9 cohorts. Consistency + survival + mechanism > one flashy AUC.
**Q: Why no machine learning?** A: Deliberate. A fixed subtraction can't overfit, any lab can compute it, and clinicians can understand it in one sentence.
**Q: Why median splits, not "optimal" cutoffs?** A: Optimal cutoffs chosen after seeing data inflate results. The median is pre-specified and unbiased; a locked clinical cutoff comes from a future prospective trial.
**Q: Why no multiple-testing penalty across cohorts?** A: We validate ONE pre-specified score (like testing one drug), not hunt among thousands — standard TRIPOD-style validation practice.
**Q: Why does stroma AUC 0.21 count as signal?** A: AUC < 0.5 = inverted prediction: LOW stroma predicts response. That's exactly the biology (walls block drugs).

---

## Chapter 9. Rerun it yourself
```bash
cd /home/user/npj_immunotherapy_paper
python3 code/01_download.py   # fetch data (needs internet)
python3 code/03_inspect.py    # verify files
python3 code/05_parse.py      # standardize cohorts
python3 code/06_imvigor.py    # fetch IMvigor210
python3 code/07_scores.py     # core results
python3 code/07b_refine.py    # CIs, ORs, adjusted models
python3 code/08_figures.py    # Figures 1–6
python3 code/10_enhance_fetch.py  # query Enrichr/STRING/DGIdb
python3 code/11_enhance_analysis.py  # Figures 7–9
python3 code/13_build_docx.py # rebuild manuscript.md + .docx
python3 code/14_audit.py      # style check
```
Each script prints what it's doing; outputs land in `data/`, `results/figures/`, `manuscript/`.

---

## Chapter 10. Glossary (30-second definitions)
**Adjuvant/neoadjuvant:** treatment after/before surgery. **AUC:** 0–1 prediction report card. **Bootstrap CI:** uncertainty estimated by resampling. **CAF:** cancer-associated fibroblast. **CD8:** killer T cell marker. **CI:** confidence interval. **Cox/HR:** survival comparison model/ratio. **CR/PR/SD/PD:** complete/partial response, stable/progressive disease. **DCB:** durable clinical benefit. **EMT:** epithelial→mesenchymal (fibroblast-like) shift. **FPKM/TPM:** RNA abundance units. **Hazard ratio:** relative death rate. **IFN-gamma:** key anti-tumor signal. **Kaplan–Meier:** survival curve. **k-means:** clustering into k groups. **Log-rank:** survival-curve comparison test. **Mann–Whitney:** rank-based group comparison. **MPR:** major pathologic response. **OS/PFS:** overall/progression-free survival. **PD-1/PD-L1:** the off-switch and its blocker target. **RECIST:** tumor-shrinkage rules. **ROC:** tradeoff curve for prediction. **RSEM:** RNA quantification unit. **Spearman:** rank correlation. **TIDE:** rival RNA predictor. **TLS:** tertiary lymphoid structure. **TMB:** tumor mutational burden. **TRRUST/miRTarBase:** regulator/microRNA databases. **z-score:** value in standard-deviation units.

---

## What's next (the road from here)
1. Prospective trial with a locked ECO cutoff (the definitive proof).
2. RNA + multiplex imaging (link ECO to TLS maturity and CAF geography).
3. Single-cell RNA (assign hubs to exact cell types).
4. More cancers/regimens (kidney, head & neck, MSI-high colorectal; CTLA-4 combos).
5. Bench-test miR-29 mimics + CXCR4 blockade in exclusion models.

*End of guide. Re-read Chapter 6 once more — it's the whole paper in 11 steps.*
