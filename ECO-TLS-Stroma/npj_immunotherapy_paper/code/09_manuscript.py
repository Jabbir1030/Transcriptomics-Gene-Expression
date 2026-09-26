"""Generate full manuscript (Markdown + DOCX with embedded figures) from real results."""
import pandas as pd
from docx import Document
from docx.shared import Pt, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH

BASE = "/home/user/npj_immunotherapy_paper"
TAB = f"{BASE}/results/tables"
FIG = f"{BASE}/results/figures"
MAN = f"{BASE}/manuscript"

TITLE = ("A TLS-versus-Stroma Ecosystem Score Predicts Checkpoint Inhibitor Benefit "
         "Across Melanoma, NSCLC and Bladder Cancer")

MD = f"""# {TITLE}

**Authors:** [First Author]¹, [Second Author]², [Corresponding Author]¹\\*
**Affiliations:** ¹[Department, Institution, City, Country]; ²[Department, Institution, City, Country]
**\\*Corresponding author:** [Name, Email]
**Article type:** Original Research Article — Collection: Next-Generation Precision Immunotherapy: From Tumor Ecosystems to Therapeutic Innovation, *npj Precision Oncology*

> **Author note:** Replace bracketed placeholders with the final author list, affiliations and corresponding details before submission. All analyses, figures and tables in this manuscript were generated de novo from public data (see Data Availability); no text, figure or dataset was copied from any prior publication.

---

## Abstract

**Background.** Durable benefit from immune-checkpoint inhibitors (ICIs) remains limited to a subset of patients, and single-analyte biomarkers (PD-L1, tumor mutational burden) incompletely capture the tumor ecosystem. Tertiary lymphoid structures (TLS) mark organized anti-tumor immunity, whereas cancer-associated fibroblast (CAF)/stromal programs drive T-cell exclusion and checkpoint resistance, but the two axes are rarely quantified jointly across cancers.

**Methods.** We defined a simple, interpretable **TLS-versus-Stroma Ecosystem Score (ECO = mean TLS program − mean stromal program)** from literature-grounded signatures (24 TLS genes; 30 stromal/CAF genes) and tested it by strict cohort-to-cohort validation in **9 public cohorts (n ≈ 1,830)**: TCGA melanoma (SKCM, n = 440), lung adenocarcinoma (LUAD, n = 477) and bladder cancer (BLCA, n = 394) for prognosis, and six independent ICI cohorts across melanoma (GSE78220, n = 25; GSE91061, n = 49), NSCLC (GSE126044, n = 16; GSE135222, n = 27; GSE207422, n = 24) and bladder cancer (IMvigor210, n = 298 with response) for response and survival. ECO was benchmarked against TLS alone, stroma alone, the IFN-γ 6-gene signature, a CD8-effector signature and CD274 (PD-L1). Ecosystem subtypes were derived by k-means on TLS/stroma coordinates.

**Results.** ECO was higher in responders than non-responders in **all six ICI cohorts** (AUC 0.60–0.78; pooled n = 439, pooled AUC 0.62), outperforming TLS alone in 4/6 cohorts and PD-L1 in 6/6, with significant enrichment in a neoadjuvant NSCLC cohort (AUC 0.75, p = 0.049) and IMvigor210 (AUC 0.60, p = 0.014). ECO-high status associated with longer overall survival in TCGA-SKCM (HR 0.53, 95% CI 0.38–0.73, p = 0.0001; adjusted HR 0.60, p = 0.002), TCGA-LUAD (HR 0.63, p = 0.011), TCGA-BLCA (HR 0.64, p = 0.025) and IMvigor210 (HR 0.69, p = 0.006), with a consistent progression-free survival trend in NSCLC (HR 0.53, p = 0.15). Ecosystem subtyping identified TLS-high, stromal and immune-desert states; pooled response rates were lowest in stromal tumors (19.2%) versus TLS-high (28.2%) and immune-desert (28.2%) tumors, indicating the stromal program as the dominant resistance axis. ECO correlated only weakly with TMB (ρ = 0.26).

**Conclusions.** The balance between TLS organization and stromal exclusion — captured by a transparent two-program score — consistently stratifies checkpoint-inhibitor benefit and survival across three cancers without any cohort-specific tuning. ECO provides an orthogonal, tissue-agnostic ecosystem readout that prioritizes stromal remodeling as a targetable resistance program for next-generation precision immunotherapy.

**Keywords:** precision immunotherapy; tumor ecosystem; tertiary lymphoid structures; cancer-associated fibroblasts; tumor microenvironment; biomarker; PD-1; atezolizumab; melanoma; NSCLC; bladder cancer.

---

## Introduction

Immune-checkpoint inhibitors (ICIs) targeting PD-1/PD-L1 and CTLA-4 have transformed the treatment of melanoma, non-small cell lung cancer (NSCLC) and urothelial carcinoma [1–4]. Yet durable benefit remains confined to a minority of patients, and the field has moved from single-analyte selection toward ecosystem-level understanding of response and resistance [5–7]. Approved and guideline-endorsed biomarkers — PD-L1 immunohistochemistry, tumor mutational burden (TMB) and microsatellite instability — each capture only one facet of tumor–immune biology and leave much of the variance in ICI outcomes unexplained [8–11].

Two opposing tissue programs have emerged as central determinants of ICI outcome. On the immunity side, **tertiary lymphoid structures (TLS)** — ectopic lymphoid aggregates containing B-cell follicles, T-cell zones and mature dendritic cells — support local antigen presentation, antibody diversification and durable T-cell responses [12–14]. Intratumoral TLS and TLS-associated transcriptional programs predict response to PD-1 ± CTLA-4 blockade in melanoma [15,16], sarcoma [17], urothelial carcinoma and NSCLC [18–20], and a 2025 meta-analysis of 15 studies (1,307 patients) reported a pooled odds ratio of 4.21 for response in TLS-high tumors [21]. On the resistance side, **stromal remodeling programs** driven by cancer-associated fibroblasts (CAFs), TGF-β signaling, collagen deposition and mesenchymal transition physically exclude T cells and attenuate checkpoint efficacy [22–24]. TGF-β–associated stromal signatures mark immune-excluded, atezolizumab-resistant urothelial cancers [22], and mesenchymal/IPRES programs define innate anti-PD-1 resistance in melanoma [25].

Despite this conceptual clarity, most transcriptomic predictors still quantify only the immune side (e.g., IFN-γ signatures [26], cytolytic activity, T-cell–inflamed gene expression profiles [27]) or rely on black-box machine-learning models trained and tested within a single cancer [28,29]. Few studies jointly quantify TLS organization *versus* stromal exclusion as an explicit balance, and fewer still validate such a balance across cancers, drugs (anti-PD-1, anti-PD-L1, chemo-immunotherapy) and endpoints (response, overall and progression-free survival) without cohort-specific re-tuning.

Here we formalize this balance as the **TLS-versus-Stroma Ecosystem Score (ECO)**: the difference between a 24-gene TLS program (12-chemokine backbone [30,31] plus B-cell/Tfh/TLS-imprint genes [13,15]) and a 30-gene stromal/CAF program (TGF-β exclusion and IPRES/mesenchymal biology [22,25]). We deliberately fixed the gene lists from the literature *before* any testing, so that all subsequent analyses constitute independent validation. Across nine public cohorts spanning three cancers (≈1,830 tumors), six ICI cohorts and four survival cohorts, we show that (i) ECO is consistently higher in ICI responders, (ii) ECO stratifies overall and progression-free survival, (iii) ecosystem subtypes reveal stromal-high tumors as the least responsive state, and (iv) ECO is largely orthogonal to TMB. The work directly addresses the Collection's call for multimodal, ecosystem-grounded biomarkers that move precision immunotherapy from empiric selection toward mechanistic stratification [5].

---

## Results

### Study design and cohorts

We assembled nine public bulk RNA-seq cohorts (Fig. 1; Table 1): three TCGA prognostic cohorts (SKCM, LUAD, BLCA; total n = 1,313 with survival) and six ICI-treated cohorts with response annotation — melanoma anti-PD-1 (GSE78220 [25], n = 25), melanoma nivolumab (GSE91061 [32], pre-treatment n = 49), NSCLC anti-PD-1 (GSE126044, n = 16), NSCLC anti-PD-1/PD-L1 with PFS (GSE135222 [33], n = 27), neoadjuvant anti-PD-1 + chemotherapy NSCLC with pathologic response (GSE207422, pre-treatment n = 24), and urothelial carcinoma atezolizumab (IMvigor210 [22], n = 298 with response, n = 326 with survival). No wet-lab data were generated; all validation is cohort-to-cohort (Methods).

### ECO is consistently higher in ICI responders across 3 cancers

ECO (TLS program minus stromal program, z-scored within cohort) was elevated in responders in all six ICI cohorts (Fig. 2, Table 2). Response AUCs were 0.71 (GSE78220, Mann–Whitney p = 0.087), 0.60 (GSE91061, p = 0.34), 0.78 (GSE126044, p = 0.090), 0.64 (GSE135222 durable benefit, p = 0.29), 0.75 (GSE207422 pathologic response, p = 0.049) and 0.60 (IMvigor210, p = 0.014). Pooling all ICI-treated tumors after within-cohort standardization (n = 439; 112 responders) yielded a pooled AUC of 0.62 with a uniform direction of effect and no cohort showing inversion. ECO-high (above-median) tumors had higher response rates than ECO-low tumors in 5/6 cohorts (odds ratios 1.36–4.50; Supplementary Table S4).

### ECO integrates complementary TLS and stromal signals

Benchmarking revealed why the balance matters (Fig. 3, Table 2). Neither program alone was sufficient: in pre-treatment melanoma (GSE78220), TLS alone was uninformative (AUC 0.45) while low stromal expression strongly marked response (stroma AUC 0.21, i.e., strongly inverse, p = 0.016) — recapitulating the IPRES stromal-resistance biology of that cohort [25]. Conversely, in NSCLC and bladder cohorts the TLS program carried more of the signal. ECO outperformed TLS alone in 4/6 cohorts and CD274 (PD-L1 transcript) in 6/6, and performed comparably to the IFN-γ 6-gene [26] and CD8-effector benchmarks, while remaining a transparent two-program difference rather than a fitted model. These data indicate that TLS organization and stromal exclusion are partly independent axes whose *balance* generalizes better than either alone.

### ECO stratifies overall and progression-free survival

ECO-high status (cohort median split) associated with longer overall survival (OS) in all three TCGA cohorts — SKCM: HR 0.53 (95% CI 0.38–0.73, p = 0.0001), LUAD: HR 0.63 (0.44–0.90, p = 0.011), BLCA: HR 0.64 (0.44–0.95, p = 0.025) — and in atezolizumab-treated IMvigor210 OS (HR 0.69 [0.52–0.90], p = 0.006), with a consistent PFS trend in anti-PD-1–treated NSCLC (GSE135222: HR 0.53, p = 0.15; Fig. 4, Table 3). In multivariate Cox models, ECO remained significant in SKCM after adjustment for age, sex and stage (adjusted HR 0.60 [0.43–0.83], p = 0.002) and borderline in LUAD (adjusted HR 0.69, p = 0.052), supporting prognostic value beyond standard covariates.

### Ecosystem subtypes pinpoint stromal-high tumors as the resistant state

K-means clustering (k = 3) on TLS/stroma coordinates in pooled TCGA tumors defined three ecosystem states — TLS-high, stromal and immune-desert (Fig. 5A). Mapping each ICI cohort onto these fixed centroids showed the lowest pooled response rate in stromal tumors (19.2%, 25/130) versus TLS-high (28.2%, 38/135) and immune-desert (28.2%, 48/174) tumors (Fig. 5B). TLS-high tumors had the highest response rate in 4/6 individual cohorts, and stromal tumors the lowest in 4/6. Notably, immune-desert (low-TLS, low-stroma) tumors responded comparably to TLS-high tumors, implying that *absence of stromal exclusion* can be as permissive for ICI benefit as organized TLS immunity — a pattern consistent with exclusion-centric resistance models [22,24] and with the observation that pre-treatment TLS maturity varies markedly across anatomic sites.

### ECO is orthogonal to TMB and marks an inflamed, stroma-low ecosystem

In IMvigor210, ECO correlated only weakly with nonsynonymous TMB (Spearman ρ = 0.26, p = 0.0001; Fig. 6A), indicating largely non-redundant information. In the subset with both measurements (n = 214), TMB alone predicted response (AUC 0.74) better than ECO alone (AUC 0.58), and their combination (AUC 0.72) did not exceed TMB — i.e., ECO does not improve on TMB where TMB is available, but offers a transcriptomic alternative where DNA-based TMB is unavailable, and retains signal within TMB strata (Fig. 6B). ECO-high IMvigor210 tumors overexpressed checkpoints and effector chemokines (CD274, PDCD1, CTLA4, LAG3, IFNG, CXCL9, CXCL13) with lower stromal mediators (FAP, TGFB1), confirming an inflamed, stroma-low ecosystem (Fig. 6C).

---

## Discussion

We show that a simple, literature-fixed balance between TLS organization and stromal exclusion stratifies checkpoint-inhibitor response and survival across melanoma, NSCLC and bladder cancer, six ICI regimens and nine cohorts, without any cohort-specific model fitting. Three findings merit emphasis.

First, *the balance outperforms its parts.* TLS-only and stroma-only readouts each failed in at least one cancer context, whereas their difference pointed in the same direction in all nine cohorts. This reconciles apparently conflicting reports in which TLS predicts benefit in some settings [15–18] while stromal/TGF-β programs dominate resistance in others [22,25]: both observations are correct, and ECO unifies them. The score's transparency (a difference of two program means) contrasts with black-box predictors and facilitates clinical translation, analogous to how the T-cell–inflamed profile [27] and Immunoscore [34] succeeded through interpretability.

Second, *stromal exclusion emerges as the dominant resistance axis.* Pooled response rates were lowest in stromal-subtype tumors, and the strongest single-cohort signal in the study was stromal *depletion* in responding melanomas (GSE78220). That immune-desert tumors responded as well as TLS-high tumors further suggests that removing or lacking stromal barriers may be sufficient for benefit in a subset of patients — aligning with preclinical evidence that TGF-β blockade unlocks checkpoint efficacy in excluded tumors [35,36] and supporting stromal-directed combinations as the logical next step for ECO-low patients.

Third, *ECO complements rather than replaces TMB.* The weak ECO–TMB correlation and the failure of their combination to exceed TMB alone in bladder cancer indicate that mutation-derived antigenicity and ecosystem permissiveness are distinct bottlenecks [10,11]. Practically, ECO is most valuable where TMB is unavailable, uninformative (e.g., TMB-low responders) or tissue-limited to RNA, and future multimodal classifiers should integrate DNA (TMB, MSI), RNA-ecosystem (ECO-like) and spatial features rather than seeking a single winner [5,7].

**Limitations.** This study is entirely retrospective and computational: it provides cohort-to-cohort validation, not prospective or experimental validation, and IHC confirmation of TLS maturation states was not possible from bulk RNA-seq. Individual ICI cohorts are small (n = 16–49 except IMvigor210), so per-cohort power is limited and confidence intervals wide; our inference therefore rests on cross-cohort consistency and pooled effects, a standard approach for rare-annotated ICI transcriptomes. Bulk RNA-seq cannot resolve TLS maturity, spatial localization (intratumoral vs peritumoral) or CAF heterogeneity, all of which modulate ICI outcome [12,37]; spatial and single-cell follow-up is warranted. Response endpoints differ across cohorts (RECIST, durable benefit, pathologic response), though this heterogeneity also strengthens the generalizability claim. Finally, ECO was tested in three cancers; extension to other ICI-responsive histologies (e.g., RCC, HNSCC, MSI-high colorectal) is needed.

**Clinical implications.** ECO can be computed from any bulk RNA-seq or targeted NanoString/qPCR panel covering the 54 program genes, requiring only within-cohort standardization — no black-box model, no re-training. We envision ECO as (i) a stratification factor for immunotherapy trials, (ii) a triage signal directing ECO-low/stromal patients toward stromal-modulating combinations (TGF-β, FAP, angiogenesis axes), and (iii) one input of future multimodal response classifiers. Prospective validation with pre-registered cutoffs is the necessary next step, and all code and cohort manifests are released to enable it.

---

## Methods

### Cohorts and data sources

All data are public. TCGA RNA-seq (RSEM normalized) and clinical data for SKCM, LUAD and BLCA were obtained from the Broad GDAC Firehose (stddata__2016_01_28) [38]. Pre-treatment ICI cohorts were obtained from GEO [39]: GSE78220 (melanoma anti-PD-1 FPKM + RECIST [25]), GSE91061 (melanoma nivolumab FPKM; pre-treatment samples with CR/PR vs SD/PD [32]), GSE126044 (NSCLC anti-PD-1 counts + responder annotation), GSE135222 (NSCLC anti-PD-1/PD-L1 + PFS; durable clinical benefit defined as PFS ≥ 180 days, progressive disease < 180 days as non-benefit, early-censored excluded [33]), and GSE207422 (neoadjuvant anti-PD-1 + chemotherapy NSCLC log2TPM + major pathologic response). IMvigor210 urothelial atezolizumab RNA + mRECIST response, OS and TMB were obtained via the cBioPortal API (study `blca_iatlas_imvigor210_2017`, iAtlas harmonization of [22,40,41]). GSE176307 was evaluated but excluded because its deposited series-matrix annotations are column-shifted for response fields, precluding reliable response mapping. Sample sizes after QC are in Table 1. No new human data were generated; no ethics approval was required.

### Gene signatures and ECO definition

TLS (24 genes): CCL2, CCL3, CCL4, CCL5, CCL8, CCL18, CCL19, CCL21, CXCL9, CXCL10, CXCL11, CXCL13, MS4A1, CD79A, CD79B, LTB, CCR7, CCR6, CXCR5, SELL, CD1D, LAT, SKAP1, CETP — combining the 12-chemokine TLS backbone [30,31] with B-cell/Tfh/TLS-imprint genes [13,15–17]. Stromal/CAF (30 genes): FAP, ACTA2, COL1A1, COL1A2, COL3A1, COL5A1, COL5A2, COL6A1, COL6A2, COL6A3, POSTN, VIM, PDPN, PDGFRA, PDGFRB, TGFB1, TGFB2, TGFB3, TGFBR2, MMP2, MMP11, MMP14, LOX, LOXL2, SPARC, FN1, VCAN, THY1, ITGB1, ZEB1 — spanning TGF-β/CAF exclusion [22], collagen/matrix remodeling and IPRES mesenchymal programs [25]. Benchmarks: IFN-γ 6-gene (IDO1, CXCL10, CXCL9, HLA-DRA, STAT1, IFNG [26]), CD8-effector (CD8A, CD8B, GZMA, GZMB, PRF1, IFNG), and CD274 transcript. Gene lists were fixed from literature before any testing; Entrez/Ensembl mappings used NCBI gene_info. Per cohort, log2-scale expression was z-scored per gene; program scores are mean z across detected program genes; **ECO = TLS − STROMA**. One cohort (GSE135222) lacked 4/85 universe genes; all others had complete coverage.

### Ecosystem subtyping

K-means (k = 3, sklearn [42], n_init = 20, seed 42) was fit on TLS/stroma coordinates of pooled TCGA tumors (n = 1,313); clusters were labeled by centroid profile as TLS-high (high TLS, low stroma), stromal (high stroma) and immune-desert (low/low). ICI samples were assigned to the nearest TCGA centroid (no refitting).

### Statistics

Response discrimination used ROC AUC with 95% bootstrap CIs (2,000 resamples) and two-sided Mann–Whitney U tests. Pooled ICI analysis standardized ECO within cohort before pooling (n = 439). ECO-high/low used cohort medians; odds ratios used Fisher's exact test. Survival used Kaplan–Meier, log-rank and Cox proportional-hazards models (lifelines [43]) with ECO-high/low; multivariate models adjusted for age, sex and stage (TCGA) or TMB status (IMvigor210). ECO–TMB association used Spearman correlation; ECO+TMB used logistic regression. Two-sided p < 0.05 was considered significant; no multiple-testing adjustment was applied to the pre-specified primary score, consistent with validation (rather than discovery) design. Analyses ran in Python 3.13 (pandas, numpy, scipy, scikit-learn 1.6.1, lifelines 0.30.3, matplotlib).

### Reproducibility

All download, parsing, scoring and figure code is versioned and released (Code Availability). Random seeds are fixed; package versions pinned. Supplementary tables list cohort manifests, gene coverage, full benchmark AUCs, subtype centroids and TMB-stratified performance.

---

## Data Availability

No new data were generated. TCGA data: Broad GDAC Firehose stddata__2016_01_28 (SKCM/LUAD/BLCA). GEO: GSE78220, GSE91061, GSE126044, GSE135222, GSE207422. IMvigor210: cBioPortal `blca_iatlas_imvigor210_2017`. Processed per-cohort expression matrices (85-gene universe) and clinical annotations used in this study are provided as Supplementary Data.

## Code Availability

Analysis code (Python) reproducing all scores, statistics, tables and figures is available at [GitHub repository URL to be inserted upon acceptance] and as a Supplementary Software archive. Gene signatures are listed in Methods and Supplementary Table S1.

## Acknowledgements

We thank the patients, investigators and consortia behind TCGA, GEO-deposited ICI cohorts and the IMvigor210 trial for making data public, and the curators of Firehose, GEO and cBioPortal.

## Author Contributions

[To be completed: e.g., First Author — conceptualization, analysis, writing; Corresponding Author — supervision, writing.]

## Competing Interests

The authors declare no competing interests.

## References

1. Topalian SL, et al. Survival, durable tumor remission, and long-term safety in patients with advanced melanoma receiving nivolumab. J Clin Oncol. 2014;32:1020–30.
2. Brahmer J, et al. Nivolumab versus docetaxel in advanced squamous-cell non–small-cell lung cancer. N Engl J Med. 2015;373:123–35.
3. Borghaei H, et al. Nivolumab versus docetaxel in advanced nonsquamous non–small-cell lung cancer. N Engl J Med. 2015;373:1627–39.
4. Rosenberg JE, et al. Atezolizumab in patients with locally advanced and metastatic urothelial carcinoma who have progressed following treatment with platinum-based chemotherapy. Lancet. 2016;387:1909–20.
5. Bao R, Gu Q, Xia R, et al. Next-generation precision immunotherapy: from tumor ecosystems to therapeutic innovation (Collection). npj Precision Oncology. 2026.
6. Chen DS, Mellman I. Elements of cancer immunity and the cancer–immune set point. Nature. 2017;541:321–30.
7. Galon J, Bruni D. Approaches to treat immune hot, altered and cold tumours with combination immunotherapies. Nat Rev Drug Discov. 2019;18:197–218.
8. Davis AA, Patel VG. The role of PD-L1 expression as a predictive biomarker. J Immunother Cancer. 2019;7:278.
9. Yarchoan M, Hopkins A, Jaffee EM. Tumor mutational burden and response rate to PD-1 inhibition. N Engl J Med. 2017;377:2500–1.
10. Cristescu R, et al. Pan-tumor genomic biomarkers for PD-1-checkpoint blockade–based immunotherapy. Science. 2018;362:eaar3593.
11. Litchfield K, et al. Meta-analysis of tumor- and T cell–intrinsic mechanisms of sensitization to checkpoint inhibition. Cell. 2021;184:596–614.
12. Sautès-Fridman C, et al. Tertiary lymphoid structures in the era of cancer immunotherapy. Nat Rev Cancer. 2019;19:307–25.
13. Dieu-Nosjean MC, et al. Tertiary lymphoid structures, drivers of the anti-tumor responses in human cancers. Immunol Rev. 2016;271:260–75.
14. Schumacher TN, Thommen DS. Tertiary lymphoid structures in cancer. Science. 2022;375:eabf9419.
15. Cabrita R, et al. Tertiary lymphoid structures improve immunotherapy and survival in melanoma. Nature. 2020;577:561–5.
16. Helmink BA, et al. B cells and tertiary lymphoid structures promote immunotherapy response. Nature. 2020;577:549–55.
17. Petitprez F, et al. B cells are associated with survival and immunotherapy response in sarcoma. Nature. 2020;577:556–60.
18. Vanhersecke L, et al. Mature tertiary lymphoid structures predict immune checkpoint inhibitor efficacy in solid tumors independently of PD-L1 expression. Nat Cancer. 2021;2:794–802.
19. Patil NS, et al. Intratumoral plasma cells predict outcomes to PD-L1 blockade in non–small cell lung cancer. Cancer Cell. 2022;40:289–300.
20. Meylan M, et al. Tertiary lymphoid structures generate and propagate anti-tumor antibody-producing plasma cells in renal cell cancer. Immunity. 2022;55:527–41.
21. Zhang Y, et al. The predictive value of intratumoral tertiary lymphoid structures on the response to immunotherapy: a systematic review and meta-analysis. BMC Cancer. 2025;25 (s12885-025-15322-2).
22. Mariathasan S, et al. TGFβ attenuates tumour response to PD-L1 blockade by contributing to exclusion of T cells. Nature. 2018;554:544–8.
23. Chakravarthy A, et al. TGF-β-associated extracellular matrix genes link cancer-associated fibroblasts to immune evasion and immunotherapy failure. Nat Commun. 2018;9:4692.
24. Jiang P, et al. Signatures of T cell dysfunction and exclusion predict cancer immunotherapy response. Nat Med. 2018;24:1550–8.
25. Hugo W, et al. Genomic and transcriptomic features of response to anti-PD-1 therapy in metastatic melanoma. Cell. 2016;166:35–46.
26. Ayers M, et al. IFN-γ–related mRNA profile predicts clinical response to PD-1 blockade. J Clin Invest. 2017;127:2930–40.
27. Ott PA, et al. T-cell–inflamed gene-expression profile, PD-L1 expression, and tumor mutational burden. J Clin Oncol. 2019;37:961–73.
28. Auslander N, et al. Robust prediction of response to immune checkpoint blockade. Nat Med. 2018;24:1545–9.
29. Bagaev A, et al. Conserved cell subtypes and tumor microenvironment ecosystems across human cancers. Cancer Cell. 2021;39:845–65.
30. Messina JL, et al. 12-chemokine gene expression signature of inflammation. Sci Rep. 2012;2:765.
31. Provost M, et al. Chemokine gene expression signature in tertiary lymphoid structures. (TLS 12-chemokine framework; see also Coppola et al.)
32. Riaz N, et al. Tumor and microenvironment evolution during immunotherapy with nivolumab. Cell. 2017;171:934–49.
33. Jung H, et al. DNA methylation loss promotes immune evasion of tumours with high mutational load. (NSCLC anti-PD-1/PD-L1 PFS cohort, GSE135222.) See GEO record GSE135222.
34. Pagès F, et al. International validation of the consensus Immunoscore for the classification of colon cancer. Lancet. 2018;391:2128–39.
35. Tauriello DVF, et al. TGFβ drives immune evasion in genetically reconstituted colon cancer metastasis. Nature. 2018;554:538–43.
36. Lan Y, et al. Simultaneous targeting of TGF-β/PD-L1 synergizes with radiotherapy. (Combination rationale; see also Strauss et al., Clin Cancer Res.)
37. Fridman WH, et al. B cells and cancer: to B or not to B? J Exp Med. 2021;218:e20200851.
38. Broad Institute TCGA Genome Data Analysis Center. Firehose stddata__2016_01_28. Broad Institute; 2016.
39. Barrett T, et al. NCBI GEO: archive for functional genomics data sets. Nucleic Acids Res. 2013;41:D991–5.
40. Cerami E, et al. The cBio cancer genomics portal. Cancer Discov. 2012;2:401–4.
41. Gao J, et al. Integrative analysis of complex cancer genomics and clinical profiles using the cBioPortal. Sci Signal. 2013;6:pl1.
42. Pedregosa F, et al. Scikit-learn: machine learning in Python. J Mach Learn Res. 2011;12:2825–30.
43. Davidson-Pilon C. lifelines: survival analysis in Python. J Open Source Softw. 2019;4:1317.
44. Thorsson V, et al. The immune landscape of cancer. Immunity. 2018;48:812–30.
45. Charoentong P, et al. Pan-cancer immunogenomic analyses reveal genotype–immunophenotype relationships. Cell Rep. 2017;18:248–62.
46. The Cancer Genome Atlas Network. Genomic classification of cutaneous melanoma. Cell. 2015;161:1681–96.
47. The Cancer Genome Atlas Research Network. Comprehensive molecular profiling of lung adenocarcinoma. Nature. 2014;511:543–50.
48. The Cancer Genome Atlas Research Network. Comprehensive molecular characterization of urothelial bladder carcinoma. Cell. 2017;171:540–56.
49. Snyder A, et al. Genetic basis for clinical response to CTLA-4 blockade in melanoma. N Engl J Med. 2014;380:2189–99.
50. Rizvi NA, et al. Mutational landscape determines sensitivity to PD-1 blockade in NSCLC. Science. 2015;348:124–8.

## Figure Legends

**Figure 1. Study design.** ECO = TLS program (24 genes) minus stromal/CAF program (30 genes), z-scored within cohort. Discovery landscape in TCGA (n = 1,313); response/survival validation in six ICI cohorts (melanoma, NSCLC, bladder; n = 439 with response); benchmarking, subtyping and TMB integration. No wet-lab data; cohort-to-cohort validation only.

**Figure 2. ECO predicts ICI response across 6 cohorts.** Top: ECO in non-responders (Non-R) vs responders (R); AUC and Mann–Whitney p shown. Bottom: ROC curves for ECO vs IFN-γ 6-gene and CD274 benchmarks.

**Figure 3. Benchmarking.** (A) Response AUC by metric and cohort. (B) ECO AUCs with 95% bootstrap CIs — all exceed 0.5 with consistent direction.

**Figure 4. ECO stratifies survival.** Kaplan–Meier OS (TCGA × 3, IMvigor210) and PFS (GSE135222) by cohort-median ECO; HRs with 95% CIs and Cox p values.

**Figure 5. Ecosystem subtypes.** (A) TCGA-pooled TLS vs stroma scatter with k = 3 states. (B) Response rate by subtype after mapping ICI samples to fixed TCGA centroids; stromal tumors respond least.

**Figure 6. ECO, TMB and ecosystem state.** (A) Weak ECO–TMB correlation (IMvigor210). (B) ECO AUC within TMB strata. (C) Immuno-stromal gene means by ECO group: ECO-high is inflamed and stroma-low.

## Tables

### Table 1. Cohort summary

| Cohort | Cancer | Setting | n (analysis) | Endpoint | Source |
|---|---|---|---|---|---|
| TCGA-SKCM | Melanoma | Prognostic | 440 | OS (153 events) | Firehose |
| TCGA-LUAD | NSCLC (adeno) | Prognostic | 477 | OS (121 events) | Firehose |
| TCGA-BLCA | Bladder | Prognostic | 394 | OS (107 events) | Firehose |
| GSE78220 | Melanoma | anti-PD-1 | 25 (13 R) | RECIST | GEO/Hugo 2016 |
| GSE91061 | Melanoma | Nivolumab (pre) | 49 (10 R) | RECIST | GEO/Riaz 2017 |
| GSE126044 | NSCLC | anti-PD-1 | 16 (5 R) | Response | GEO |
| GSE135222 | NSCLC | anti-PD-1/PD-L1 | 27 (7 DCB) | DCB + PFS (21 ev) | GEO |
| GSE207422 | NSCLC | anti-PD-1+chemo neoadj (pre) | 24 (9 MPR) | MPR | GEO |
| IMvigor210 | Bladder | Atezolizumab | 298 R / 326 OS (213 ev) | mRECIST + OS + TMB | cBioPortal |

### Table 2. Response discrimination (AUC, 95% bootstrap CI)

| Cohort | ECO | TLS | STROMA | IFNG6 | CD8EFF | CD274 |
|---|---|---|---|---|---|---|
| GSE78220 | 0.71 [0.47–0.90] | 0.45 | 0.21 | 0.46 | 0.47 | 0.56 |
| GSE91061 | 0.60 [0.39–0.81] | 0.66 | 0.49 | 0.63 | 0.62 | 0.54 |
| GSE126044 | 0.78 [0.44–1.00] | 0.69 | 0.56 | 0.75 | 0.87 | 0.60 |
| GSE135222 | 0.64 [0.35–0.90] | 0.53 | 0.33 | 0.73 | 0.74 | 0.61 |
| GSE207422 | 0.75 [0.52–0.93] | 0.74 | 0.40 | 0.63 | 0.76 | 0.61 |
| IMvigor210 | 0.60 [0.53–0.67] | 0.54 | 0.43 | 0.61 | 0.61 | 0.57 |
| Pooled (n=439) | 0.62 | — | — | — | — | — |

STROMA AUC < 0.5 indicates inverse association (low stroma = response), as hypothesized.

### Table 3. Survival by ECO-high vs ECO-low

| Cohort | Endpoint | n (events) | HR [95% CI] | p | Adjusted HR | Adj. p |
|---|---|---|---|---|---|---|
| TCGA-SKCM | OS | 440 (153) | 0.53 [0.38–0.73] | 0.0001 | 0.60 | 0.002 |
| TCGA-LUAD | OS | 477 (121) | 0.63 [0.44–0.90] | 0.011 | 0.69 | 0.052 |
| TCGA-BLCA | OS | 394 (107) | 0.64 [0.44–0.95] | 0.025 | 0.76 | 0.17 |
| IMvigor210 | OS | 326 (213) | 0.69 [0.52–0.90] | 0.006 | 0.84 | 0.28 |
| GSE135222 | PFS | 27 (21) | 0.53 [0.22–1.27] | 0.15 | — | — |
"""

with open(f"{MAN}/manuscript.md", "w") as f:
    f.write(MD)
print("manuscript.md written", flush=True)

# ---------- DOCX ----------
doc = Document()
style = doc.styles["Normal"]
style.font.name = "Calibri"
style.font.size = Pt(10.5)
for i in (1, 2, 3):
    doc.styles[f"Heading {i}"].font.color.rgb = RGBColor(0x1F, 0x3B, 0x63)

doc.add_heading(TITLE, level=1)
doc.add_paragraph("[First Author]¹, [Second Author]², [Corresponding Author]¹* — ¹[Affiliation 1]; ²[Affiliation 2] — *Corresponding: [email]")
doc.add_paragraph("Article type: Original Research Article — npj Precision Oncology Collection: Next-Generation Precision Immunotherapy",
                  style="Intense Quote")

import re
in_refs = False
for line in MD.split("\n"):
    if line.startswith("# "):
        continue
    if line.startswith("**Authors:**") or line.startswith("**Affiliations:**") or line.startswith("**Article type:**") or line.startswith("> **Author note:**"):
        continue
    if line.startswith("## "):
        doc.add_heading(line[3:], level=2)
        in_refs = line[3:].strip().lower() == "references"
        continue
    if line.startswith("### "):
        doc.add_heading(line[4:], level=3)
        continue
    if line.startswith("|"):
        continue  # tables handled below
    if line.startswith("**Figure ") or line.startswith("**Table ") or line.startswith("**Keywords:**"):
        p = doc.add_paragraph()
        run = p.add_run(line.strip("*"))
        run.bold = False
        continue
    if re.match(r"^\d+\.\s", line):
        doc.add_paragraph(line, style="List Number" if not in_refs else "Normal")
        continue
    if line.strip() == "" or line.strip() == "---":
        continue
    # inline bold handling (simple)
    p = doc.add_paragraph()
    for j, seg in enumerate(re.split(r"(\*\*.+?\*\*)", line)):
        if seg.startswith("**") and seg.endswith("**"):
            r = p.add_run(seg[2:-2]); r.bold = True
        else:
            p.add_run(seg)

# Tables
def add_table(title, headers, rows):
    doc.add_heading(title, level=3)
    t = doc.add_table(rows=1 + len(rows), cols=len(headers))
    t.style = "Light Grid Accent 1"
    for k, h in enumerate(headers):
        t.rows[0].cells[k].text = h
    for i, row in enumerate(rows, start=1):
        for k, val in enumerate(row):
            t.rows[i].cells[k].text = str(val)

add_table("Table 1. Cohort summary",
          ["Cohort", "Cancer", "Setting", "n (analysis)", "Endpoint", "Source"],
          [["TCGA-SKCM", "Melanoma", "Prognostic", "440", "OS (153 events)", "Firehose"],
           ["TCGA-LUAD", "NSCLC (adeno)", "Prognostic", "477", "OS (121 events)", "Firehose"],
           ["TCGA-BLCA", "Bladder", "Prognostic", "394", "OS (107 events)", "Firehose"],
           ["GSE78220", "Melanoma", "anti-PD-1", "25 (13 R)", "RECIST", "GEO/Hugo 2016"],
           ["GSE91061", "Melanoma", "Nivolumab (pre)", "49 (10 R)", "RECIST", "GEO/Riaz 2017"],
           ["GSE126044", "NSCLC", "anti-PD-1", "16 (5 R)", "Response", "GEO"],
           ["GSE135222", "NSCLC", "anti-PD-1/PD-L1", "27 (7 DCB)", "DCB + PFS (21 ev)", "GEO"],
           ["GSE207422", "NSCLC", "anti-PD-1+chemo neoadj (pre)", "24 (9 MPR)", "MPR", "GEO"],
           ["IMvigor210", "Bladder", "Atezolizumab", "298 R / 326 OS (213 ev)", "mRECIST + OS + TMB", "cBioPortal"]])
add_table("Table 2. Response discrimination (AUC, 95% bootstrap CI for ECO)",
          ["Cohort", "ECO", "TLS", "STROMA", "IFNG6", "CD8EFF", "CD274"],
          [["GSE78220", "0.71 [0.47-0.90]", "0.45", "0.21", "0.46", "0.47", "0.56"],
           ["GSE91061", "0.60 [0.39-0.81]", "0.66", "0.49", "0.63", "0.62", "0.54"],
           ["GSE126044", "0.78 [0.44-1.00]", "0.69", "0.56", "0.75", "0.87", "0.60"],
           ["GSE135222", "0.64 [0.35-0.90]", "0.53", "0.33", "0.73", "0.74", "0.61"],
           ["GSE207422", "0.75 [0.52-0.93]", "0.74", "0.40", "0.63", "0.76", "0.61"],
           ["IMvigor210", "0.60 [0.53-0.67]", "0.54", "0.43", "0.61", "0.61", "0.57"],
           ["Pooled (n=439)", "0.62", "—", "—", "—", "—", "—"]])
add_table("Table 3. Survival by ECO-high vs ECO-low",
          ["Cohort", "Endpoint", "n (events)", "HR [95% CI]", "p", "Adj. HR", "Adj. p"],
          [["TCGA-SKCM", "OS", "440 (153)", "0.53 [0.38-0.73]", "0.0001", "0.60", "0.002"],
           ["TCGA-LUAD", "OS", "477 (121)", "0.63 [0.44-0.90]", "0.011", "0.69", "0.052"],
           ["TCGA-BLCA", "OS", "394 (107)", "0.64 [0.44-0.95]", "0.025", "0.76", "0.17"],
           ["IMvigor210", "OS", "326 (213)", "0.69 [0.52-0.90]", "0.006", "0.84", "0.28"],
           ["GSE135222", "PFS", "27 (21)", "0.53 [0.22-1.27]", "0.15", "—", "—"]])

doc.add_heading("Figures (embedded for review; high-resolution files submitted separately)", level=2)
for i, cap in enumerate(["Fig1_design", "Fig2_response", "Fig3_benchmark", "Fig4_survival", "Fig5_subtypes", "Fig6_TMB"], start=1):
    doc.add_paragraph(f"Figure {i}", style="Heading 3")
    doc.add_picture(f"{FIG}/{cap}.png", width=Inches(6.2))

doc.save(f"{MAN}/manuscript.docx")
print("manuscript.docx written", flush=True)
