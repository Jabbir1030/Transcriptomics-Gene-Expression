# npj Precision Oncology — Next-Generation Precision Immunotherapy
## Full Paper Protocol: 100% Computational, Submission-Ready, Original Analysis

**Collection:** Next-Generation Precision Immunotherapy: From Tumor Ecosystems to Therapeutic Innovation
**Deadline:** 16 June 2027 (~9 months — comfortable timeline)
**Article type:** Original Research Article
**Status:** PROTOCOL v1.0 — awaiting user confirmation of cancer focus

---

## 1. Why this Collection is 100% feasible computationally

Unlike CAR-T (needs cell engineering + trials) or ctDNA MRD (needs prospective plasma + ultra-deep sequencing),
this Collection explicitly invites:

- "AI-based biomarker discovery"
- "Computational pathology, radiology"
- "Integrated multi-omic patient stratification"
- "Mechanisms of response and resistance / tumor-immune-stromal crosstalk"

All of these can be done rigorously with **public RNA-seq + clinical + mutation data + single-cell validation**.
No wet lab, no patient recruitment, no ethics approval required.

---

## 2. Proposed paper (recommended option)

### Working title
**An Explainable Multimodal Tumor Ecosystem Classifier Integrating Stromal Remodeling and Tertiary Lymphoid Structure Programs Predicts Immunotherapy Response and Survival in Melanoma with Cross-Cancer Validation in NSCLC and Urothelial Cancer**

### Core research question
Can tumor-immune-stromal ecosystem states — quantified from bulk RNA-seq and decoded with explainable machine learning — predict anti-PD-1/CTLA-4 response better than single analytes (PD-L1, TMB alone), reveal a targetable stromal-exclusion resistance axis, and generalize across cancers?

### Why it is novel (gap analysis, Sep 2026)
- Most published signatures are single-cohort, single-gene-list, black-box ML, no ecosystem logic.
- This paper: (a) ecosystem-subtype framework (TLS-high vs Stromal-Excluded vs Desert), (b) tumor-stromal ligand–receptor crosstalk, (c) explainable AI (SHAP) + nomogram, (d) discovery in melanoma + independent validation in 5 immunotherapy cohorts across 3 cancers, (e) scRNA-seq orthogonal validation, (f) in silico drug repurposing for the resistant subtype → matches "From Tumor Ecosystems to Therapeutic Innovation" exactly.
- Recent npj Precision Oncology precedent (Yu et al. 2026, ovarian ML + scRNA + drug design) shows this exact computational style is accepted.

### Alternative titles (pick one)
1. Tumor Ecosystem Subtypes Define Immunotherapy Response and Stromal Resistance Programs in Melanoma
2. A TLS-versus-Stroma Ecosystem Score Predicts Checkpoint Inhibitor Benefit Across Melanoma, NSCLC and Bladder Cancer

---

## 3. Data sources — all public, free, citable (no plagiarism, no fabrication)

| Cohort | N (approx) | Data | Access | Use |
|--------|-----------|------|--------|-----|
| TCGA-SKCM | ~470 | RNA-seq, mutation, clinical | GDC Portal, open | Discovery, subtyping, survival |
| GSE78220 (Hugo et al., anti-PD-1 melanoma) | 28 | RNA-seq + response | GEO, open | Response validation |
| GSE91061 (Riaz et al., nivolumab melanoma) | 51 | RNA-seq + response | GEO, open | Response validation |
| PRJEB23709 / Liu-Gide (melanoma, anti-PD-1 ± anti-CTLA-4) | ~144 RNA | RNA-seq + response/survival | ENA / cBioPortal, open | Main immunotherapy test set |
| GSE120575 (melanoma scRNA-seq, anti-PD-1/CTLA-4) | 48 samples, 16k cells | scRNA-seq + response | GEO, open | Single-cell validation |
| GSE115978 (melanoma scRNA-seq ecosystem) | 31 tumors | scRNA-seq | GEO, open | Malignant resistance program validation |
| TCGA-LUAD + LUSC + GSE126044 + GSE135222 (NSCLC anti-PD-1) | ~1000 + 16 + 27 | RNA-seq + response/PFS | GDC + GEO, open | Cross-cancer validation |
| IMvigor210 (urothelial, atezolizumab) | ~350 | RNA-seq + response/survival | R package IMvigor210CoreBiologies | Cross-cancer validation |
| GDSC/CTRP (cell lines) | — | Expression + drug IC50 | DepMap / GSCA, open | Drug sensitivity prediction for resistant subtype |

Total: >2,000 bulk samples + ~2 scRNA cohorts. All accession numbers will be cited. No restricted (dbGaP) data.

---

## 4. Full analysis pipeline (every step re-runnable, code saved)

1. **Data acquisition & QC** — scripted download, gene harmonization, batch inspection, response harmonization (CR/PR vs SD/PD).
2. **TME deconvolution** — ESTIMATE, MCP-counter, xCell, EPIC, quanTIseq, CIBERSORT-LM22 (open signatures) → consensus immune/stromal fractions.
3. **Program scoring** — TLS signature, stromal/CAF/EMT, T-cell dysfunction/exclusion (TIDE-style open logic), IFN-γ, antigen presentation, wound-healing/angiogenesis scores via ssGSEA/GSVA.
4. **Ecosystem subtyping** — consensus clustering on ecosystem features → 3 subtypes; characterization + Sankey/cluster plots.
5. **Clinical association** — response rate per subtype (chi-square), PFS/OS Kaplan-Meier + Cox (uni/multivariate), ROC/AUC vs PD-L1/TMB baselines.
6. **Differential expression + enrichment** — DESeq2/limma, GO/KEGG/Hallmark GSEA, volcano/heatmaps.
7. **Genomic correlates** — TMB, MSI proxy, driver mutation enrichment per subtype (maftools).
8. **Crosstalk analysis** — ligand–receptor prioritization (bulk co-expression + scRNA CellChat-style) focusing on CAF→T exclusion axes (TGF-β, CXCL12–CXCR4, VEGF, IL-10).
9. **Explainable ML predictor** — LASSO → RF/XGBoost/logistic/SVM benchmark with 5-fold CV on discovery, locked model tested on 5 independent cohorts; SHAP values, calibration, decision-curve analysis, nomogram combining ecosystem score + TMB + stage.
10. **Single-cell validation** — map hub programs to malignant/CAF/T-cell compartments in GSE120575/GSE115978 (Seurat, UMAP, violin/feature plots).
11. **Therapeutic innovation** — drug sensitivity imputation (oncoPredict/GDSC) + candidate repurposing for Stromal-Excluded subtype; immune-checkpoint & TGF-β/CAF-targeting hypothesis.
12. **Reproducibility** — fixed seeds, sessionInfo, R/Python notebooks, GitHub-ready repo structure, supplementary tables.

---

## 5. Figure plan (all generated de novo — zero copied figures)

- Fig 1: Study design + ecosystem subtyping (consensus matrix, PCA/UMAP, heatmap of features).
- Fig 2: Subtypes vs response + survival (stacked bars, KM curves, forest plots).
- Fig 3: Biology of subtypes (GSEA ridges, volcano, pathway heatmap).
- Fig 4: Genomics + crosstalk (oncoplot per subtype, ligand–receptor network).
- Fig 5: ML predictor (ROC curves across 6 cohorts, SHAP beeswarm, calibration, nomogram).
- Fig 6: scRNA validation (UMAP, program mapping, cell-type violins).
- Fig 7: Therapeutic implications (drug sensitivity boxplots, proposed stratified strategy schematic — drawn original).

---

## 6. Originality / no-plagiarism guarantees

- All text written fresh for this study; references formatted; no text recycling.
- All figures rendered from our own code + public data — no copied images.
- All data reuse cited with accession numbers + primary citations (Hugo, Riaz, Liu, Jerby-Arnon, etc.).
- No fabricated patients, no simulated response labels — only real response/survival annotations from source studies.
- Methods transparency: code + version pins so editors/reviewers can reproduce.
- iThenticate-safe by construction: original narrative + proper quotation-free paraphrase + citations.

---

## 7. How we handle "experimental validation preferred"

The Collection does not mandate wet-lab validation for computational biomarker Articles, but reviewers often ask.
Mitigation built into design:
- 5+ independent immunotherapy cohorts (stronger than single-center qPCR).
- Orthogonal scRNA-seq validation in independent tumors.
- Benchmarking against established biomarkers (PD-L1 signature, TMB, TIDE, IMPRES) to prove added value.
- Clear Limitations paragraph + proposed prospective/IHC validation as future work.
- OPTIONAL (if user has lab access): low-cost IHC (CD20/CD8/αSMA TLS vs stroma) or qPCR of top 5 hub genes in local FFPE cohort — I will write the validation protocol, but paper is submission-ready without it.

---

## 8. Submission package I will deliver

- Main manuscript (.docx + .md): Title, Abstract, Intro, Methods, Results, Discussion, Limitations, Data/Code Availability, References, Figure legends.
- Figures: high-res PNG + PDF + source code.
- Supplementary: Tables S1–S8 (cohort summary, gene lists, model coefficients, performance metrics), Methods S1.
- Cover letter tailored to Collection editors (Bao, Gu, Xia) + suggested reviewers.
- Reproducibility repo: /code, /data_manifest, README, environment file.
- Plagiarism self-check report: data provenance table + originality notes.

Target length: 4,500–6,000 words, 60–80 refs, 6–7 figures — npj Precision Oncology standard.

---

## 9. Timeline (with 16 June 2027 deadline)

- Week 1–2: Data download + QC + deconvolution + subtyping (Fig 1–2 draft).
- Week 3–4: Enrichment + genomics + crosstalk (Fig 3–4).
- Week 5–6: ML + SHAP + nomogram + multi-cohort validation (Fig 5).
- Week 7–8: scRNA + drug prediction (Fig 6–7) + full manuscript draft v1.
- Week 9+: Internal QC, reference audit, cover letter, submission checklist.

We are ~9 months ahead — ample buffer for revisions and your authorship inputs.

---

## 10. Needed from you to start full pipeline

1. Confirm cancer focus: Melanoma-discovery + NSCLC/bladder validation (recommended) OR single-cancer alternative.
2. Author list + affiliations + corresponding author email (for title page).
3. Any local cohort for optional validation? (Yes/No — No is fine.)
4. Target journal only npj Precision Oncology, or backup journal too?

Once confirmed, I launch Step 1 (scripted data download + QC) in this workspace.
