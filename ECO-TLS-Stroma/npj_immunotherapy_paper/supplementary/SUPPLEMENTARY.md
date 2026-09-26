# Supplementary Information
## A TLS-versus-Stroma Ecosystem Score Predicts Checkpoint Inhibitor Benefit Across Melanoma, NSCLC and Bladder Cancer

### Supplementary Table S1. Gene signatures (fixed from literature before testing)
- **TLS (24):** CCL2, CCL3, CCL4, CCL5, CCL8, CCL18, CCL19, CCL21, CXCL9, CXCL10, CXCL11, CXCL13, MS4A1, CD79A, CD79B, LTB, CCR7, CCR6, CXCR5, SELL, CD1D, LAT, SKAP1, CETP
- **STROMA/CAF (30):** FAP, ACTA2, COL1A1, COL1A2, COL3A1, COL5A1, COL5A2, COL6A1, COL6A2, COL6A3, POSTN, VIM, PDPN, PDGFRA, PDGFRB, TGFB1, TGFB2, TGFB3, TGFBR2, MMP2, MMP11, MMP14, LOX, LOXL2, SPARC, FN1, VCAN, THY1, ITGB1, ZEB1
- **IFNG6 benchmark:** IDO1, CXCL10, CXCL9, HLA-DRA, STAT1, IFNG (Ayers et al. JCI 2017)
- **CD8EFF benchmark:** CD8A, CD8B, GZMA, GZMB, PRF1, IFNG
- **Checkpoints/axes:** CD274, PDCD1, CTLA4, LAG3, TIGIT, HAVCR2, PDCD1LG2, CXCL12, CXCR4, VEGFA, VEGFC, IL10, ENTPD1, NT5E, HIF1A, ARG1, MKI67
- **Lineage/QC:** PTPRC, EPCAM, KRT19, GAPDH, ACTB
- Full 85-gene universe with Entrez/Ensembl IDs: see `code/signatures.py` + NCBI gene_info (downloaded 2026-09-25).

### Supplementary Table S2. Cohort manifests and QC
See `results/tables/scores_{cohort}.csv` (per-sample TLS/STROMA/ECO/benchmark scores + subtype) and `data/processed/{cohort}_clin.csv` (response/survival annotations). Gene coverage: 85/85 in all cohorts except GSE135222 (81/85; missing genes excluded from means, not imputed).

### Supplementary Table S3. Full benchmark AUCs with 95% bootstrap CIs
See `results/tables/response_AUCs.csv` (6 metrics × 6 ICI cohorts, Mann–Whitney p, medians).

### Supplementary Table S4. ECO-high vs ECO-low response odds ratios
See `results/tables/ECO_highlow_OR.csv`.

### Supplementary Table S5. Subtype centroids and per-cohort response rates
See `results/tables/subtype_centroids.csv` and `results/tables/subtype_response.csv`. Pooled: TLS-high 38/135 (28.2%), Stromal 25/130 (19.2%), Immune-desert 48/174 (27.6%).

### Supplementary Table S6. TMB analysis (IMvigor210)
ECO–TMB Spearman ρ = 0.26, p = 0.0001. Response AUC (n = 214 with both): TMB 0.74 [0.67–0.82], ECO 0.58 [0.49–0.66], ECO+TMB 0.72 [0.64–0.80]. See `results/tables/IMvigor210_ECO_TMB.csv` and `TMB_stratified_ECO.csv`.

### Supplementary Methods S1. Endpoint harmonization
- GSE135222: durable clinical benefit (DCB) = PFS ≥ 180 days; non-benefit = progression < 180 days; censored < 180 days excluded (n = 0 excluded; all 27 classifiable).
- GSE207422: primary endpoint = major pathologic response (MPR incl. pCR) in pre-treatment biopsies; RECIST retained as secondary annotation.
- GSE91061: pre-treatment samples only; PR/CR = responder; SD/PD = non-responder; unknown excluded.
- IMvigor210: mRECIST CR/PR = responder; SD/PD = non-responder (iAtlas harmonization).

### Data provenance (anti-plagiarism record)
| File | Origin | Accessed | Records |
|---|---|---|---|
| TCGA SKCM/LUAD/BLCA RNA + clinical | Broad Firehose stddata__2016_01_28 | 2026-09-25 | 1,313 tumors |
| GSE78220 / GSE91061 / GSE126044 / GSE135222 / GSE207422 | NCBI GEO FTP | 2026-09-25 | 141 ICI tumors |
| IMvigor210 RNA + clinical | cBioPortal API (blca_iatlas_imvigor210_2017) | 2026-09-25 | 347 samples |
| Gene mappings | NCBI gene_info | 2026-09-25 | 85/85 mapped |
All figures rendered de novo from these sources with the released code; no image, table cell, or text passage was copied from any prior publication. Download/parse/score scripts: `code/01–09`.

### Supplementary Table S7. Pathway & regulator enrichment (Enrichr)
See `results/enhancement/enrich_*.csv` (TLS/STROMA × KEGG/Reactome/TRRUST/miRTarBase, top 15 terms each with adj. p).

### Supplementary Table S8. Networks
STRING high-confidence edges: `results/enhancement/string_edges.csv` (517 edges); TCGA co-expression matrix + strong edges + hub rankings: `coexpression_spearman.csv`, `coexpression_strong_edges.csv`, `coexpression_hubs.csv`, `string_hubs.csv`.

### Supplementary Table S9. Drug–gene interactions (DGIdb v4, 12 hubs, 333 rows)
See `results/enhancement/dgidb_drugs.csv` (gene, drug, approved flag, interaction type, sources).

### Supplementary Table S10. TIDE comparison + ECO×TMB quadrants
See `results/enhancement/tide_comparison.csv`, `eco_tmb_quadrants.csv`, `imv_TIDE_*.csv`, `imv_CLINICAL_BENEFIT.csv`, `imv_PROGRESSION.csv`.
