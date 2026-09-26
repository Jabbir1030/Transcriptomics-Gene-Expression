# ECO: TLS-versus-Stroma Ecosystem Score (npj Precision Oncology submission)

One number per tumor — immune organization **minus** stromal blockage — predicting
checkpoint-inhibitor benefit and survival across melanoma, NSCLC and bladder cancer.

- 9 public cohorts, ~1,830 tumors; 6 immunotherapy cohorts (n = 439, AUC 0.60–0.78)
- Survival validated in TCGA SKCM/LUAD/BLCA + IMvigor210 (HR 0.53–0.69)
- Mechanism: pathways, regulators, miR-29 brake, STRING networks, DGIdb drug map, TIDE beaten

## Contents (`npj_immunotherapy_paper/`)
- `code/` — pipeline 01 (download) → 14 (audit); run in numeric order
- `manuscript/` — `manuscript.md` + `manuscript.docx` (submittable), cover letter, backup journals
- `results/figures/` — Fig1–Fig9 (PNG + PDF, 300 DPI)
- `supplementary/` — Tables S1–S10 + Methods S1
- `data/` — raw + processed public data (re-downloadable via `code/01_download.py`)
- `BEGINNER_GUIDE.md` — step-by-step course explaining every analysis

## Quick start
```bash
cd npj_immunotherapy_paper
python3 code/07_scores.py       # core results (needs data/ present or re-downloaded)
python3 code/13_build_docx.py   # rebuild manuscript
```

## Status
Manuscript formatted per npj Precision Oncology Article guide (58/60 refs, unstructured
150-word abstract). Placeholders to fill before submission: author list, funding line, this repo URL.
