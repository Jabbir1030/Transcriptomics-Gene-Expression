---

## Discussion

This project began with an almost embarrassingly plain idea. Pick two published gene lists, subtract one average from the other, and check whether the difference foresees immunotherapy benefit. No training set tricks. No tuning. No black box. Across nine cohorts and some 1,830 tumors, the difference worked for response in six immunotherapy cohorts and for survival in five, pointing the same way without exception. Three lessons feel worth keeping.

First, the pair beats each partner. We kept watching single-sided scores stumble somewhere. TLS genes alone went mute in pre-treatment melanoma. Stromal genes alone faded in parts of lung. Their difference never pointed backward. That result stitches together two fields that seldom reference each other: TLS research showing organized immunity foresees benefit [17–22], and stroma research showing fibroblast and TGF-beta programs foresee resistance [23,24,34]. Both camps are right. Tumors escape checkpoint drugs through two exits (too little immune organization, or too much stromal barricade), and ECO simply weighs both exits on one scale. Transparency counts as a strength here rather than a weakness. Scores a clinician can grasp in a sentence, like the Immunoscore [31] and the T-cell–inflamed profile [30], travel to the bedside faster than 200-gene black boxes.

Second, stromal exclusion reads as the chief resistance axis. Look back at the loudest signal in this study: not any TLS gene, but the missing stromal program in responding melanomas. The subtype work concurred, with stromal tumors trailing the pooled response table and most individual cohorts. The desert finding seals the interpretation better than any statistic. Tumors nearly empty of immune signal answered as well as TLS-rich ones, which is precisely what exclusion biology predicts [23,28,41]. If this view holds, the highest-value future drugs in immunotherapy might not be the next T-cell checkpoint but the stromal normalizers: TGF-beta blockade [23,25], CXCR4 antagonists [39], anti-angiogenics, adenosine-axis drugs and FAP-directed agents. ECO-low status offers a pragmatic way to fill those trials with the right patients.

Third, ECO complements mutation testing instead of replacing it. ECO and TMB correlate at 0.26, which for practical purposes means they report separate bottlenecks. TMB asks whether targets exist. ECO asks whether T cells can arrive and organize. A tumor needs yes twice. In bladder cancer TMB won the duel outright, and we would rather print that than bury it. ECO's job there is backup duty when TMB is absent, plus one vote inside future multi-part classifiers fusing DNA, RNA and imaging [11,32]. Solo-hero narratives serve nobody.

Regulator analysis supplies the mechanism behind cross-cancer generalization, and it is what makes the result believable. TLS genes sit beneath NF-kB and RELA plus IRF1 and IRF3 command: old, tissue-agnostic software for building lymphoid tissue and sensing interferons [15,16]. Stromal genes sit beneath EMT factors (TWIST2, ETV4, SP1) with TGF-beta and SMAD signaling [25,26]: equally old software for wound repair and scarring, which cancers co-opt. Since neither circuit belongs to one tumor type, their balance legibly transfers from melanoma to lung to bladder. The miR-29 result adds a control lever. That microRNA family pins down collagens, SPARC, LOX enzymes and TGF-beta ligands, and its disappearance is a known trigger of fibrotic states [38]. An ECO-low tumor may thus be a miR-29–deficient tumor at the molecular level, which nominates miR-29 mimics as a sensitizer to test in models. The CXCL12 and CXCR4 hub pushes the same theme through another door, because fibroblast CXCL12 parks T cells in stroma and receptor blockade releases them [39].

What could clinicians plausibly do with this? We will stay concrete while refusing to overclaim. ECO demands nothing exotic: any bulk RNA readout covering the 54 program genes (RNA-seq, NanoString, targeted PCR) plus a within-cohort standardization any hospital lab can script. Three near-term roles look realistic. Trial stratification first, so stromal-heavy and TLS-heavy cases balance across arms instead of quietly confounding outcomes. Triage second, steering ECO-low and stromal patients toward stromal-combination studies (Table 4) rather than monotherapy likely to fail. Monitoring third, re-scoring repeat biopsies while ecosystems churn under therapeutic pressure [12,35], which no static baseline assay can follow. Fig. 9C sketches that future. It is a hypothesis awaiting a prospective trial with locked cutoffs, not a guideline, and we present it as such.

A few design choices harden the claims. Freezing gene lists from literature before touching outcomes shuts the overfitting route that sinks most signature papers. Validation crosses three cancers, five treatment contexts and four endpoints with zero cohort-specific adjustment. Benchmarks are the genuine field standards (IFN-gamma [29], TIDE [41], PD-L1, TMB), not weak foils. Mechanism layers rest on independent databases, so the biology narrative never leans on our cohorts alone. Data manifests, code and per-sample scores all ship for re-analysis.

Several limitations frame what comes next. Everything is retrospective and computational, which buys cohort-to-cohort validation but neither prospective proof nor bench experiments. Bulk RNA stays blind to TLS maturity (a loose cluster is no mature follicle [16,20]), blind to TLS address (intratumoral versus peritumoral placement flips meaning in some cancers), and blind to CAF flavors [27]. Spatial and single-cell work is not optional garnish; the field needs it. Outside IMvigor210 our immunotherapy cohorts run 16 to 49 patients, so single-cohort intervals gape and the argument leans on cross-cohort agreement. That is standard practice for these scarce annotated datasets, and standard does not mean ideal. Endpoints vary by cohort (RECIST, durable benefit, pathologic response), which argues for breadth while complicating pooling. Three cancers are covered; kidney, head and neck, MSI-high colorectal and the rest wait their turn. Finally the algorithm and drug pairings are proposals. The IDO1 story [40] stands as the permanent warning that elegant rationale plus elegant biomarkers can still detonate in phase III.

The natural next steps order themselves. A prospective study with a locked ECO cutoff comes first, preferably pan-cancer and multi-arm, asking whether ECO-low patients gain from stromal combinations. Paired RNA with multiplex imaging should tie ECO states to TLS maturity and CAF geography. Single-cell work should assign hub genes to cell types. ECO deserves testing in further cancers and regimens, including chemo-immunotherapy and neoadjuvant settings [5]. And the miR-29 and CXCR4 hypotheses deserve bench testing in exclusion models [38,39]. Every item is doable with tools already on shelves.

Checkpoint therapy dragged oncology from targeting the cancer cell to treating its ecosystem [13,14]. ECO is a modest, usable step along that path: a single number from two programs, biology you can explain, and a straight line from measurement to the next experiment worth running.

---

## Methods

### Data sources

Public sources only. No new patients and no new sequencing. TCGA RNA (RSEM-normalized) with clinical follow-up for melanoma (SKCM), lung adenocarcinoma (LUAD) and bladder cancer (BLCA) came from the Broad GDAC Firehose frozen release of January 2016 [33], which standardizes data from The Cancer Genome Atlas [42]. Immunotherapy cohorts came from GEO [43]. GSE78220 gave melanoma anti-PD-1 FPKM with RECIST response [34]. GSE91061 gave melanoma nivolumab data; we kept pre-treatment samples and scored CR/PR against SD/PD [35]. GSE126044 gave NSCLC anti-PD-1 counts with responder tags. GSE135222 gave NSCLC anti-PD-1/PD-L1 data with progression-free survival. GSE207422 gave neoadjuvant anti-PD-1 plus chemotherapy NSCLC log2TPM with major pathologic response. Bladder atezolizumab RNA, modified RECIST response, survival, TMB and precomputed TIDE calls arrived via the cBioPortal interface [36] for the iAtlas-harmonized IMvigor210 study [23]. We inspected GSE176307 and set it aside: its deposited response columns are shifted across samples, which makes response mapping untrustworthy. Future users of that file deserve the warning. Table 1 lists final sample sizes after quality checks.

### Ethics and consent

All data analyzed here were previously published, de-identified and publicly available; no new participants were recruited and no new sequencing was performed. Ethics approval and participant consent for the original studies were obtained by the original investigators through their own institutional review boards [23,33–36,42,43]. Because no identifiable human material or data were accessed, no additional consent was required for this secondary analysis, which was conducted in accordance with the principles of the Declaration of Helsinki.

### Gene lists and the ECO score

The design fits in two sentences. Two published gene teams vote against each other, and ECO is the margin. The TLS roster (24 genes: CCL2, CCL3, CCL4, CCL5, CCL8, CCL18, CCL19, CCL21, CXCL9, CXCL10, CXCL11, CXCL13, MS4A1, CD79A, CD79B, LTB, CCR7, CCR6, CXCR5, SELL, CD1D, LAT, SKAP1, CETP) merges the classic 12-chemokine TLS signal [44] with B-cell, Tfh and TLS-imprint genes from the defining TLS papers [15–19]. The stromal roster (30 genes: FAP, ACTA2, COL1A1, COL1A2, COL3A1, COL5A1, COL5A2, COL6A1, COL6A2, COL6A3, POSTN, VIM, PDPN, PDGFRA, PDGFRB, TGFB1, TGFB2, TGFB3, TGFBR2, MMP2, MMP11, MMP14, LOX, LOXL2, SPARC, FN1, VCAN, THY1, ITGB1, ZEB1) spans TGF-beta and CAF exclusion [23,24], collagen and matrix remodeling, plus mesenchymal resistance programs [34]. We committed both rosters to writing before seeing outcomes. Comparators were the IFN-gamma 6-gene set [29], a CD8-effector set (CD8A, CD8B, GZMA, GZMB, PRF1, IFNG), and CD274 transcript solo. Identifier mapping used NCBI gene_info. Per cohort we standardized log2 expression per gene into z-scores. Each program score is the mean z across its detected genes. **ECO = TLS mean minus stroma mean.** Coverage hit 85 of 85 universe genes everywhere except GSE135222 (81 of 85); absent genes dropped out of means rather than suffering imputation.

### Response endpoints

Trials define success differently, so we translated each to a binary label honoring its original publication (Supplementary Methods S1). GSE135222 lacked RECIST tags, hence durable clinical benefit: progression-free survival reaching 180 days counts as benefit, progression before that as non-benefit, and censoring before 180 days as unknowable (excluded). GSE207422 used major pathologic response including complete pathologic response in pre-treatment biopsies. GSE91061 used pre-treatment samples (CR/PR versus SD/PD, unknowns out). IMvigor210 used modified RECIST (CR/PR versus SD/PD).

### Ecosystem subtypes

For natural groupings we clustered the 1,313 pooled TCGA tumors on TLS/stroma coordinates with k-means (k = 3) [45]. Names follow profiles: TLS-high (rich TLS, calm stroma), stromal (rich stroma), immune-desert (scant both). Every immunotherapy tumor then joined its nearest fixed TCGA centroid. The immunotherapy data never touched the fitting.

### Pathways, regulators, networks, drugs and TIDE

Program biology got annotated through gene-set enrichment [46,47] spanning KEGG 2021, Reactome 2022, TRRUST transcription factors [48] and miRTarBase microRNAs [49], with standard multiple-testing correction. STRING v12 supplied protein links across the 85 genes at high confidence (score at least 0.7) [50]. Co-expression is Spearman correlation over pooled TCGA tumors (n = 1,363 samples with full data), with hubs ranked by network degree [51]. Twelve hubs (FAP, TGFB1, POSTN, VEGFA, CXCL12, CXCR4, CD274, CTLA4, IDO1, ENTPD1, NT5E, MMP14) went to the Drug–Gene Interaction Database v4 GraphQL service [52]; approval flags are DGIdb-curated, and Table 4 spotlights agents with direct mechanistic fit. TIDE responder tags were the precomputed iAtlas annotations, never recomputed by us; single-point AUC equals (sensitivity + specificity)/2 [41]. Marker-based interpretation follows established deconvolution practice [53,54].

### Statistics and reproducibility

Response prediction used ROC curves with AUC, 95% percentile bootstrap intervals from 2,000 resamples, and two-sided Mann–Whitney tests. We chose Mann–Whitney because expression scores cannot be assumed normal, particularly in cohorts of 16 to 49 tumors. Pooled immunotherapy analysis standardized ECO within cohorts first, then merged (n = 439). High/low cuts always sat at the cohort median; odds ratios used Fisher's exact test. Survival analysis used Kaplan–Meier curves [55], log-rank tests and Cox models [56] contrasting ECO-high with ECO-low, plus multivariate runs adjusting for age, sex and stage (TCGA) or TMB status (IMvigor210). ECO–TMB linkage used Spearman's rank test; their joint model was ordinary logistic regression. Significance sat at two-sided p below 0.05, and we report exact p values throughout rather than thresholds alone. Computing ran on Python 3.13 (pandas, numpy, scipy, scikit-learn 1.6.1 [45], lifelines 0.30.3 [57], matplotlib, networkx [51]). Since ECO was one pre-specified score under validation rather than a discovery trawl, cohorts carry no multiple-testing penalty, a conventional validation setup consistent with TRIPOD-style reporting [58]. Download, parsing, scoring, statistics and figures exist as seeded scripts (Code Availability) with pinned package versions. Supplementary tables hold cohort manifests, gene coverage, full benchmarks, subtype centroids and TMB-stratified performance.

### Use of AI assistance

Large language models assisted with language editing and code drafting during preparation of this manuscript. All data analyses were executed and checked by the authors, every reported number traces to the released code and public datasets, and the authors take full responsibility for the content. This disclosure follows Springer Nature editorial policy on AI use.

---

## Data Availability

No new data were generated. TCGA data: Broad GDAC Firehose stddata__2016_01_28 (SKCM/LUAD/BLCA). GEO: GSE78220, GSE91061, GSE126044, GSE135222, GSE207422. IMvigor210: cBioPortal study `blca_iatlas_imvigor210_2017`. Processed per-cohort expression matrices (85-gene universe) and clinical annotations used in this study are provided as Supplementary Data.

## Code Availability

Analysis code (Python) reproducing all scores, statistics, tables and figures is available at [GitHub repository URL to be inserted upon acceptance] and as a Supplementary Software archive. Gene signatures are listed in Methods and Supplementary Table S1.

## Acknowledgements

We thank the patients, investigators and consortia behind TCGA, the GEO-deposited immunotherapy cohorts and the IMvigor210 trial for making their data public, and the curators of Firehose, GEO and cBioPortal for keeping it usable. We are grateful to the open-source developers of the scientific Python ecosystem. This study received no funding. [Confirm or replace with funder and grant details before submission.]

## Author Contributions

[To be completed: e.g., First Author — conceptualization, analysis, writing; Second Author — validation, review; Corresponding Author — supervision, writing. Use initials for each author's contribution.]

## Competing Interests

All authors declare no financial or non-financial competing interests.

## References

1. Ribas, A. & Wolchok, J. D. Cancer immunotherapy using checkpoint blockade. *Science* **359**, 1350–1355 (2018).
2. Hodi, F. S. et al. Improved survival with ipilimumab in patients with metastatic melanoma. *N. Engl. J. Med.* **363**, 711–723 (2010).
3. Robert, C. et al. Nivolumab in previously untreated melanoma without BRAF mutation. *N. Engl. J. Med.* **372**, 320–330 (2015).
4. Reck, M. et al. Pembrolizumab versus chemotherapy for PD-L1–positive non-small-cell lung cancer. *N. Engl. J. Med.* **375**, 1823–1833 (2016).
5. Gandhi, L. et al. Pembrolizumab plus chemotherapy in metastatic non-small-cell lung cancer. *N. Engl. J. Med.* **378**, 2078–2092 (2018).
6. Rosenberg, J. E. et al. Atezolizumab in patients with locally advanced and metastatic urothelial carcinoma who have progressed following treatment with platinum-based chemotherapy. *Lancet* **387**, 1909–1920 (2016).
7. Davis, A. A. & Patel, V. G. The role of PD-L1 expression as a predictive biomarker: an analysis of published clinical trials. *J. Immunother. Cancer* **7**, 278 (2019).
8. Yarchoan, M., Hopkins, A. & Jaffee, E. M. Tumor mutational burden and response rate to PD-1 inhibition. *N. Engl. J. Med.* **377**, 2500–2501 (2017).
9. Samstein, R. M. et al. Tumor mutational load predicts survival after immunotherapy across multiple cancer types. *Nat. Genet.* **51**, 202–206 (2019).
10. Le, D. T. et al. Mismatch repair deficiency predicts response of solid tumors to PD-1 blockade. *Science* **357**, 409–413 (2017).
11. Litchfield, K. et al. Meta-analysis of tumor- and T cell–intrinsic mechanisms of sensitization to checkpoint inhibition. *Cell* **184**, 596–614 (2021).
12. Sharma, P., Hu-Lieskovan, S., Wargo, J. A. & Ribas, A. Primary, adaptive, and acquired resistance to cancer immunotherapy. *Cell* **168**, 707–723 (2017).
13. Chen, D. S. & Mellman, I. Elements of cancer immunity and the cancer–immune set point. *Nature* **541**, 321–330 (2017).
14. Galon, J. & Bruni, D. Approaches to treat immune hot, altered and cold tumours with combination immunotherapies. *Nat. Rev. Drug Discov.* **18**, 197–218 (2019).
15. Sautès-Fridman, C. et al. Tertiary lymphoid structures in the era of cancer immunotherapy. *Nat. Rev. Cancer* **19**, 307–325 (2019).
16. Schumacher, T. N. & Thommen, D. S. Tertiary lymphoid structures in cancer. *Science* **375**, eabf9419 (2022).
17. Cabrita, R. et al. Tertiary lymphoid structures improve immunotherapy and survival in melanoma. *Nature* **577**, 561–565 (2020).
18. Helmink, B. A. et al. B cells and tertiary lymphoid structures promote immunotherapy response. *Nature* **577**, 549–555 (2020).
19. Petitprez, F. et al. B cells are associated with survival and immunotherapy response in sarcoma. *Nature* **577**, 556–560 (2020).
20. Vanhersecke, L. et al. Mature tertiary lymphoid structures predict immune checkpoint inhibitor efficacy in solid tumors independently of PD-L1 expression. *Nat. Cancer* **2**, 794–802 (2021).
21. Patil, N. S. et al. Intratumoral plasma cells predict outcomes to PD-L1 blockade in non–small cell lung cancer. *Cancer Cell* **40**, 289–300 (2022).
22. Li, J., Qi, W., Ma, L. et al. The predictive value of intratumoral tertiary lymphoid structures on the response to immunotherapy in cancer patients: a systematic review and meta-analysis. *BMC Cancer* **25**, 1935 (2025).
23. Mariathasan, S. et al. TGFβ attenuates tumour response to PD-L1 blockade by contributing to exclusion of T cells. *Nature* **554**, 544–548 (2018).
24. Chakravarthy, A. et al. TGF-β-associated extracellular matrix genes link cancer-associated fibroblasts to immune evasion and immunotherapy failure. *Nat. Commun.* **9**, 4692 (2018).
25. Batlle, E. & Massagué, J. Transforming growth factor-β signaling in immunity and cancer. *Immunity* **50**, 924–940 (2019).
26. Sahai, E. et al. A framework for advancing our understanding of cancer-associated fibroblasts. *Nat. Rev. Cancer* **20**, 174–186 (2020).
27. Öhlund, D. et al. Distinct populations of inflammatory fibroblasts and myofibroblasts in pancreatic cancer. *J. Exp. Med.* **214**, 579–596 (2017).
28. Turley, S. J., Cremasco, V. & Astarita, J. L. Immunological hallmarks of stromal cells in the tumour microenvironment. *Nat. Rev. Immunol.* **15**, 669–682 (2015).
29. Ayers, M. et al. IFN-γ–related mRNA profile predicts clinical response to PD-1 blockade. *J. Clin. Invest.* **127**, 2930–2940 (2017).
30. Ott, P. A. et al. T-cell–inflamed gene-expression profile, programmed death ligand 1 expression, and tumor mutational burden predict efficacy in patients treated with pembrolizumab across 20 cancers: KEYNOTE-028. *J. Clin. Oncol.* **37**, 318–327 (2019).
31. Pagès, F. et al. International validation of the consensus Immunoscore for the classification of colon cancer. *Lancet* **391**, 2128–2139 (2018).
32. Bao, R., Gu, Q., Xia, R. et al. Next-generation precision immunotherapy: from tumor ecosystems to therapeutic innovation (Collection). *npj Precis. Oncol.* (2026).
33. Broad Institute TCGA Genome Data Analysis Center. Firehose stddata__2016_01_28. Broad Institute (2016).
34. Hugo, W. et al. Genomic and transcriptomic features of response to anti-PD-1 therapy in metastatic melanoma. *Cell* **166**, 35–46 (2016).
35. Riaz, N. et al. Tumor and microenvironment evolution during immunotherapy with nivolumab. *Cell* **171**, 934–949 (2017).
36. Cerami, E. et al. The cBio cancer genomics portal: an open platform for exploring multidimensional cancer genomics data. *Cancer Discov.* **2**, 401–404 (2012).
37. Wherry, E. J. & Kurachi, M. Molecular and cellular insights into T cell exhaustion. *Nat. Rev. Immunol.* **15**, 486–492 (2015).
38. Rupaimoole, R. & Slack, F. J. MicroRNA therapeutics: towards a new era for the management of cancer and other diseases. *Nat. Rev. Drug Discov.* **16**, 203–222 (2017).
39. Feig, C. et al. Targeting CXCL12 from FAP-expressing carcinoma-associated fibroblasts synergizes with anti–PD-L1 immunotherapy in pancreatic cancer. *Proc. Natl Acad. Sci. USA* **110**, 20711–20716 (2013).
40. Long, G. V. et al. Epacadostat plus pembrolizumab versus placebo plus pembrolizumab in patients with unresectable or metastatic melanoma (ECHO-301/KEYNOTE-252). *Lancet Oncol.* **20**, 1083–1097 (2019).
41. Jiang, P. et al. Signatures of T cell dysfunction and exclusion predict cancer immunotherapy response. *Nat. Med.* **24**, 1550–1558 (2018).
42. Weinstein, J. N. et al. The Cancer Genome Atlas Pan-Cancer analysis project. *Nat. Genet.* **45**, 1113–1120 (2013).
43. Barrett, T. et al. NCBI GEO: archive for functional genomics data sets—update. *Nucleic Acids Res.* **41**, D991–D995 (2013).
44. Messina, J. L. et al. 12-chemokine gene expression signature of inflammation: a new biomarker for cancer prognosis and immunotherapy. *Sci. Rep.* **2**, 765 (2012).
45. Pedregosa, F. et al. Scikit-learn: machine learning in Python. *J. Mach. Learn. Res.* **12**, 2825–2830 (2011).
46. Subramanian, A. et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *Proc. Natl Acad. Sci. USA* **102**, 15545–15550 (2005).
47. Kuleshov, M. V. et al. Enrichr: a comprehensive gene set enrichment analysis web server 2016 update. *Nucleic Acids Res.* **44**, W90–W97 (2016).
48. Han, H. et al. TRRUST v2: an expanded reference database of human and mouse transcriptional regulatory interactions. *Nucleic Acids Res.* **46**, D380–D386 (2018).
49. Chou, C. H. et al. miRTarBase update 2018: a resource for experimentally validated microRNA–target interactions. *Nucleic Acids Res.* **46**, D296–D302 (2018).
50. Szklarczyk, D. et al. The STRING database in 2023: protein–protein association networks and functional enrichment analyses. *Nucleic Acids Res.* **51**, D638–D646 (2023).
51. Hagberg, A. A., Schult, D. A. & Swart, P. J. Exploring network structure, dynamics, and function using NetworkX. In *Proc. 7th Python in Science Conference* 11–15 (2008).
52. Freshour, S. L. et al. Integration of the Drug–Gene Interaction database (DGIdb 4.0) with open crowdsource efforts. *Nucleic Acids Res.* **49**, D1144–D1151 (2021).
53. Newman, A. M. et al. Robust enumeration of cell subsets from tissue expression profiles. *Nat. Methods* **12**, 453–457 (2015).
54. Yoshihara, K. et al. Inferring tumour purity and stromal and immune cell admixture from expression data. *Nat. Commun.* **4**, 2612 (2013).
55. Kaplan, E. L. & Meier, P. Nonparametric estimation from incomplete observations. *J. Am. Stat. Assoc.* **53**, 457–481 (1958).
56. Cox, D. R. Regression models and life-tables. *J. R. Stat. Soc. Series B* **34**, 187–220 (1972).
57. Davidson-Pilon, C. lifelines: survival analysis in Python. *J. Open Source Softw.* **4**, 1317 (2019).
58. Collins, G. S. et al. Transparent reporting of a multivariable prediction model for individual prognosis or diagnosis (TRIPOD): the TRIPOD statement. *BMJ* **350**, g7594 (2015).

## Figure Legends

**Figure 1. Study design.** ECO is the TLS program (24 genes) minus the stromal/CAF program (30 genes), standardized within each cohort. We charted the overview in TCGA (n = 1,313) and validated response and survival in six immunotherapy cohorts (melanoma, NSCLC, bladder; n = 439 with response). Benchmarking, subtyping, mechanism work and drug mapping followed. Public data only; validation runs strictly from cohort to cohort.

**Figure 2. ECO and response in six cohorts.** Top row: ECO scores in non-responders versus responders, with AUC and Mann–Whitney p-values. Bottom row: ROC curves for ECO against the IFN-gamma and PD-L1 benchmarks.

**Figure 3. Benchmarking.** (A) Response AUC for each metric in each cohort. (B) ECO AUCs with 95% bootstrap confidence intervals — above 0.5 in all six cohorts, always in the same direction.

**Figure 4. ECO and survival.** Kaplan–Meier curves for overall survival (three TCGA cohorts plus IMvigor210) and progression-free survival (GSE135222), split at each cohort's median ECO. Hazard ratios with 95% confidence intervals and Cox p-values are shown.

**Figure 5. Ecosystem states.** (A) Pooled TCGA tumors plotted by TLS versus stroma activity, grouped into three states. (B) Response rates by state after placing each immunotherapy tumor into the nearest fixed TCGA state — stromal tumors respond least.

**Figure 6. ECO, TMB and the inflamed tumor.** (A) ECO correlates only weakly with TMB in IMvigor210. (B) Response AUCs for TMB, ECO and the two combined. (C) Average gene activity by ECO group: ECO-high tumors look inflamed and stroma-poor.

**Figure 7. How the two programs are wired.** (A) Pathway analysis: TLS genes belong to chemokine, cytokine and Toll-like receptor signaling; stromal genes belong to matrix, adhesion and collagen handling. (B) Regulators: NF-kB and IRF factors drive TLS genes, while TWIST2, ETV4 and SP1 drive stromal genes. (C) Protein interaction network and (D) TCGA co-expression network both split into a TLS module (teal) and a stromal module (rust), joined by checkpoint genes.

**Figure 8. Control switches and drug targets.** (A) The regulation story in one picture: interferon and NF-kB signaling switch TLS chemokines on; EMT factors and TGF-beta switch stromal genes on, while the miR-29 family holds 15 stromal genes down — losing miR-29 opens the door to exclusion. (B) Drug–gene network linking 12 hub genes to approved and experimental drugs (DGIdb).

**Figure 9. Toward the clinic.** (A) ECO beats precomputed TIDE calls in IMvigor210. (B) Response rates by combined ECO and TMB status. (C) A proposed ECO-guided treatment algorithm — ECO-high patients toward checkpoint blockade, ECO-low/stromal patients toward trials that add a stromal-axis drug, with ECO re-checked at progression. A hypothesis for prospective testing, not a guideline.

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

STROMA AUC below 0.5 means an inverse link (low stroma goes with response), as expected.

### Table 3. Survival by ECO-high vs ECO-low

| Cohort | Endpoint | n (events) | HR [95% CI] | p | Adjusted HR | Adj. p |
|---|---|---|---|---|---|---|
| TCGA-SKCM | OS | 440 (153) | 0.53 [0.38–0.73] | 0.0001 | 0.60 | 0.002 |
| TCGA-LUAD | OS | 477 (121) | 0.63 [0.44–0.90] | 0.011 | 0.69 | 0.052 |
| TCGA-BLCA | OS | 394 (107) | 0.64 [0.44–0.95] | 0.025 | 0.76 | 0.17 |
| IMvigor210 | OS | 326 (213) | 0.69 [0.52–0.90] | 0.006 | 0.84 | 0.28 |
| GSE135222 | PFS | 27 (21) | 0.53 [0.22–1.27] | 0.15 | — | — |

### Table 4. Treatment ideas: ECO-network hubs, matching drug classes, and trial concepts for ECO-low/stromal patients

| Hub gene(s) | Program | Drug class (examples) | Proposed concept |
|---|---|---|---|
| TGFB1, TGFBR2 | Stroma | Pirfenidone; PD-L1×TGF-β bispecifics | Checkpoint + TGF-β blockade in TGFB/COL-high tumors |
| CXCL12, CXCR4 | Stroma–immune bridge | Plerixafor, motixafortide, mavorixafor | Checkpoint + CXCR4 blocker in CXCL12-high exclusion |
| VEGFA | Stroma/blood vessels | Bevacizumab, aflibercept, pazopanib | Checkpoint + VEGF blockade in angiogenic stroma |
| FAP | CAF | FAP-directed (CAR, IL2v, radioligands; experimental) | Checkpoint + FAP targeting in FAP-heavy CAF |
| ENTPD1, NT5E (CD39, CD73) | Adenosine | CD39/CD73 blockers (experimental) | Checkpoint + adenosine blockade in NT5E-high tumors |
| IDO1 | Metabolic | Epacadostat (caution: failed phase III) | Only in biomarker-selected re-testing |
| POSTN, MMP14, collagens | Matrix | Antifibrotics; miR-29 mimics (preclinical) | Stromal softening + checkpoint |
| CD274, CTLA4, LAG3 | TLS/effector | Approved checkpoint drugs | ECO-high: checkpoint therapy per label |
