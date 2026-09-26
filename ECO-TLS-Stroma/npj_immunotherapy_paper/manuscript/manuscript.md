# A TLS-versus-stroma ecosystem score predicts checkpoint inhibitor benefit across three cancers

**Authors:** [First Author]¹, [Second Author]², [Corresponding Author]¹\*
**Affiliations:** ¹[Department, Institution, City, Country]; ²[Department, Institution, City, Country]
**\*Corresponding author:** [Name, Email]
**Article type:** Article — Collection: Next-Generation Precision Immunotherapy: From Tumor Ecosystems to Therapeutic Innovation, *npj Precision Oncology*

> **Author note:** Replace bracketed placeholders with the final author list, affiliations and corresponding details before submission. All analyses, figures and tables in this manuscript were generated de novo from public data (see Data Availability); no text, figure or dataset was copied from any prior publication.

---

## Abstract

Checkpoint inhibitors benefit only a minority of patients with melanoma, lung or bladder cancer, and single-axis biomarkers incompletely predict response. We asked whether the balance between organized immunity and stromal exclusion governs benefit. We defined ECO, the TLS-versus-stroma Ecosystem Score, as mean expression of 24 tertiary-lymphoid-structure genes minus mean expression of 30 stromal and fibroblast genes, with both lists fixed from the literature before any outcome analysis. Across six immunotherapy cohorts (n = 439), responders scored higher than non-responders in every cohort (AUC 0.60 to 0.78; pooled 0.62). ECO-high patients lived longer in three TCGA cohorts (hazard ratios 0.53 to 0.64) and in IMvigor210 (0.69). Stromal-subtype tumors responded worst (19.2%). ECO correlated weakly with tumor mutational burden, outperformed TIDE in IMvigor210, and resolved to druggable stromal hubs. One pre-specified balance between immune organization and stromal exclusion predicts immunotherapy benefit and survival across three cancers.

---

## Introduction

Antibodies against CTLA-4, PD-1 and PD-L1 have produced remissions lasting years in melanoma [1–3], lung cancer [4,5], bladder cancer [6] and beyond [1]. Yet most patients still draw no lasting benefit. Response rates to PD-1 blockade alone hover around 20 to 45 percent depending on tumor type, and predicting who lands on which side remains unreliable. That uncertainty has real costs. Non-responders suffer side effects, lose months, and miss their window for something else.

Clinic-ready biomarkers help at the margins but no further. PD-L1 staining is the everyday test, and its predictive power is modest at best, shifting with each antibody, cutoff and cancer type [7]. Tumor mutational burden predicts response in several settings [8,9], as does mismatch-repair deficiency [10]. Still, TMB-high tumors often resist and TMB-low tumors often respond. Clonal architecture, antigen quality and the immune neighborhood around the tumor weigh just as heavily as raw mutation counts, and pan-cancer synthesis keeps confirming that no single axis suffices [11]. Research into resistance keeps adding fragments, yet no DNA-only test sees the full picture [12].

A decade of work has pushed the field toward a wider view: the tumor as an ecosystem rather than a lump of malignant cells [13]. Killer T cells, B cells, macrophages, dendritic cells, vessels and fibroblasts share one crowded tissue and negotiate constantly through chemokines, cytokines and contact. Tumors run hot, altered or cold, and treatment logic increasingly follows that geography rather than mutation lists alone [14]. Within that crowded tissue, two features now tower over the rest.

The first is the tertiary lymphoid structure. TLS are small lymph-node-like hubs that assemble inside or beside tumors. A mature one carries B-cell follicles, T-cell zones and specialized dendritic cells, and it functions as a local training ground where immune cells learn the tumor [15,16]. Then came 2020, when three landmark papers tied B cells and TLS to checkpoint response in melanoma and sarcoma [17–19]. Follow-up studies carried the finding into lung cancer and broader solid tumors [20,21]. A 2025 meta-analysis across 15 studies and 1,307 patients put the effect near fourfold higher response in TLS-high disease (odds ratio 4.21) [22]. Organized local immunity, in short, matters.

The second feature works against treatment: stroma. Cancer-associated fibroblasts spin collagen and matrix, pour out TGF-beta and growth factors, and can wall T cells out of the tumor core [23,24]. These fibroblasts are not one cell type; inflammatory and muscle-like subtypes do different jobs, and the field has begun to systematize them [25–27]. Stromal programs already predict bad outcomes in several cancers. In bladder cancer, stromal TGF-beta signaling marks tumors that anti-PD-L1 cannot enter [23]. Fibroblasts and immune cells talk nonstop [28], and when fibroblasts run the conversation, checkpoint drugs tend to fail.

Something is missing between these two literatures: they rarely meet inside a single test. Nearly every published RNA predictor scores immunity alone (interferon signatures [29], T-cell–inflamed profiles [30], cytolytic scores) or trains an opaque model within one cancer. The Immunoscore proved that a readable readout can reach the clinic [31]. But a joint TLS-versus-stroma balance, fixed before testing and validated across cancers without re-tuning, did not exist. Real tumors carry both programs at once. Treatment decisions depend on their net sum. So why does nobody measure the sum?

Here, we define that sum. We call it ECO, the TLS-versus-stroma Ecosystem Score: mean activity of 24 TLS genes minus mean activity of 30 stromal and fibroblast genes, read straight from bulk tumor RNA. No machine learning. No fitted weights. No cutoff hunting. Both gene lists came from the published literature and were frozen before we examined a single outcome, so everything below is genuine validation. This work was developed for the Collection on next-generation precision immunotherapy [32].

We test ECO in nine public cohorts and roughly 1,830 tumors from melanoma, lung and bladder cancer. Response to anti-PD-1, anti-PD-L1 and chemo-immunotherapy. Overall and progression-free survival. Benchmarks against the standards. Then we dig into mechanism: pathways, transcription factors, microRNAs, gene networks, matching drugs, and a direct contest with an established predictor. We close with a concrete proposal for putting ECO to work. Public data throughout, and code released so anyone can redo every step.

---

## Results

### Study design and cohorts

Fig. 1 lays out the plan. For the overview and prognosis we used three TCGA cohorts: skin melanoma (SKCM, 440 tumors with survival), lung adenocarcinoma (LUAD, 477) and bladder urothelial carcinoma (BLCA, 394), all from the Broad GDAC Firehose frozen release [33]. That is 1,313 tumors (Table 1). For response prediction we collected six independent pre-treatment immunotherapy cohorts with RNA plus response labels. Two melanoma cohorts on anti-PD-1 drugs (GSE78220 [34], n = 25; GSE91061 [35], pre-treatment n = 49). Three NSCLC cohorts (GSE126044, n = 16; GSE135222, n = 27 with progression-free survival; GSE207422, pre-treatment neoadjuvant n = 24 with pathologic response). And the large bladder atezolizumab trial IMvigor210 (n = 298 with response, n = 326 with survival, TMB in a subset) [23], accessed via cBioPortal [36]. Altogether 439 treated tumors with response labels, spanning three cancers and five treatment settings. We generated no new patient data. Each validation step runs from one independent cohort to the next (Methods).

### ECO separates responders from non-responders in six immunotherapy cohorts

The headline finding needs little adornment (Fig. 2, Table 2). Responders carried higher ECO scores than non-responders in every one of the six immunotherapy cohorts. AUCs ran 0.71 in GSE78220, 0.60 in GSE91061, 0.78 in GSE126044, 0.64 in GSE135222, 0.75 in GSE207422 and 0.60 in IMvigor210. Pooling all 439 tumors after standardizing ECO inside each cohort gave a combined AUC of 0.62. Formal significance arrived in two cohorts (GSE207422, p = 0.049; IMvigor210, p = 0.014) and approached it in two more (GSE78220, p = 0.087; GSE126044, p = 0.090). The last two small cohorts pointed the same direction without the sample size to prove it. For us the consistency carries the argument: six cohorts, three cancers, five regimens, zero reversals. Split each cohort at its own median and ECO-high tumors responded more often in five of six, with odds ratios from 1.36 to 4.50 (Supplementary Table S4).

### Benchmarking against single-program and established signatures

Benchmarking explains why subtraction beats single-sided scoring (Fig. 3, Table 2). Take pre-treatment melanoma (GSE78220). TLS genes alone were dead weight there (AUC 0.45). Stromal genes told the opposite story and screamed it: low stromal activity marked responders with an AUC of 0.21 (far under 0.5, so strongly inverted, p = 0.016). Anyone who knows that cohort will recognize the IPRES stromal resistance program from its original publication [34]. Elsewhere the roles reversed and TLS did the heavy lifting. Across the full set, ECO topped TLS alone in four cohorts and PD-L1 transcript in all six, while running even with the IFN-gamma [29] and CD8-effector standards. Our reading is plain. Tumors dodge immunotherapy two ways: too little immune organization, or too much stromal fortification. A test watching one door misses half the escapes. ECO watches both.

### ECO stratifies overall and progression-free survival

Response matters, but patients live or die by survival. At each cohort's median ECO, the upper half outlived the lower half everywhere testable (Fig. 4, Table 3). TCGA melanoma showed a hazard ratio of 0.53 (95% CI 0.38 to 0.73, p = 0.0001), which means roughly half the death rate at any moment. Lung came in at 0.63 (0.44 to 0.90, p = 0.011) and bladder at 0.64 (0.44 to 0.95, p = 0.025). Treated patients behaved the same: atezolizumab bladder cancer gave 0.69 (0.52 to 0.90, p = 0.006), and anti-PD-1 lung cancer trended identically for progression-free survival (HR 0.53, p = 0.15, n = 27). Because age and stage can masquerade as biology, we adjusted for them. ECO stayed significant in melanoma with age, sex and stage in the model (adjusted HR 0.60, p = 0.002) and borderline in lung (adjusted HR 0.69, p = 0.052).

### Ecosystem subtypes and response

Tumors refuse to spread smoothly across the TLS–stroma map. They clump. Clustering the 1,313 TCGA tumors on the two coordinates yielded three natural states (Fig. 5A). TLS-high tumors pair organized immunity with quiet stroma. Stromal tumors run busy fibroblasts and heavy matrix. Immune-desert tumors show little of either. We then dropped each immunotherapy tumor into its nearest fixed TCGA state. No re-clustering, no outcome peeking. Pooled response rates separated cleanly (Fig. 5B). Stromal tumors answered least (25 of 130, 19.2%). TLS-high (38 of 135, 28.1%) and desert (49 of 174, 28.2%) answered about equally. Cohort by cohort, TLS-high led in four of six and stromal trailed in four of six.

The desert result startles newcomers, so it deserves a straight explanation. Tumors with barely any immune signal responded as well as TLS-rich ones. Under an exclusion model of resistance [23,28], that makes sense. What usually defeats checkpoint drugs is less the absence of immunity than the presence of a stromal wall. Desert tumors are quiet but unguarded: once blockade wakes T cells, nothing holds them back. Stromal tumors fight back actively. ECO subtracts stroma for exactly this reason.

### ECO complements tumor mutational burden and marks inflamed tumors

TMB is the DNA rival, so we staged a direct contest in IMvigor210 where both measurements exist. ECO and TMB barely tracked each other (Spearman 0.26, p = 0.0001; Fig. 6A). They report on different bottlenecks. Where both were available (n = 214), TMB alone predicted response better (AUC 0.74 against ECO's 0.58), and the pair together (0.72) could not top TMB solo (Fig. 6B). We state that plainly. In bladder cancer ECO does not upgrade TMB when TMB is already on the table. Its worth lies elsewhere: a transcriptomic fallback when DNA-based TMB is missing, slow or silent, and a signal that holds within TMB strata. Gene-level detail backs the inflamed picture. ECO-high bladder tumors cranked up checkpoints and effector signals (CD274, PDCD1, CTLA4, LAG3, IFNG, CXCL9, CXCL13) while stromal drivers (FAP, TGFB1) sat low (Fig. 6C).

### Pathway and regulatory architecture of the two programs

Which biology do the lists capture? Pathway analysis split them neatly (Fig. 7A–B). TLS genes belonged to chemokine signaling (adjusted p = 7.9×10⁻²⁴), cytokine–receptor crosstalk (3.1×10⁻²³) and Toll-like receptor alarms: the vocabulary of immune-cell homing. Stromal genes belonged to matrix organization (4.1×10⁻²⁸), focal adhesion and collagen processing: the vocabulary of scar construction. Regulators split just as sharply. TLS genes take orders from NF-kB and RELA plus interferon sensors IRF1 and IRF3, the classic lymphoid-organizing switches [15,16]. Stromal genes take orders from EMT drivers TWIST2, ETV4, SP1, KLF8 and SRF, the fibroblast-activation crew [25,26]. So ECO pits two independently governed machines against each other: an NF-kB and IRF immune-organizing machine and an EMT and TGF-beta scar-building machine. Both machines predate any one cancer type, which helps explain why their balance reads out across melanoma, lung and bladder disease [25,28].

### Network hubs and miR-29 control of the stromal program

Genes work in crews, so we charted the wiring (Fig. 7C–D). STRING supplied 517 high-confidence protein links among our genes, and the 1,363 pooled TCGA tumors supplied co-expression across the same set. Both maps took the same shape. A tight TLS and chemokine cluster. A tight stromal and collagen cluster. A bridge of checkpoint and effector genes spanning them. Connection counts crowned CD8A, IFNG and IL10, the homing axis of CXCL12, CXCR4, CCR7 and CXCL9, TGFB1 on the stromal flank, and exhaustion nodes PDCD1, TIGIT and HAVCR2 [37]. Hubs earn attention because so much biology routes through them, which makes them the obvious intervention shortlist.

Regulation delivered one gem. Mining validated microRNA–target records, the miR-29 family (miR-29b at p = 2.9×10⁻¹⁸, with miR-29a and miR-29c close behind) directly suppresses 15 of the 30 stromal genes: TGFB1, TGFB2, TGFB3, multiple collagens, SPARC, LOX enzymes, MMP2. Nothing comparable governed the TLS list. This echoes established fibrosis biology, where miR-29 acts as the natural brake on scarring and its loss licenses desmoplasia [38]. A stromal-high tumor may therefore be, at heart, a miR-29–deficient tumor. Restoring that brake with miR-29 mimics becomes a testable sensitizing move (Fig. 8A).

### Drug mapping, comparison with TIDE and a proposed treatment framework

Hubs invite drugs. Querying 12 of them against the Drug–Gene Interaction Database returned 333 compound links (Fig. 8B, Table 4). Approved checkpoint agents appear (atezolizumab, durvalumab, ipilimumab, tremelimumab). So do approved stromal-axis drugs: the antifibrotic pirfenidone around TGF-beta signaling, bevacizumab and aflibercept around VEGF, plerixafor, motixafortide and mavorixafor around CXCR4. Experimental slots cover FAP, the CD39 and CD73 adenosine axis, IDO1 and MMP14. Table 4 converts these into trial sketches for ECO-low patients. CXCL12-heavy exclusion suggests checkpoint plus CXCR4 blockade [39]. Collagen and TGF-beta heaviness suggests checkpoint plus TGF-beta co-blockade [23]. History counsels humility here: IDO1 inhibition once looked equally destined and then collapsed in phase III [40]. Every pairing needs biomarker-stratified trials rather than blind faith.

Against TIDE [41], the reigning transcriptomic predictor, ECO won on TIDE's own turf. Using independently precomputed TIDE calls for IMvigor210, ECO reached AUC 0.60 against TIDE's 0.49 (n = 298; combined 0.62; Fig. 9A). Joint ECO–TMB grouping showed TMB-high patients responding at 36 to 44 percent whatever their ECO, and TMB-low near 10 percent whatever their ECO (Fig. 9B). That cements ECO's positioning in bladder cancer: a TMB substitute, not a TMB booster. Fig. 9C assembles the pieces into a treatment sketch. Sequence bulk RNA, compute ECO. ECO-high and TLS-high go to checkpoint blockade. ECO-low and stromal go to trials pairing checkpoints with the matching stromal-axis drug. Re-biopsy at progression and re-score, because ecosystems shift under therapy. A hypothesis, labeled as one, built for prospective testing.
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
