"""Literature-grounded gene signatures for TLS-vs-Stroma Ecosystem Score."""

TLS = ["CCL2","CCL3","CCL4","CCL5","CCL8","CCL18","CCL19","CCL21","CXCL9","CXCL10","CXCL11","CXCL13",
       "MS4A1","CD79A","CD79B","LTB","CCR7","CCR6","CXCR5","SELL","CD1D","LAT","SKAP1","CETP"]
# TLS: 12-chemokine signature (Messina/Provost; Cabrita et al. Nature 2020) + B-cell/Tfh/TLS-imprint
# genes (Sautes-Fridman et al.; Dieu-Nosjean et al.)

STROMA = ["FAP","ACTA2","COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL6A1","COL6A2","COL6A3",
          "POSTN","VIM","PDPN","PDGFRA","PDGFRB","TGFB1","TGFB2","TGFB3","TGFBR2","MMP2","MMP11",
          "MMP14","LOX","LOXL2","SPARC","FN1","VCAN","THY1","ITGB1","ZEB1"]
# STROMA/CAF: Mariathasan et al. Nature 2018 (TGFb/CAF exclusion), Hugo et al. Cell 2016 (IPRES/mesenchymal)

IFNG6 = ["IDO1","CXCL10","CXCL9","HLA-DRA","STAT1","IFNG"]  # Ayers et al. JCI 2017
CD8EFF = ["CD8A","CD8B","GZMA","GZMB","PRF1","IFNG"]
CHECKPOINTS = ["CD274","PDCD1","CTLA4","LAG3","TIGIT","HAVCR2","PDCD1LG2"]
AXES = ["CXCL12","CXCR4","VEGFA","VEGFC","IL10","ENTPD1","NT5E","HIF1A","ARG1","MKI67"]
LINEAGE = ["PTPRC","EPCAM","KRT19","GAPDH","ACTB"]

def union():
    seen, out = set(), []
    for g in TLS + STROMA + IFNG6 + CD8EFF + CHECKPOINTS + AXES + LINEAGE:
        if g not in seen:
            seen.add(g); out.append(g)
    return out
