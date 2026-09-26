"""Fetch IMvigor210 (bladder, atezolizumab) RNA + clinical via cBioPortal API (GET endpoints)."""
import sys, time, gzip, requests
import pandas as pd
import numpy as np
sys.path.insert(0, "/home/user/npj_immunotherapy_paper/code")
from signatures import union

BASE = "https://www.cbioportal.org/api"
SID = "blca_iatlas_imvigor210_2017"
MP = "blca_iatlas_imvigor210_2017_rna_seq_mrna"
SL = SID + "_all"
PROC = "/home/user/npj_immunotherapy_paper/data/processed"

GENES = union()
gmap = {}
with gzip.open("/home/user/npj_immunotherapy_paper/Homo_sapiens.gene_info.gz", "rt", errors="replace") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) > 2 and p[0] == "9606" and p[2] in GENES and p[2] not in gmap:
            gmap[p[2]] = int(p[1])
print(f"genes: {len(gmap)}", flush=True)

# expression per gene
cols = {}
for i, (sym, ent) in enumerate(gmap.items()):
    r = requests.get(f"{BASE}/molecular-profiles/{MP}/molecular-data",
                     params={"sampleListId": SL, "entrezGeneId": ent, "projection": "SUMMARY"}, timeout=120)
    r.raise_for_status()
    d = r.json()
    cols[sym] = {x["sampleId"]: x["value"] for x in d}
    if (i + 1) % 20 == 0:
        print(f"  {i+1}/{len(gmap)} genes", flush=True)
expr = pd.DataFrame(cols).T
expr = expr.apply(pd.to_numeric, errors="coerce")
print(f"expr shape: {expr.shape}, range: {np.nanmin(expr.values):.2f}..{np.nanmax(expr.values):.2f}", flush=True)

# clinical per attribute
attrs = ["RESPONSE", "OS_MONTHS", "OS_STATUS", "TMB_NONSYNONYMOUS", "SEX", "CLINICAL_BENEFIT"]
clin = {}
for a in attrs:
    r = requests.get(f"{BASE}/studies/{SID}/clinical-data",
                     params={"attributeId": a, "clinicalDataType": "SAMPLE", "projection": "SUMMARY",
                             "pageSize": 1000}, timeout=120)
    r.raise_for_status()
    clin[a] = {x["sampleId"]: x["value"] for x in r.json()}
cdf = pd.DataFrame(clin)
print("RESPONSE:", cdf["RESPONSE"].value_counts().to_dict(), flush=True)
print("OS_STATUS:", cdf["OS_STATUS"].value_counts().to_dict(), flush=True)

resp_map = {"Complete Response": 1, "Partial Response": 1, "Stable Disease": 0, "Progressive Disease": 0,
            "CR": 1, "PR": 1, "SD": 0, "PD": 0}
common = [s for s in expr.columns if s in cdf.index]
expr = expr[common]
cdf = cdf.loc[common]
out = pd.DataFrame({"sample_id": common})
out["response_label"] = [str(cdf.loc[s].get("RESPONSE", "")) for s in common]
out["response01"] = [resp_map.get(x.strip(), np.nan) for x in out["response_label"]]
out["time"] = pd.to_numeric(cdf["OS_MONTHS"], errors="coerce") * 30.44
out["event"] = [1 if ("DECEASED" in str(x).upper() or str(x).strip() == "1") else 0 for x in cdf["OS_STATUS"]]
out["tmb"] = pd.to_numeric(cdf["TMB_NONSYNONYMOUS"], errors="coerce")
if np.nanmax(expr.values) > 50:
    expr = np.log2(expr.astype(np.float32) + 1)
    print("applied log2(x+1)", flush=True)
expr.to_csv(f"{PROC}/IMvigor210_expr.csv")
out.to_csv(f"{PROC}/IMvigor210_clin.csv", index=False)
print(f"saved: expr={expr.shape}, n_resp_known={int(out.response01.notna().sum())}, resp_rate={out.response01.mean():.3f}", flush=True)
print(out.head(6).to_string(), flush=True)
