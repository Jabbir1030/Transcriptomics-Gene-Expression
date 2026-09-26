"""Enhancement fetches: Enrichr pathways/TF/miRNA, STRING PPI, DGIdb drugs, TIDE comparison data."""
import sys, json, time
import requests
import pandas as pd
import numpy as np
sys.path.insert(0, "/home/user/npj_immunotherapy_paper/code")
from signatures import TLS, STROMA

TAB = "/home/user/npj_immunotherapy_paper/results/tables"
ENH = "/home/user/npj_immunotherapy_paper/results/enhancement"
import os
os.makedirs(ENH, exist_ok=True)

# ---------- 1. Enrichr: pathways + TF + miRNA ----------
def enrichr(genes, libs):
    r = requests.post("https://maayanlab.cloud/Enrichr/addList",
                      files={"list": (None, "\n".join(genes)), "description": (None, "ECO")}, timeout=60)
    r.raise_for_status()
    uid = r.json()["userListId"]
    out = {}
    for lib in libs:
        g = requests.get("https://maayanlab.cloud/Enrichr/enrich",
                         params={"userListId": uid, "backgroundType": lib}, timeout=60)
        g.raise_for_status()
        rows = []
        for e in (g.json().get(lib) or [])[:15]:
            rows.append({"term": e[1], "p": e[2], "adj_p": e[6], "odds": e[7] if len(e) > 7 else np.nan,
                         "genes": ";".join(e[5])})
        out[lib] = pd.DataFrame(rows)
        time.sleep(0.5)
    return out

LIBS = ["KEGG_2021_Human", "Reactome_2022", "TRRUST_Transcription_Factors_2019", "miRTarBase_2017"]
print("Enrichr TLS...", flush=True)
tls_en = enrichr(TLS, LIBS)
print("Enrichr STROMA...", flush=True)
str_en = enrichr(STROMA, LIBS)
for lib in LIBS:
    tls_en[lib].to_csv(f"{ENH}/enrich_TLS_{lib.split('_')[0]}.csv", index=False)
    str_en[lib].to_csv(f"{ENH}/enrich_STROMA_{lib.split('_')[0]}.csv", index=False)
    print(f"--- TLS {lib} top3 ---")
    print(tls_en[lib][["term", "adj_p"]].head(3).to_string(index=False))
    print(f"--- STROMA {lib} top3 ---")
    print(str_en[lib][["term", "adj_p"]].head(3).to_string(index=False))

# ---------- 2. STRING PPI ----------
print("\nSTRING PPI...", flush=True)
genes85 = sorted(set(TLS + STROMA + ["CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "IFNG", "IDO1", "CD274",
                                      "PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "PDCD1LG2", "CXCL12", "CXCR4",
                                      "VEGFA", "VEGFC", "IL10", "ENTPD1", "NT5E", "HIF1A", "ARG1", "MKI67",
                                      "PTPRC", "EPCAM", "KRT19"]))
r = requests.post("https://string-db.org/api/json/network",
                  data={"identifiers": "\n".join(genes85), "species": 9606,
                        "required_score": 700, "network_type": "functional"}, timeout=120)
r.raise_for_status()
edges = r.json()
ed = pd.DataFrame([{"A": e["preferredName_A"], "B": e["preferredName_B"], "score": float(e["score"])} for e in edges])
ed.to_csv(f"{ENH}/string_edges.csv", index=False)
print(f"STRING edges (score>=700): {len(ed)}", flush=True)
print(ed.head(8).to_string(index=False))

# ---------- 3. DGIdb drugs for hub genes ----------
print("\nDGIdb...", flush=True)
HUBS = ["FAP", "TGFB1", "POSTN", "VEGFA", "CXCL12", "CXCR4", "CD274", "CTLA4", "IDO1", "ENTPD1", "NT5E", "MMP14"]
q = """{ genes(names: [%s]) { nodes { name drugInteractions(first: 60) { nodes {
  interactionTypes { type } sources { sourceDbName } drug { name approved } } } } } }""" % ",".join(f'"{g}"' for g in HUBS)
r = requests.post("https://dgidb.org/api/v2/graphql", json={"query": q}, timeout=120)
print("DGIdb status:", r.status_code, flush=True)
recs = []
try:
    d = r.json()
    for node in d["data"]["genes"]["nodes"]:
        g = node["name"]
        for di in (node.get("drugInteractions") or {}).get("nodes", []):
            drug = (di.get("drug") or {}).get("name", "?")
            appr = (di.get("drug") or {}).get("approved", "")
            itypes = ",".join(t["type"] for t in (di.get("interactionTypes") or []))
            srcs = ",".join(s["sourceDbName"] for s in (di.get("sources") or []))
            recs.append({"gene": g, "drug": drug, "approved": appr, "type": itypes, "sources": srcs})
except Exception as e:
    print("DGIdb parse issue:", str(e)[:300], flush=True)
    print(r.text[:500], flush=True)
dg = pd.DataFrame(recs).drop_duplicates()
dg.to_csv(f"{ENH}/dgidb_drugs.csv", index=False)
print(f"DGIdb interactions: {len(dg)}", flush=True)
if len(dg):
    print(dg.groupby("gene").size().to_string(), flush=True)

# ---------- 4. TIDE + clinical benefit (IMvigor210, iAtlas precomputed) ----------
print("\ncBioPortal TIDE...", flush=True)
BASE = "https://www.cbioportal.org/api"
SID = "blca_iatlas_imvigor210_2017"
for a in ["TIDE_RESPONDER", "TIDE_NO_BENEFITS", "CLINICAL_BENEFIT", "PROGRESSION"]:
    r = requests.get(f"{BASE}/studies/{SID}/clinical-data",
                     params={"attributeId": a, "clinicalDataType": "SAMPLE", "projection": "SUMMARY", "pageSize": 1000}, timeout=120)
    d = {x["sampleId"]: x["value"] for x in r.json()}
    pd.DataFrame(list(d.items()), columns=["sample_id", a]).to_csv(f"{ENH}/imv_{a}.csv", index=False)
    vals = pd.Series(list(d.values())).value_counts().to_dict()
    print(f"{a}: n={len(d)} values={vals}", flush=True)
print("DONE enhance_fetch", flush=True)
