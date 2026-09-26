"""Parse all cohorts into standardized processed files.
expr: genes (rows, HGNC symbol) x samples (cols), log2-scale values.
clin: sample_id, response01 (1=responder), response_label, time, event, extra.
"""
import os, gzip, re
import numpy as np
import pandas as pd

RAW = "/home/user/npj_immunotherapy_paper/data/raw"
PROC = "/home/user/npj_immunotherapy_paper/data/processed"
os.makedirs(PROC, exist_ok=True)

def parse_series_matrix(gse):
    """Return dict: sample_ids (GSM order), titles, characteristics list of dicts."""
    path = f"{RAW}/{gse}_series_matrix.txt.gz"
    sample_ids, titles, chars, extra = [], [], [], {}
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Series_sample_id"):
                sample_ids = line.strip().split('"')[1].split()
            elif line.startswith("!Sample_title"):
                titles = re.findall(r'"([^"]*)"', line)
            elif line.startswith("!Sample_characteristics_ch1"):
                vals = re.findall(r'"([^"]*)"', line)
                chars.append(vals)
            elif line.startswith("!Sample_geo_accession"):
                pass
    return sample_ids, titles, chars

def chars_to_df(sample_ids, chars):
    """Each chars row: list like 'key: value' per sample -> DataFrame samples x keys."""
    recs = []
    for j in range(len(sample_ids)):
        d = {}
        for row in chars:
            if j < len(row) and ":" in row[j]:
                k, v = row[j].split(":", 1)
                d[k.strip()] = v.strip()
        recs.append(d)
    return pd.DataFrame(recs, index=sample_ids)

# ============ 1. TCGA Firehose ============
def parse_tcga(code):  # SKCM / LUAD / BLCA
    rna_path = f"{RAW}/firehose_rna/{code}.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt"
    print(f"[{code}] reading RNA...", flush=True)
    df = pd.read_csv(rna_path, sep="\t", low_memory=False)
    genes = df["Hybridization REF"].astype(str)
    sym = genes.str.split("|").str[0]
    df = df.drop(columns=["Hybridization REF"])
    df.index = sym
    df = df[~df.index.isin(["?", "gene_id"])]
    df = df.groupby(df.index).max()  # collapse duplicates
    # keep tumor samples (01 primary, 06 metastatic)
    keep = [c for c in df.columns if c[13:15] in ("01", "06")]
    df = df[keep]
    expr = np.log2(df.astype(np.float32) + 1)
    expr.columns = [c[:12] for c in expr.columns]  # patient barcode
    expr = expr.loc[:, ~expr.columns.duplicated(keep="first")]
    # clinical
    clin = pd.read_csv(f"{RAW}/firehose_clin/{code}.clin.merged.txt", sep="\t", low_memory=False)
    clin = clin.set_index(clin.columns[0])
    attrs = clin.index.str.lower()
    def find(pat):
        m = [i for i, a in enumerate(attrs) if pat in a]
        return clin.index[m[0]] if m else None
    k_death = find("days_to_death")
    k_fu = find("days_to_last_followup") or find("days_to_last_follow_up")
    k_vs = [a for a in clin.index if "vital_status" in a.lower()]
    print(f"[{code}] clin keys: death={k_death}, fu={k_fu}, vital={k_vs[:2] if k_vs else None}")
    k_vs = k_vs[0] if k_vs else None
    pats = [c for c in clin.columns]
    recs = []
    for p in pats:
        try: d = float(clin.loc[k_death, p]) if k_death else np.nan
        except: d = np.nan
        try: f = float(clin.loc[k_fu, p]) if k_fu else np.nan
        except: f = np.nan
        vs = str(clin.loc[k_vs, p]).lower() if k_vs else ""
        event = 1 if "dead" in vs else 0
        t = d if event == 1 and not np.isnan(d) else f
        recs.append((p, t, event))
    cdf = pd.DataFrame(recs, columns=["sample_id", "time", "event"])
    common = [s for s in expr.columns if s in set(cdf.sample_id)]
    expr = expr[common]
    cdf = cdf.set_index("sample_id").loc[common].reset_index()
    expr.to_csv(f"{PROC}/TCGA-{code}_expr.csv")
    cdf.to_csv(f"{PROC}/TCGA-{code}_clin.csv", index=False)
    print(f"[{code}] expr={expr.shape}, clin n={len(cdf)}, events={cdf.event.sum()}")

for code in ["SKCM", "LUAD", "BLCA"]:
    parse_tcga(code)

# ============ 2. GSE78220 melanoma anti-PD-1 ============
print("[GSE78220] parsing...", flush=True)
sids, titles, chars = parse_series_matrix("GSE78220")
cdf = chars_to_df(sids, chars)
cdf["gsm"] = sids
print("  char keys:", list(cdf.columns))
xl = pd.read_excel(f"{RAW}/GSE78220_PatientFPKM.xlsx", sheet_name="FPKM")
xl = xl.set_index("Gene")
xl.columns = [c.replace(".baseline", "") for c in xl.columns]  # Pt1...
mp = dict(zip(cdf["patient id"], cdf["anti-pd-1 response"]))
resp = {"Complete Response": 1, "Partial Response": 1, "Progressive Disease": 0}
cols, labels = [], []
for pt in xl.columns:
    if pt in mp and mp[pt] in resp:
        cols.append(pt); labels.append((pt, resp[mp[pt]], mp[pt]))
expr = np.log2(xl[cols].astype(np.float32) + 1)
pd.DataFrame(labels, columns=["sample_id", "response01", "response_label"]).to_csv(f"{PROC}/GSE78220_clin.csv", index=False)
expr.to_csv(f"{PROC}/GSE78220_expr.csv")
print(f"  expr={expr.shape}, responders={sum(l[1] for l in labels)}/{len(labels)}")

# ============ 3. GSE91061 melanoma nivolumab (Pre only) ============
print("[GSE91061] parsing...", flush=True)
sids, titles, chars = parse_series_matrix("GSE91061")
cdf = chars_to_df(sids, chars)
cdf["gsm"] = sids
print("  char keys:", list(cdf.columns), "| titles sample:", titles[:3])
fp = pd.read_csv(f"{RAW}/GSE91061_fpkm.csv.gz")
gcol = fp.columns[0]
fp = fp.set_index(gcol)
fp.index = fp.index.astype(str).str.split(".").str[0]  # strip version? inspect gene format
print("  gene id sample:", list(fp.index[:3]), "| col sample:", list(fp.columns[:2]))
# need gene symbols: check if index already symbols
# align by order: fp columns order should match series matrix order
print("  n_fpkm_cols:", len(fp.columns), "n_gsm:", len(sids))
order_visit = []
for c in fp.columns:
    m = re.search(r"_(Pre|On)_", c)
    order_visit.append(m.group(1) if m else "?")
sm_visit = [(v or "?") for v in cdf.get("visit (pre or on treatment)", ["?"] * len(sids))]
match = sum(1 for a, b in zip(order_visit, sm_visit) if a.lower() == b.lower())
print(f"  Pre/On order match: {match}/{len(sids)}")
cdf["fpkm_col"] = list(fp.columns)
pre = cdf[(cdf["visit (pre or on treatment)"] == "Pre") & (cdf["response"].isin(["PRCR", "PD", "SD"]))]
print(f"  Pre with known response: {len(pre)}")
rmap = {"PRCR": 1, "PD": 0, "SD": 0}
# gene symbol mapping: index may be hg19KnownGene IDs -> need conversion; try direct symbol match first
expr = np.log2(fp[[c for c in pre["fpkm_col"]]].astype(np.float32) + 1)
expr.columns = [f"Pt{i}" for i in range(expr.shape[1])]
out = pd.DataFrame({"sample_id": expr.columns, "response01": [rmap[r] for r in pre["response"]], "response_label": list(pre["response"])})
out.to_csv(f"{PROC}/GSE91061_clin.csv", index=False)
expr.to_csv(f"{PROC}/GSE91061_expr.csv")
print(f"  expr={expr.shape}, responders={(out.response01==1).sum()}/{len(out)}")

# ============ 4. GSE126044 NSCLC ============
print("[GSE126044] parsing...", flush=True)
sids, titles, chars = parse_series_matrix("GSE126044")
cdf = chars_to_df(sids, chars)
cdf["gsm"] = sids
print("  titles:", titles)
ct = pd.read_csv(f"{RAW}/GSE126044_counts.txt.gz", sep="\t")
print("  gene col sample:", list(ct.iloc[:2, 0]), "| cols:", list(ct.columns))
gcol = ct.columns[0]
ct = ct.set_index(gcol)
cpm = ct.div(ct.sum(axis=0), axis=1) * 1e6
expr = np.log2(cpm.astype(np.float32) + 1)
# map Dis_XX columns to GSM via titles
col2gsm = dict(zip(titles, sids))
mapped = [(c, col2gsm.get(c, None)) for c in expr.columns]
print("  col->GSM mapped:", sum(1 for _, g in mapped if g is not None), "/", len(mapped))
r = dict(zip(sids, cdf["patient response"]))
labels = []
keep = []
for c, g in mapped:
    if g in r and r[g] in ("responder", "non-responder"):
        keep.append(c); labels.append((c, 1 if r[g] == "responder" else 0, r[g]))
expr = expr[keep]
pd.DataFrame(labels, columns=["sample_id", "response01", "response_label"]).to_csv(f"{PROC}/GSE126044_clin.csv", index=False)
expr.to_csv(f"{PROC}/GSE126044_expr.csv")
print(f"  expr={expr.shape}, responders={sum(l[1] for l in labels)}/{len(labels)}")

# ============ 5. GSE135222 NSCLC (PFS + response?) ============
print("[GSE135222] parsing...", flush=True)
sids, titles, chars = parse_series_matrix("GSE135222")
cdf = chars_to_df(sids, chars)
cdf["gsm"] = sids
print("  char keys:", list(cdf.columns))
print("  titles sample:", titles[:3])
print(cdf.iloc[:, :8].to_string()[:1500])
