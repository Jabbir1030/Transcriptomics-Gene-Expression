"""Memory-safe parser: grep-extract signature genes from big files, parse all cohorts.
Outputs data/processed/{cohort}_expr.csv (symbol x sample, log2) + {cohort}_clin.csv."""
import os, sys, gzip, re, subprocess
import numpy as np
import pandas as pd
sys.path.insert(0, "/home/user/npj_immunotherapy_paper/code")
from signatures import union

RAW = "/home/user/npj_immunotherapy_paper/data/raw"
WORK = "/home/user/npj_immunotherapy_paper/data/work"
PROC = "/home/user/npj_immunotherapy_paper/data/processed"
for d in [WORK, PROC]:
    os.makedirs(d, exist_ok=True)

GENES = union()
print(f"Gene universe: {len(GENES)} genes", flush=True)

# ---------- gene maps from NCBI gene_info ----------
print("Building symbol->Entrez/Ensembl maps...", flush=True)
sym2ent, sym2ens = {}, {}
with gzip.open("/home/user/npj_immunotherapy_paper/Homo_sapiens.gene_info.gz", "rt", errors="replace") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 10 or p[0] != "9606":
            continue
        ent, sym, xrefs = p[1], p[2], p[5]
        if sym in GENES and sym not in sym2ent:
            sym2ent[sym] = ent
        if sym in GENES and sym not in sym2ens:
            m = re.search(r"Ensembl:(ENSG\d+)", xrefs)
            if m:
                sym2ens[sym] = m.group(1)
print(f"  mapped Entrez: {len(sym2ent)}/{len(GENES)}; Ensembl: {len(sym2ens)}/{len(GENES)}", flush=True)
print("  missing Entrez:", [g for g in GENES if g not in sym2ent], flush=True)
print("  missing Ensembl:", [g for g in GENES if g not in sym2ens], flush=True)
ent2sym = {v: k for k, v in sym2ent.items()}
ens2sym = {v: k for k, v in sym2ens.items()}

def run(cmd):
    subprocess.run(cmd, shell=True, check=True, cwd=RAW)

# ---------- 1. TCGA Firehose (grep extract) ----------
for code in ["SKCM", "LUAD", "BLCA"]:
    src = f"firehose_rna/{code}.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt"
    pat = f"{WORK}/{code}.pat"
    with open(f"{RAW}/{pat}" if False else pat, "w") as f:
        for g in GENES:
            f.write(f"{g}|\n")
    out = f"{WORK}/{code}_sub.tsv"
    run(f"(head -1 {src} > {out}) && grep -F -f {pat} {src} >> {out} && wc -l {out}")
    df = pd.read_csv(out, sep="\t", low_memory=False)
    sym = df["Hybridization REF"].astype(str).str.split("|").str[0]
    df = df.drop(columns=["Hybridization REF"])
    df.index = sym
    df = df[df.index.isin(GENES)].groupby(level=0).max()
    keep = [c for c in df.columns if len(c) >= 15 and c[13:15] in ("01", "06")]
    df = df[keep]
    expr = np.log2(df.astype(np.float32) + 1)
    expr.columns = [c[:12] for c in expr.columns]
    expr = expr.loc[:, ~expr.columns.duplicated(keep="first")]
    # clinical (grep only needed rows; file has no header, barcodes in a row)
    with open(f"{WORK}/{code}_clin.pat", "w") as f:
        for a in ["patient.bcr_patient_barcode", "patient.days_to_death", "patient.days_to_last_followup",
                  "patient.vital_status", "patient.age_at_initial_pathologic_diagnosis", "patient.gender",
                  "patient.stage_event.pathologic_stage"]:
            f.write(a + "\n")
    run(f"grep -F -f {WORK}/{code}_clin.pat firehose_clin/{code}.clin.merged.txt > {WORK}/{code}_clin_sub.tsv")
    clin = pd.read_csv(f"{WORK}/{code}_clin_sub.tsv", sep="\t", header=None, low_memory=False)
    clin = clin.set_index(0)
    clin.columns = [str(x).upper() for x in clin.loc["patient.bcr_patient_barcode"]]
    clin = clin.drop(index="patient.bcr_patient_barcode")
    print(f"[{code}] clin rows={list(clin.index)} n_pats={clin.shape[1]}", flush=True)
    def row(key):
        return clin.loc[key] if key in clin.index else pd.Series(np.nan, index=clin.columns)
    dth, fup, vs = row("patient.days_to_death"), row("patient.days_to_last_followup"), row("patient.vital_status")
    age, sex = row("patient.age_at_initial_pathologic_diagnosis"), row("patient.gender")
    stage = row("patient.stage_event.pathologic_stage")
    recs = []
    for p in clin.columns:
        try: dd = float(dth[p])
        except: dd = np.nan
        try: ff = float(fup[p])
        except: ff = np.nan
        ev = 1 if str(vs[p]).lower() == "dead" else 0
        t = dd if (ev == 1 and not np.isnan(dd)) else ff
        recs.append((p, t, ev, age[p], sex[p], stage[p]))
    cdf = pd.DataFrame(recs, columns=["sample_id", "time", "event", "age", "sex", "stage"])
    cdf = cdf.dropna(subset=["time"])
    common = [s for s in expr.columns if s in set(cdf.sample_id)]
    expr = expr[common]
    cdf = cdf.set_index("sample_id").loc[common].reset_index()
    expr.to_csv(f"{PROC}/TCGA-{code}_expr.csv")
    cdf.to_csv(f"{PROC}/TCGA-{code}_clin.csv", index=False)
    print(f"[{code}] expr={expr.shape} genes_found={len(df)} clin_n={len(cdf)} events={int(cdf.event.sum())}", flush=True)

# ---------- helpers for series matrix ----------
def parse_series_matrix(gse):
    path = f"{RAW}/{gse}_series_matrix.txt.gz"
    sids, titles, chars = [], [], []
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Series_sample_id"):
                sids = line.strip().split('"')[1].split()
            elif line.startswith("!Sample_title"):
                titles = re.findall(r'"([^"]*)"', line)
            elif line.startswith("!Sample_characteristics_ch1"):
                chars.append(re.findall(r'"([^"]*)"', line))
    return sids, titles, chars

def chars_to_df(sids, chars):
    recs = []
    for j in range(len(sids)):
        d = {}
        for r in chars:
            if j < len(r) and ":" in r[j]:
                k, v = r[j].split(":", 1)
                d[k.strip()] = v.strip()
        recs.append(d)
    return pd.DataFrame(recs, index=sids)

# ---------- 2. GSE78220 ----------
print("[GSE78220] parsing...", flush=True)
sids, titles, chars = parse_series_matrix("GSE78220")
cdf = chars_to_df(sids, chars)
xl = pd.read_excel(f"{RAW}/GSE78220_PatientFPKM.xlsx", sheet_name="FPKM").set_index("Gene")
xl.columns = [c.replace(".baseline", "") for c in xl.columns]
xl = xl.loc[[g for g in GENES if g in xl.index]]
mp = dict(zip(cdf["patient id"], cdf["anti-pd-1 response"]))
rmap = {"Complete Response": 1, "Partial Response": 1, "Progressive Disease": 0}
cols = [c for c in xl.columns if c in mp and mp[c] in rmap]
expr = np.log2(xl[cols].astype(np.float32) + 1)
pd.DataFrame({"sample_id": cols, "response01": [rmap[mp[c]] for c in cols],
              "response_label": [mp[c] for c in cols]}).to_csv(f"{PROC}/GSE78220_clin.csv", index=False)
expr.to_csv(f"{PROC}/GSE78220_expr.csv")
print(f"  expr={expr.shape} n={len(cols)} resp={sum(rmap[mp[c]] for c in cols)}", flush=True)

# ---------- 3. GSE91061 (grep by Entrez) ----------
print("[GSE91061] extracting...", flush=True)
ents = [sym2ent[g] for g in GENES if g in sym2ent]
run(f"zcat GSE91061_fpkm.csv.gz | head -1 > {WORK}/GSE91061_sub.csv && zcat GSE91061_fpkm.csv.gz | grep -E '^(\"?)({"|".join(ents)}),\"' >> {WORK}/GSE91061_sub.csv; zcat GSE91061_fpkm.csv.gz | grep -E '^({"|".join(ents)}),' >> {WORK}/GSE91061_sub.csv; wc -l {WORK}/GSE91061_sub.csv")
fp = pd.read_csv(f"{WORK}/GSE91061_sub.csv")
fp.iloc[:, 0] = fp.iloc[:, 0].astype(str).str.replace('"', "", regex=False)
fp = fp[fp.iloc[:, 0].isin(ent2sym)].copy()
fp.iloc[:, 0] = fp.iloc[:, 0].map(ent2sym)
fp = fp.groupby(fp.columns[0]).max()
sids, titles, chars = parse_series_matrix("GSE91061")
cdf = chars_to_df(sids, chars)
cdf["title"] = titles
cdf = cdf.set_index("title")
rmap = {"PRCR": 1, "PD": 0, "SD": 0}
pre = cdf[(cdf["visit (pre or on treatment)"] == "Pre") & (cdf["response"].isin(rmap))]
cols = [c for c in fp.columns if c in pre.index]
expr = np.log2(fp[cols].astype(np.float32) + 1)
expr.columns = [f"S{i}" for i in range(len(cols))]
pd.DataFrame({"sample_id": expr.columns, "response01": [rmap[pre.loc[c, "response"]] for c in cols],
              "response_label": [pre.loc[c, "response"] for c in cols]}).to_csv(f"{PROC}/GSE91061_clin.csv", index=False)
expr.to_csv(f"{PROC}/GSE91061_expr.csv")
print(f"  expr={expr.shape} n={len(cols)} resp={sum(rmap[pre.loc[c,'response']] for c in cols)}", flush=True)

# ---------- 4. GSE126044 ----------
print("[GSE126044] parsing...", flush=True)
sids, titles, chars = parse_series_matrix("GSE126044")
cdf = chars_to_df(sids, chars)
ct = pd.read_csv(f"{RAW}/GSE126044_counts.txt.gz", sep="\t", index_col=0)
ct = ct.loc[[g for g in GENES if g in ct.index]]
cpm = ct.div(ct.sum(axis=0), axis=1) * 1e6
expr = np.log2(cpm.astype(np.float32) + 1)
t2g = {t.replace("RNA-seq_", ""): s for t, s in zip(titles, sids)}
r = dict(zip(sids, cdf["patient response"]))
keep = [c for c in expr.columns if c in t2g and r.get(t2g[c]) in ("responder", "non-responder")]
expr = expr[keep]
pd.DataFrame({"sample_id": keep, "response01": [1 if r[t2g[c]] == "responder" else 0 for c in keep],
              "response_label": [r[t2g[c]] for c in keep]}).to_csv(f"{PROC}/GSE126044_clin.csv", index=False)
expr.to_csv(f"{PROC}/GSE126044_expr.csv")
print(f"  expr={expr.shape} n={len(keep)}", flush=True)

# ---------- 5. GSE135222 (Ensembl; DCB from PFS) ----------
print("[GSE135222] parsing...", flush=True)
sids, titles, chars = parse_series_matrix("GSE135222")
cdf = chars_to_df(sids, chars)
cdf["title"] = [t.replace(" ", "") for t in titles]
df = pd.read_csv(f"{RAW}/GSE135222_exp.tsv.gz", sep="\t")
df["symbol"] = df["gene_id"].astype(str).str.split(".").str[0].map(ens2sym)
df = df.dropna(subset=["symbol"]).groupby("symbol").max(numeric_only=True)
df = df.loc[[g for g in GENES if g in df.index]]
expr = np.log2(df.astype(np.float32) + 1)
t2row = cdf.set_index("title")
keep, labs = [], []
for c in expr.columns:
    if c not in t2row.index:
        continue
    ev = int(t2row.loc[c, "progression-free survival (pfs)"])
    tm = float(t2row.loc[c, "pfs.time"])
    if ev == 1 and tm < 180:
        lab = (c, 0, "NDB", tm, ev)
    elif tm >= 180:
        lab = (c, 1, "DCB", tm, ev)
    else:
        continue  # censored <180d: unknown
    keep.append(c); labs.append(lab)
expr = expr[keep]
pd.DataFrame(labs, columns=["sample_id", "response01", "response_label", "time", "event"]).to_csv(f"{PROC}/GSE135222_clin.csv", index=False)
expr.to_csv(f"{PROC}/GSE135222_expr.csv")
print(f"  expr={expr.shape} n={len(keep)} DCB={sum(l[1] for l in labs)}", flush=True)

# ---------- 6. GSE207422 (MPR + RECIST) ----------
print("[GSE207422] parsing...", flush=True)
meta = pd.read_excel(f"{RAW}/GSE207422_metadata.xlsx", sheet_name="sheet1")
df = pd.read_csv(f"{RAW}/GSE207422_log2TPM.txt.gz", sep="\t", index_col=0)
df = df.loc[[g for g in GENES if g in df.index]]
pre = meta[meta["Resource"].str.contains("Pre", na=False)]
rmap = {"MPR": 1, "MPR (pCR)": 1, "NMPR": 0}
pre = pre[pre["Pathologic Response"].isin(rmap)]
cols = [c for c in df.columns if c in set(pre["Sample"])]
expr = df[cols].astype(np.float32)  # already log2TPM
m = pre.set_index("Sample")
pd.DataFrame({"sample_id": cols, "response01": [rmap[m.loc[c, "Pathologic Response"]] for c in cols],
              "response_label": [m.loc[c, "Pathologic Response"] for c in cols],
              "recist": [m.loc[c, "RECIST"] for c in cols]}).to_csv(f"{PROC}/GSE207422_clin.csv", index=False)
expr.to_csv(f"{PROC}/GSE207422_expr.csv")
print(f"  expr={expr.shape} n={len(cols)} MPR={sum(rmap[m.loc[c,'Pathologic Response']] for c in cols)}", flush=True)

print("\nALL PARSED. Processed files:", flush=True)
for f in sorted(os.listdir(PROC)):
    print("  ", f, flush=True)
