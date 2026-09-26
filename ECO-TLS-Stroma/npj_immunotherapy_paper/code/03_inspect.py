"""Inspect every downloaded file: shapes, headers, clinical characteristics."""
import os, gzip, tarfile, glob
import pandas as pd

RAW = "/home/user/npj_immunotherapy_paper/data/raw"
os.chdir(RAW)

# extract firehose data files
for c in ["SKCM", "LUAD", "BLCA"]:
    with tarfile.open(f"FIREHOSE_{c}_RNA.tar.gz") as t:
        for m in t.getmembers():
            if m.name.endswith(".data.txt"):
                m.name = os.path.basename(m.name)
                t.extract(m, path="firehose_rna")
    with tarfile.open(f"FIREHOSE_{c}_CLIN.tar.gz") as t:
        for m in t.getmembers():
            if m.name.endswith("clin.merged.txt"):
                m.name = f"{c}.clin.merged.txt"
                t.extract(m, path="firehose_clin")

print("=== Firehose RNA ===")
for f in sorted(glob.glob("firehose_rna/*.txt")):
    df = pd.read_csv(f, sep="\t", nrows=3)
    print(f, df.shape, list(df.columns[:4]), "| gene col sample:", df.iloc[0, 0], "|", df.iloc[1, 0])

print("\n=== Firehose clinical (transposed) ===")
for f in sorted(glob.glob("firehose_clin/*.txt")):
    df = pd.read_csv(f, sep="\t", header=None, nrows=6)
    print(f, "rows/cols:", df.shape)
    print("  attrs:", [str(x)[:45] for x in df.iloc[:, 0].tolist()])

def series_meta(gse):
    print(f"\n=== {gse} series_matrix characteristics ===")
    n = 0
    with gzip.open(f"{gse}_series_matrix.txt.gz", "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_characteristics_ch1"):
                print("  " + line.strip()[:600]); n += 1
                if n >= 6: break
    # sample list
    with gzip.open(f"{gse}_series_matrix.txt.gz", "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Series_sample_id"):
                ids = line.strip().split('"')[1].split()
                print(f"  n_samples={len(ids)} first3={ids[:3]}")
                break

for gse in ["GSE78220", "GSE91061", "GSE126044", "GSE135222", "GSE176307", "GSE207422"]:
    series_meta(gse)

print("\n=== Expression supplements ===")
print("--- GSE91061 fpkm ---")
df = pd.read_csv("GSE91061_fpkm.csv.gz", nrows=2)
print(df.shape, list(df.columns[:5]))
print("--- GSE126044 counts ---")
df = pd.read_csv("GSE126044_counts.txt.gz", sep="\t", nrows=2)
print(df.shape, list(df.columns[:6]))
print("--- GSE135222 exp ---")
df = pd.read_csv("GSE135222_exp.tsv.gz", sep="\t", nrows=2)
print(df.shape, list(df.columns[:6]))
print("--- GSE176307 logRNA ---")
df = pd.read_csv("GSE176307_logRNA.csv.gz", nrows=2)
print(df.shape, list(df.columns[:6]))
print("--- GSE176307 key ---")
df = pd.read_csv("GSE176307_key.csv.gz", nrows=3)
print(df.shape, list(df.columns))
print(df.to_string()[:800])
print("--- GSE207422 log2TPM ---")
df = pd.read_csv("GSE207422_log2TPM.txt.gz", sep="\t", nrows=2)
print(df.shape, list(df.columns[:6]))

print("\n--- GSE78220 xlsx sheets ---")
xl = pd.ExcelFile("GSE78220_PatientFPKM.xlsx")
print(xl.sheet_names)
for s in xl.sheet_names[:4]:
    d = xl.parse(s, nrows=3)
    print(f" sheet={s} shape~{d.shape} cols={list(d.columns[:6])}")

print("\n--- GSE207422 metadata xlsx ---")
xl = pd.ExcelFile("GSE207422_metadata.xlsx")
print(xl.sheet_names)
for s in xl.sheet_names:
    d = xl.parse(s, nrows=4)
    print(f" sheet={s} shape~{d.shape} cols={list(d.columns)}")
    print(d.to_string()[:1200])
