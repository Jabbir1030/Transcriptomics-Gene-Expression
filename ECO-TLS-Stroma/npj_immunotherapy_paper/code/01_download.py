"""Step 1: Download all public cohorts (TCGA via UCSC Xena, ICI cohorts via GEO).
Real data only. Prints inventory with shapes. No fabrication."""
import os, gzip, shutil, requests, io

RAW = "/home/user/npj_immunotherapy_paper/data/raw"
os.makedirs(RAW, exist_ok=True)

def dl(url, out):
    if os.path.exists(out) and os.path.getsize(out) > 1000:
        print(f"SKIP (exists): {out} [{os.path.getsize(out)/1e6:.1f} MB]")
        return True
    print(f"GET {url}")
    try:
        r = requests.get(url, timeout=120, stream=True)
        if r.status_code != 200:
            print(f"  FAIL http={r.status_code}")
            return False
        with open(out, "wb") as f:
            for chunk in r.iter_content(chunk_size=1 << 20):
                if chunk: f.write(chunk)
        print(f"  OK -> {out} [{os.path.getsize(out)/1e6:.1f} MB]")
        return True
    except Exception as e:
        print(f"  ERROR {e}")
        return False

TASKS = []
# --- TCGA via UCSC Xena GDC hub ---
for cohort in ["TCGA-SKCM", "TCGA-LUAD", "TCGA-BLCA"]:
    TASKS.append((f"https://gdc-hub.s3.us-east-1.amazonaws.com/download/{cohort}.htseq_fpkm.tsv.gz",
                  f"{RAW}/{cohort}.htseq_fpkm.tsv.gz"))
    TASKS.append((f"https://gdc-hub.s3.us-east-1.amazonaws.com/download/{cohort}.survival.tsv",
                  f"{RAW}/{cohort}.survival.tsv"))
    TASKS.append((f"https://gdc-hub.s3.us-east-1.amazonaws.com/download/{cohort}.GDC_phenotype.tsv.gz",
                  f"{RAW}/{cohort}.GDC_phenotype.tsv.gz"))

# --- GEO series matrices (clinical + microarray expression where applicable) ---
def geo_matrix(gse):
    n = gse[3:]
    prefix = f"GSE{n[:2]}nnn" if len(n) >= 4 else "GSE112nnn"
    # generic: first digits
    if gse.startswith("GSE78"): prefix = "GSE78nnn"
    elif gse.startswith("GSE91"): prefix = "GSE91nnn"
    elif gse.startswith("GSE126"): prefix = "GSE126nnn"
    elif gse.startswith("GSE135"): prefix = "GSE135nnn"
    elif gse.startswith("GSE176"): prefix = "GSE176nnn"
    elif gse.startswith("GSE207"): prefix = "GSE207nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse}/matrix/{gse}_series_matrix.txt.gz"

for gse in ["GSE78220", "GSE91061", "GSE126044", "GSE135222", "GSE176307", "GSE207422"]:
    TASKS.append((geo_matrix(gse), f"{RAW}/{gse}_series_matrix.txt.gz"))

# --- GEO supplementary processed tables (RNA-seq cohorts) ---
TASKS += [
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78220/suppl/GSE78220_FPKM_Table.txt.gz",
     f"{RAW}/GSE78220_FPKM_Table.txt.gz"),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/suppl/GSE91061_COUNTs_oct2019.txt.gz",
     f"{RAW}/GSE91061_counts.txt.gz"),
]

ok, fail = 0, []
for url, out in TASKS:
    if dl(url, out): ok += 1
    else: fail.append(url)

print(f"\n==== DONE: {ok}/{len(TASKS)} ok ====")
if fail:
    print("FAILED:")
    for u in fail: print("  ", u)

# quick inventory
print("\n--- inventory ---")
for f in sorted(os.listdir(RAW)):
    p = os.path.join(RAW, f)
    print(f"{f:45s} {os.path.getsize(p)/1e6:8.1f} MB")
