"""Step 1b: Download GEO supplementary expression tables + TCGA via cBioPortal."""
import os, requests

RAW = "/home/user/npj_immunotherapy_paper/data/raw"
os.makedirs(RAW, exist_ok=True)

def dl(url, out):
    if os.path.exists(out) and os.path.getsize(out) > 1000:
        print(f"SKIP: {os.path.basename(out)} [{os.path.getsize(out)/1e6:.2f} MB]")
        return True
    print(f"GET {os.path.basename(out)} ...", flush=True)
    try:
        r = requests.get(url, timeout=180, stream=True)
        if r.status_code != 200:
            print(f"  FAIL http={r.status_code}"); return False
        with open(out, "wb") as f:
            for ch in r.iter_content(chunk_size=1 << 20):
                if ch: f.write(ch)
        print(f"  OK [{os.path.getsize(out)/1e6:.2f} MB]")
        return True
    except Exception as e:
        print(f"  ERROR {e}"); return False

G = "https://ftp.ncbi.nlm.nih.gov/geo/series"
TASKS = [
    (f"{G}/GSE78nnn/GSE78220/suppl/GSE78220_PatientFPKM.xlsx", f"{RAW}/GSE78220_PatientFPKM.xlsx"),
    (f"{G}/GSE91nnn/GSE91061/suppl/GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz", f"{RAW}/GSE91061_fpkm.csv.gz"),
    (f"{G}/GSE126nnn/GSE126044/suppl/GSE126044_counts.txt.gz", f"{RAW}/GSE126044_counts.txt.gz"),
    (f"{G}/GSE135nnn/GSE135222/suppl/GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz", f"{RAW}/GSE135222_exp.tsv.gz"),
    (f"{G}/GSE176nnn/GSE176307/suppl/GSE176307_BACI_log_trans_normalized_RNAseq.csv.gz", f"{RAW}/GSE176307_logRNA.csv.gz"),
    (f"{G}/GSE176nnn/GSE176307/suppl/GSE176307_BACI_Omniseq_Sample_Name_Key_submitted_GEO_v2.csv.gz", f"{RAW}/GSE176307_key.csv.gz"),
    (f"{G}/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz", f"{RAW}/GSE207422_log2TPM.txt.gz"),
    (f"{G}/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx", f"{RAW}/GSE207422_metadata.xlsx"),
]
CB = "https://raw.githubusercontent.com/cBioPortal/datahub/master/public"
for study in ["skcm_tcga", "luad_tcga", "blca_tcga"]:
    TASKS += [
        (f"{CB}/{study}/data_RNA_Seq_v2_expression_median.txt", f"{RAW}/{study}_rna_median.txt"),
        (f"{CB}/{study}/data_clinical_patient.txt", f"{RAW}/{study}_clinical_patient.txt"),
        (f"{CB}/{study}/data_clinical_sample.txt", f"{RAW}/{study}_clinical_sample.txt"),
    ]

ok = sum(1 for u, o in TASKS if dl(u, o))
print(f"\n==== {ok}/{len(TASKS)} ok ====")
print("\n--- inventory ---")
for f in sorted(os.listdir(RAW)):
    p = os.path.join(RAW, f)
    print(f"{f:45s} {os.path.getsize(p)/1e6:8.2f} MB")
