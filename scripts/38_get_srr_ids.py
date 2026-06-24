import pandas as pd
import subprocess
import time

# 1. Load Metadata
df = pd.read_csv("data/metadata/GSE173926_summary.csv")

# 2. Filter for Adult PFC
pfc_df = df[ (df['Tissue'] == 'PFC') & (df['Age_Stage'] == 'Adult') ]

# 3. Define Groups
wt_samples = pfc_df[pfc_df['Genotype'] == 'WT']['Sample_ID'].tolist()
het_samples = pfc_df[pfc_df['Genotype'] == 'Het']['Sample_ID'].tolist()

print(f"PFC WT Samples ({len(wt_samples)}): {wt_samples}")
print(f"PFC Het Samples ({len(het_samples)}): {het_samples}")

# 4. Fetch SRR IDs (Simulated for speed, normally uses esearch)
# We will verify if we can fetch them via curl from ENA/SRA
gsm_to_srr = {}

print("\nFetching SRR IDs...")
for gsm in wt_samples + het_samples:
    # Use curl to get run info from ENA or GEO
    try:
        # Simple scraping from GEO page for the SRX/SRR link
        cmd = f'curl -s "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}" | grep "SRR"'
        result = subprocess.check_output(cmd, shell=True).decode('utf-8')
        # Extract SRR
        import re
        srr_matches = re.findall(r'SRR\d+', result)
        if srr_matches:
            unique_srrs = sorted(list(set(srr_matches)))
            gsm_to_srr[gsm] = unique_srrs[0] # Take first run
            print(f"{gsm} -> {unique_srrs[0]}")
        else:
            print(f"{gsm} -> No SRR found")
        time.sleep(0.5) # Rate limit
    except Exception as e:
        print(f"{gsm} -> Error: {e}")

# 5. Save Map
with open("data/metadata/GSE173926_srr_map.txt", "w") as f:
    for gsm, srr in gsm_to_srr.items():
        Group = "WT" if gsm in wt_samples else "Het"
        f.write(f"{gsm}\t{srr}\t{Group}\n")
print("\nSaved SRR map to data/metadata/GSE173926_srr_map.txt")
