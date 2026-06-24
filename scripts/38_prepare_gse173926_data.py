import pandas as pd
import requests
import io
import os

# 1. Fetch ENA Report
print("Fetching ENA Report for SRP319365...")
url = "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP319365&result=read_run&fields=run_accession,sample_alias,sample_title,fastq_ftp&format=tsv"
s = requests.get(url).content
df_ena = pd.read_csv(io.StringIO(s.decode('utf-8')), sep='\t')

print(f"Loaded {len(df_ena)} runs from ENA.")

# 2. Load our Metadata to identify relevant GSMs
input_csv = "data/metadata/GSE173926_summary.csv"
if not os.path.exists(input_csv):
    print("Metadata CSV missing!")
    exit(1)
    
df_meta = pd.read_csv(input_csv)
# Filter for PFC Adult
target_df = df_meta[ (df_meta['Tissue'] == 'PFC') & (df_meta['Age_Stage'] == 'Adult') ]
target_gsms = target_df['Sample_ID'].tolist()
print(f"Targeting {len(target_gsms)} PFC samples.")

# 3. Join ENA and Metadata
# ENA 'sample_alias' corresponds to GSM ID
merged = pd.merge(target_df, df_ena, left_on='Sample_ID', right_on='sample_alias', how='inner')

if len(merged) == 0:
    print("Merge failed! Check GSM IDs in ENA report.")
    print("ENA Sample Aliases:", df_ena['sample_alias'].head().tolist())
    exit(1)

# 4. Generate Outputs
# Map File: GSM, SRR, Group, FTP
map_file = "data/metadata/GSE173926_srr_map.txt"
url_file = "data/metadata/GSE173926_urls.txt"

with open(map_file, 'w') as f_map, open(url_file, 'w') as f_url:
    f_map.write("Sample_ID\tRun_ID\tGroup\tGenotype\n")
    
    for idx, row in merged.iterrows():
        gsm = row['Sample_ID']
        srr = row['run_accession']
        genotype = row['Genotype']
        group = "WT" if genotype == "WT" else "Het" 
        ftps = row['fastq_ftp'].split(';')
        
        f_map.write(f"{gsm}\t{srr}\t{group}\t{genotype}\n")
        
        for ftp in ftps:
            f_url.write(f"http://{ftp}\n")

print(f"Generated {map_file} and {url_file}")
print(f"Total Runs: {len(merged)}")
