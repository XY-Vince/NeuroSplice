import re
import urllib.request
import csv

url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173926"
print(f"Fetching {url}...")

with urllib.request.urlopen(url) as response:
    html = response.read().decode('utf-8')

# Regex to find GSM and Title
# Pattern: <a href="...acc=GSM5271282"...>GSM5271282</a>...<td valign="top">Adult_PFC_WT_Female_04 [RNA-Seq]</td>
# We'll look for the GSM ID and the "Sample name [RNA-Seq]" pattern which is robust here.

# Find all (GSM, Title) pairs
# We'll identify them by searching for the GSM link followed by the title cell
matches = re.findall(r'acc=(GSM\d+)".*?>\1</a></td>\s*<td valign="top">(.*?) \[RNA-Seq\]</td>', html, re.DOTALL)

print(f"Found {len(matches)} samples.")

csv_file = "data/metadata/GSE173926_summary.csv"
with open(csv_file, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Sample_ID', 'Sample_Name', 'Tissue', 'Genotype', 'Sex', 'Age_Stage'])
    
    for gsm, full_name in matches:
        # Parse name: Adult_PFC_WT_Female_04 or CTX_Hom_Female_01
        # Common parts: Tissue, Genotype, Sex
        
        parts = full_name.split('_')
        # Heuristic parsing
        tissue = "Unknown"
        genotype = "Unknown"
        sex = "Unknown"
        age = "Unknown"
        
        if "Adult_PFC" in full_name:
            age = "Adult"
            tissue = "PFC"
            # Format: Adult_PFC_WT_Female_04
            if len(parts) >= 4:
                genotype = parts[2]
                sex = parts[3]
        elif "CTX" in full_name:
            # CTX usually E14.5 based on study desc, but we'll mark as Embryonic/CTX
            age = "E14.5 (Inferred)"
            tissue = "Cortex"
            # Format: CTX_Hom_Female_01
            if len(parts) >= 3:
                genotype = parts[1]
                sex = parts[2]
                
        writer.writerow([gsm, full_name, tissue, genotype, sex, age])

print(f"Saved metadata to {csv_file}")
