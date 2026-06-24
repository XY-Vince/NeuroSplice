import pandas as pd
import os

PROJECT_DIR = "/Volumes/Untitled/NeuroSplice"
RAW_LIST = os.path.join(PROJECT_DIR, "results/string/significant_genes.txt")
FINAL_CSV = os.path.join(PROJECT_DIR, "results/final_gene_list.csv")
RMATS_SE = os.path.join(PROJECT_DIR, "results/rmats_hippocampus/SE.MATS.JCEC.txt") # Checking Hippo for Dctn1

print("--- 1. Reconciling Counts ---")
# Load Raw List
with open(RAW_LIST, 'r') as f:
    raw_genes = set(line.strip() for line in f if line.strip())
print(f"Raw List Count: {len(raw_genes)}")

# Load Final CSV
final_df = pd.read_csv(FINAL_CSV)
final_genes = set(final_df['GeneSymbol'])
print(f"Final CSV Count: {len(final_genes)}")

# Find Missing
missing = raw_genes - final_genes
print(f"Missing Genes ({len(missing)}): {list(missing)[:10]}...")

# Explanation: The cleaning script enforced |dPSI| > 0.1 AND FDR < 0.05.
# The raw list might have come from a less stringent filter or union of raw files.

print("\n--- 2. Dctn1 Coordinate Extraction ---")
# Extract Dctn1 coordinates from Hippocampus SE
try:
    se_df = pd.read_csv(RMATS_SE, sep='\t')
    dctn1 = se_df[se_df['geneSymbol'] == 'Dctn1']
    if not dctn1.empty:
        # Get the most significant event
        best_event = dctn1.sort_values('FDR').iloc[0]
        print(f"Gene: Dctn1")
        print(f"Event: SE")
        print(f"ID: {best_event['ID']}")
        print(f"Chr: {best_event['chr']}")
        print(f"Exon Start (0-based): {best_event['exonStart_0base']}")
        print(f"Exon End: {best_event['exonEnd']}")
        print(f"Upstream Exon End: {best_event['upstreamES']}")
        print(f"Downstream Exon Start: {best_event['downstreamES']}")
        print(f"FDR: {best_event['FDR']}")
        print(f"dPSI: {best_event['IncLevelDifference']}")
    else:
        print("Dctn1 not found in Hippocampus SE file.")
except Exception as e:
    print(f"Error reading rMATS file: {e}")

print("\n--- 3. DE Overlap Check ---")
# Check DESeq2 results
de_dir = os.path.join(PROJECT_DIR, "results/expression/deseq2_results")
if os.path.exists(de_dir):
    print(f"Checking {de_dir}...")
    files = [f for f in os.listdir(de_dir) if f.endswith('.txt') or f.endswith('.csv')]
    if not files:
        print("No DE result files found (consistent with 0 DE genes).")
    else:
        for f in files:
            print(f"Found: {f}")
            # Try to read and count sig genes
            # Assuming typical deseq2 output columns
else:
    print("DESeq2 results directory not found.")
