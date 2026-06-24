import pandas as pd
import os
import glob

PROJECT_DIR = "/Volumes/Untitled/NeuroSplice"
OUT_FILE = os.path.join(PROJECT_DIR, "results", "significant_genes_metadata.csv")

# Regions and their file paths
regions = {
    "Hippocampus": "rmats_hippocampus",
    "Prefrontal_Cortex": "rmats_prefrontal_cortex", 
    "Striatum": "rmats_striatum"
}

all_data = []

for region_name, folder in regions.items():
    # Check all rMATS event types
    for event in ["SE", "MXE", "A3SS", "A5SS", "RI"]:
        path = os.path.join(PROJECT_DIR, "results", folder, f"{event}.MATS.JCEC.txt")
        if os.path.exists(path):
            df = pd.read_csv(path, sep='\t')
            # Filter Significant
            sig = df[(df['FDR'] < 0.05) & (abs(df['IncLevelDifference']) > 0.1)].copy()
            
            if not sig.empty:
                # Select columns
                cols = sig[['geneSymbol', 'IncLevelDifference', 'FDR', 'PValue']].copy()
                cols['Region'] = region_name
                cols['Event_Type'] = event
                all_data.append(cols)

# Combine
if all_data:
    final_df = pd.concat(all_data)
    # Sort by FDR
    final_df = final_df.sort_values('FDR')
    
    # Save
    final_df.to_csv(OUT_FILE, index=False)
    print(f"Saved {len(final_df)} events (covering {final_df['geneSymbol'].nunique()} unique genes) to {OUT_FILE}")
else:
    print("No significant genes found to export.")
