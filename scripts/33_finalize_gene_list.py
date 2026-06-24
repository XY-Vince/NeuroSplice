import pandas as pd
import os
import numpy as np

PROJECT_DIR = "/Volumes/Untitled/NeuroSplice"
OUT_FILE = os.path.join(PROJECT_DIR, "results", "final_gene_list.csv")

# Regions and their file paths
regions = {
    "Hippocampus": "rmats_hippocampus",
    "Prefrontal_Cortex": "rmats_prefrontal_cortex", 
    "Striatum": "rmats_striatum"
}

all_events = []

print("Parsing rMATS files...")
for region_name, folder in regions.items():
    for event_type in ["SE", "MXE", "A3SS", "A5SS", "RI"]:
        path = os.path.join(PROJECT_DIR, "results", folder, f"{event_type}.MATS.JCEC.txt")
        if os.path.exists(path):
            df = pd.read_csv(path, sep='\t')
            # Filter Significant
            sig = df[(df['FDR'] < 0.05) & (abs(df['IncLevelDifference']) > 0.1)].copy()
            
            if not sig.empty:
                sig['Region'] = region_name
                sig['Event_Type'] = event_type
                # Ensure GeneID is string (Ensembl ID)
                sig['GeneID'] = sig['GeneID'].astype(str)
                all_events.append(sig[['GeneID', 'geneSymbol', 'IncLevelDifference', 'FDR', 'Region', 'Event_Type']])

if not all_events:
    print("No significant events found.")
    exit(1)

# Combine all events
full_df = pd.concat(all_events)

# 1. Analyze Directionality per Gene
print("Analyzing directionality...")
gene_stats = []

for gene, group in full_df.groupby('geneSymbol'):
    ensembl_id = group['GeneID'].iloc[0] # Take first ID
    n_events = len(group)
    mean_dpsi = group['IncLevelDifference'].mean()
    
    # Check if all events go in same direction
    all_pos = (group['IncLevelDifference'] > 0).all()
    all_neg = (group['IncLevelDifference'] < 0).all()
    
    classification = "Complex/Mixed"
    if all_pos:
        classification = "Unidirectional (Increased Retention)"
    elif all_neg:
        classification = "Unidirectional (Increased Skipping)"
        
    regions_present = ",".join(group['Region'].unique())
    
    gene_stats.append({
        'GeneSymbol': gene,
        'EnsemblID': ensembl_id,
        'Classification': classification,
        'N_Events': n_events,
        'Mean_dPSI': round(mean_dpsi, 3),
        'Regions': regions_present,
        'Min_FDR': group['FDR'].min()
    })

# Create Gene-Level Metadata DataFrame
gene_df = pd.DataFrame(gene_stats)
gene_df = gene_df.sort_values('Min_FDR')

# Save Master List
gene_df.to_csv(OUT_FILE, index=False)
print(f"Saved Master Gene List to {OUT_FILE}")
print(f"Total Unique Genes: {len(gene_df)}")

# Print Stats
print("\nClassification Summary:")
print(gene_df['Classification'].value_counts())

# Print Core 4 Check
print("\nCore 4 Status:")
core_genes = ["Crem", "Pdlim7", "Dctn1", "Spata5"]
print(gene_df[gene_df['GeneSymbol'].isin(core_genes)][['GeneSymbol', 'Classification', 'N_Events', 'Mean_dPSI']])
