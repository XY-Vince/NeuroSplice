import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import sys

# Configuration
PROJECT_DIR = "/Volumes/Untitled/NeuroSplice"
RMATS_DIR = os.path.join(PROJECT_DIR, "results") # Parent dir for regional folders
OUT_DIR = os.path.join(PROJECT_DIR, "results", "visualizations")
REGIONS = ["Hippocampus", "Prefrontal_Cortex", "Striatum"]
CORE_GENES = ["Dctn1", "Crem", "Pdlim7", "Spata5"]

def ensure_dir(d):
    if not os.path.exists(d):
        os.makedirs(d)

def load_rmats_data(region):
    # Try different naming conventions
    path = os.path.join(RMATS_DIR, f"rmats_{region.lower()}", "SE.MATS.JCEC.txt")
    if not os.path.exists(path):
        path = os.path.join(RMATS_DIR, f"rmats_{region.lower()}_filtered", "SE.MATS.JCEC.txt")
    
    if not os.path.exists(path):
        print(f"Warning: Could not find SE file for {region}")
        return None
        
    df = pd.read_csv(path, sep='\t')
    return df

def plot_volcano(df, region):
    plt.figure(figsize=(10, 8))
    
    # Clean data
    df['dPSI'] = df['IncLevelDifference']
    # Cap extremely small FDRs to avoid compressing the plot
    # -log10(1e-50) = 50. Anything significantly higher suppresses other points.
    df['logFDR'] = -np.log10(df['FDR'] + 1e-300)
    
    # Clip Y-axis outliers
    Y_LIMIT = 20
    df['logFDR_plot'] = df['logFDR'].clip(upper=Y_LIMIT)
    
    # Color coding
    df['color'] = 'grey'
    # Significant
    sig_mask = (df['FDR'] < 0.05) & (abs(df['dPSI']) > 0.1)
    df.loc[sig_mask, 'color'] = 'red'
    
    # Core Genes
    core_mask = df['geneSymbol'].isin(CORE_GENES)
    df.loc[core_mask, 'color'] = 'blue'
    
    sns.scatterplot(data=df, x='dPSI', y='logFDR_plot', c=df['color'], alpha=0.6, s=30)
    
    # Mark clipped points (Outliers)
    outliers = df[df['logFDR'] > Y_LIMIT]
    if not outliers.empty:
        # Plot markers at the limit
        plt.scatter(outliers['dPSI'], [Y_LIMIT] * len(outliers), 
                   marker='D', c='#8B008B', s=60, label='Outliers (>20)', zorder=10) # DarkMagenta
        
        # Annotate each outlier with Gene and actual value
        for _, row in outliers.iterrows():
            plt.text(row['dPSI'], Y_LIMIT + 0.2, 
                     f"{row['geneSymbol']}\n(10^-{int(row['logFDR'])})", 
                     fontsize=9, ha='center', va='bottom', color='#8B008B', fontweight='bold')
    
    # Add a dashed line at the cap
    plt.axhline(Y_LIMIT, linestyle=':', color='#8B008B', alpha=0.5)
    
    # Label Core Genes
    core_hits = df[core_mask]
    for _, row in core_hits.iterrows():
        plt.text(row['dPSI'], row['logFDR_plot'], row['geneSymbol'], 
                 fontsize=11, fontweight='bold', color='darkblue')
                 
    plt.title(f"Volcano Plot: {region} (SE Events)\nRed: FDR<0.05, |dPSI|>0.1", fontsize=14)
    plt.xlabel("Delta PSI (KO - WT)", fontsize=12)
    plt.ylabel("-log10(FDR)", fontsize=12)
    plt.axhline(-np.log10(0.05), linestyle='--', color='black', alpha=0.5)
    plt.axvline(0.1, linestyle='--', color='black', alpha=0.5)
    plt.axvline(-0.1, linestyle='--', color='black', alpha=0.5)
    
    out_path = os.path.join(OUT_DIR, f"volcano_{region}.png")
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Saved {out_path}")

def plot_heatmap(df, region):
    # Filter for significant genes
    sig_df = df[(df['FDR'] < 0.05) & (abs(df['IncLevelDifference']) > 0.1)].copy()
    
    if len(sig_df) > 100:
        # Take top 100 by dPSI magnitude
        sig_df['abs_dPSI'] = abs(sig_df['IncLevelDifference'])
        sig_df = sig_df.sort_values('abs_dPSI', ascending=False).head(100)
    
    if len(sig_df) == 0:
        print(f"No significant genes for heatmap in {region}")
        return

    # Extract PSI values
    # IncLevel1 is WT (comma separated), IncLevel2 is KO
    # We need to parse these into a matrix
    # Rows: Genes, Cols: Samples
    
    matrix = []
    labels = []
    
    for _, row in sig_df.iterrows():
        wt_psis = [float(x) for x in row['IncLevel1'].split(',')]
        ko_psis = [float(x) for x in row['IncLevel2'].split(',')]
        
        # Combine
        # Note: We assume sample order is consistent. 
        # rMATS output preserves input order.
        full_row = wt_psis + ko_psis
        matrix.append(full_row)
        labels.append(row['geneSymbol'])
        
    matrix = np.array(matrix)
    
    # Create sample labels
    n_wt = len(wt_psis)
    n_ko = len(ko_psis)
    sample_labels = [f"WT_{i+1}" for i in range(n_wt)] + [f"KO_{i+1}" for i in range(n_ko)]
    
    plt.figure(figsize=(12, max(8, len(labels)*0.2)))
    df_plot = pd.DataFrame(matrix, index=labels, columns=sample_labels)
    
    sns.heatmap(df_plot, cmap="vlag", center=0.5, yticklabels=True)
    plt.title(f"Top Splicing Events: {region} (PSI Values)", fontsize=14)
    
    out_path = os.path.join(OUT_DIR, f"heatmap_{region}.png")
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved {out_path}")

if __name__ == "__main__":
    ensure_dir(OUT_DIR)
    for region in REGIONS:
        print(f"Processing {region}...")
        df = load_rmats_data(region)
        if df is not None:
            plot_volcano(df, region)
            plot_heatmap(df, region)
