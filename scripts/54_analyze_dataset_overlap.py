#!/usr/bin/env python3
"""
Deep-dive into GSE173926 vs GSE117357 overlap
==============================================
Understand why Category A genes don't replicate and what this means.
"""

import pandas as pd
import os

GSE173926_RMATS = "/Volumes/Untitled/NeuroSplice/results/gse173926_rmats_rerun"
GSE117357_ANALYSIS = "/Volumes/Untitled/NeuroSplice/results/rmats_gse117357_analysis"
OUTPUT_DIR = "/Volumes/Untitled/NeuroSplice/results/cross_dataset_validation"
EVENT_TYPES = ["SE", "A5SS", "A3SS", "MXE", "RI"]

FDR_THRESHOLD = 0.05
DPSI_THRESHOLD = 0.1

def load_all_genes(rmats_dir, dataset_name):
    """Load all significant genes from a dataset."""
    all_events = []
    
    for event_type in EVENT_TYPES:
        filepath = os.path.join(rmats_dir, f"{event_type}.MATS.JCEC.txt")
        
        if not os.path.exists(filepath):
            continue
        
        df = pd.read_csv(filepath, sep="\t")
        
        if 'geneSymbol' in df.columns:
            df['geneSymbol'] = df['geneSymbol'].str.replace('"', '')
        
        sig_df = df[
            (df['FDR'] < FDR_THRESHOLD) & 
            (df['IncLevelDifference'].abs() > DPSI_THRESHOLD)
        ].copy()
        
        if not sig_df.empty:
            sig_df['event_type'] = event_type
            sig_df['dataset'] = dataset_name
            all_events.append(sig_df)
    
    if all_events:
        return pd.concat(all_events, ignore_index=True)
    return pd.DataFrame()

def analyze_overlap():
    """Analyze gene overlap between datasets."""
    
    print("="*70)
    print("Deep-Dive: GSE173926 (MYT1L) vs GSE117357 (Rbfox1)")
    print("="*70)
    
    # Load both datasets
    print("\nLoading GSE173926 (MYT1L+/-)...")
    gse173926 = load_all_genes(GSE173926_RMATS, "MYT1L")
    genes_173926 = set(gse173926['geneSymbol'].unique()) if not gse173926.empty else set()
    print(f"  Significant genes: {len(genes_173926)}")
    
    print("\nLoading GSE117357 (Rbfox1 KO) - all tissues combined...")
    tissues = ["hippocampus", "prefrontal", "striatum"]
    all_117357 = []
    for tissue in tissues:
        tissue_dir = os.path.join("/Volumes/Untitled/NeuroSplice/results/rmats_gse117357", tissue)
        tissue_events = load_all_genes(tissue_dir, f"Rbfox1_{tissue}")
        if not tissue_events.empty:
            all_117357.append(tissue_events)
    
    gse117357 = pd.concat(all_117357, ignore_index=True) if all_117357 else pd.DataFrame()
    genes_117357 = set(gse117357['geneSymbol'].unique()) if not gse117357.empty else set()
    print(f"  Significant genes (any tissue): {len(genes_117357)}")
    
    # Calculate overlap
    overlap_genes = genes_173926 & genes_117357
    print(f"\n  Overlap: {len(overlap_genes)} genes ({len(overlap_genes)/len(genes_117357)*100:.1f}% of Rbfox1 genes)")
    
    # Interpretation
    print("\n" + "="*70)
    print("INTERPRETATION")
    print("="*70)
    
    overlap_pct = len(overlap_genes)/len(genes_117357)*100 if genes_117357 else 0
    
    if overlap_pct < 10:
        print("✓ LOW OVERLAP (<10%): Rbfox1 and MYT1L regulate DISTINCT gene sets")
        print("  → Category A genes are Rbfox1-SPECIFIC")
        print("  → Strengthens the biological relevance of your findings")
    elif overlap_pct < 30:
        print("⚠ MODEST OVERLAP (10-30%): Some shared targets")
        print("  → Partial regulatory overlap")
    else:
        print("⚠ HIGH OVERLAP (>30%): Significant shared regulation")
        print("  → Would expect Category A genes to replicate")
    
    # Top overlapping genes
    if overlap_genes:
        print(f"\nTop 10 overlapping genes (both datasets):")
        overlap_df = gse173926[gse173926['geneSymbol'].isin(overlap_genes)]
        top_genes = overlap_df.groupby('geneSymbol')['IncLevelDifference'].apply(
            lambda x: abs(x).mean()
        ).sort_values(ascending=False).head(10)
        
        for gene, dpsi in top_genes.items():
            print(f"  {gene:15s}: Mean |ΔΨ|={dpsi:.3f}")
    
    # Check Category A specifically
    cat_a_genes = ['Pts', 'Myo9b', 'Lrp8']
    print(f"\nCategory A genes in datasets:")
    for gene in cat_a_genes:
        in_173926 = gene in genes_173926
        in_117357 = gene in genes_117357
        print(f"  {gene:10s}: Rbfox1={in_117357}, MYT1L={in_173926}")
    
    # Top MYT1L genes
    print(f"\nTop 10 MYT1L-specific genes (for comparison):")
    mytl_specific = genes_173926 - genes_117357
    if mytl_specific and not gse173926.empty:
        mytl_df = gse173926[gse173926['geneSymbol'].isin(mytl_specific)]
        top_mytl = mytl_df.groupby('geneSymbol')['IncLevelDifference'].apply(
            lambda x: abs(x).mean()
        ).sort_values(ascending=False).head(10)
        
        for gene, dpsi in top_mytl.items():
            print(f"  {gene:15s}: Mean |ΔΨ|={dpsi:.3f}")
    
    # Save results
    summary = {
        'Dataset': ['GSE173926_MYT1L', 'GSE117357_Rbfox1', 'Overlap'],
        'Genes': [len(genes_173926), len(genes_117357), len(overlap_genes)],
        'Overlap_Pct': [
            len(overlap_genes)/len(genes_173926)*100 if genes_173926 else 0,
            len(overlap_genes)/len(genes_117357)*100 if genes_117357 else 0,
            100.0
        ]
    }
    
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(os.path.join(OUTPUT_DIR, "dataset_overlap_summary.csv"), index=False)
    
    if overlap_genes:
        overlap_list = pd.DataFrame({'Gene': sorted(list(overlap_genes))})
        overlap_list.to_csv(os.path.join(OUTPUT_DIR, "overlapping_genes.csv"), index=False)
    
    print(f"\nResults saved to: {OUTPUT_DIR}")

if __name__ == "__main__":
    analyze_overlap()
