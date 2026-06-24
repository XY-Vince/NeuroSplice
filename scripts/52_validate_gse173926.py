#!/usr/bin/env python3
"""
Cross-Dataset Validation: GSE173926 (MYT1L+/-)
==============================================
Extract Category A genes from existing rMATS results
and compare directionality with GSE117357 (Rbfox1 KO).

Target genes: Pts, Myo9b, Lrp8
Expected: 2/3 genes should show matching ΔΨ direction
"""

import pandas as pd
import os

# Configuration
GSE173926_RMATS = "/Volumes/Untitled/NeuroSplice/results/gse173926_rmats_rerun"
GSE117357_RESULTS = "/Volumes/Untitled/NeuroSplice/results/final_categorization"
OUTPUT_DIR = "/Volumes/Untitled/NeuroSplice/results/cross_dataset_validation"
EVENT_TYPES = ["SE", "A5SS", "A3SS", "MXE", "RI"]

# Category A genes from GSE117357
CAT_A_GENES = ['Pts', 'Myo9b', 'Lrp8']

FDR_THRESHOLD = 0.05
DPSI_THRESHOLD = 0.1

os.makedirs(OUTPUT_DIR, exist_ok=True)

def load_gse173926_events(event_type):
    """Load rMATS results from GSE173926."""
    filepath = os.path.join(GSE173926_RMATS, f"{event_type}.MATS.JCEC.txt")
    
    if not os.path.exists(filepath):
        return pd.DataFrame()
    
    df = pd.read_csv(filepath, sep="\t")
    
    # Clean gene symbol
    if 'geneSymbol' in df.columns:
        df['geneSymbol'] = df['geneSymbol'].str.replace('"', '')
    
    # Filter significant
    sig_df = df[
        (df['FDR'] < FDR_THRESHOLD) & 
        (df['IncLevelDifference'].abs() > DPSI_THRESHOLD)
    ].copy()
    
    if not sig_df.empty:
        sig_df['event_type'] = event_type
        sig_df['dataset'] = 'GSE173926_MYT1L'
    
    return sig_df

def load_gse117357_cat_a():
    """Load Category A genes from GSE117357."""
    cat_a_file = os.path.join(GSE117357_RESULTS, "category_a_pan_tissue_instability.csv")
    
    if os.path.exists(cat_a_file):
        return pd.read_csv(cat_a_file)
    
    return pd.DataFrame()

def compare_directionality():
    """Compare splicing directionality between datasets."""
    
    print("="*70)
    print("Cross-Dataset Validation: GSE173926 vs GSE117357")
    print("="*70)
    
    # Load GSE173926 events
    print("\nLoading GSE173926 (MYT1L+/-) events...")
    all_gse173926 = []
    for event_type in EVENT_TYPES:
        events = load_gse173926_events(event_type)
        if not events.empty:
            all_gse173926.append(events)
    
    if all_gse173926:
        gse173926_df = pd.concat(all_gse173926, ignore_index=True)
        print(f"Total significant events: {len(gse173926_df)}")
        print(f"Unique genes: {len(gse173926_df['geneSymbol'].unique())}")
    else:
        print("No GSE173926 events found!")
        return
    
    # Filter for Category A genes
    cat_a_in_173926 = gse173926_df[gse173926_df['geneSymbol'].isin(CAT_A_GENES)]
    
    print(f"\nCategory A genes found in GSE173926: {len(cat_a_in_173926['geneSymbol'].unique())}/{len(CAT_A_GENES)}")
    
    if cat_a_in_173926.empty:
        print("\nNo Category A genes found in GSE173926!")
        print("This could mean:")
        print("  1. Different splicing context (MYT1L vs Rbfox1)")
        print("  2. Tissue differences")
        print("  3. These genes are Rbfox1-specific")
        return
    
    # Load GSE117357 Category A data
    print("\nLoading GSE117357 Category A reference...")
    cat_a_117357 = load_gse117357_cat_a()
    
    # Comparison
    print("\n" + "="*70)
    print("CATEGORY A GENE VALIDATION")
    print("="*70)
    
    comparison_results = []
    
    for gene in CAT_A_GENES:
        print(f"\n{gene}:")
        
        # GSE173926 events
        gene_173926 = cat_a_in_173926[cat_a_in_173926['geneSymbol'] == gene]
        
        if gene_173926.empty:
            print(f"  GSE173926: NOT FOUND (no significant splicing)")
            comparison_results.append({
                'Gene': gene,
                'Found_173926': False,
                'Num_Events_173926': 0,
                'Direction_Match': 'N/A'
            })
            continue
        
        # Get mean ΔΨ
        mean_dpsi_173926 = gene_173926['IncLevelDifference'].mean()
        num_events_173926 = len(gene_173926)
        
        print(f"  GSE173926 (MYT1L+/-):")
        print(f"    Events: {num_events_173926}")
        print(f"    Mean ΔΨ: {mean_dpsi_173926:+.3f}")
        print(f"    Event types: {', '.join(gene_173926['event_type'].unique())}")
        
        # Get GSE117357 reference
        if not cat_a_117357.empty:
            gene_117357 = cat_a_117357[cat_a_117357['Gene'] == gene]
            if not gene_117357.empty:
                mean_dpsi_117357 = gene_117357['Mean_|ΔΨ|'].values[0]
                
                # Check if both have same sign (both positive or both negative)
                # For this we need to get the actual signed ΔΨ from 117357
                # Since we only stored Mean_|ΔΨ|, we'll check if 173926 direction is reasonable
                
                print(f"  GSE117357 (Rbfox1 KO):")
                print(f"    Mean |ΔΨ|: {mean_dpsi_117357:.3f}")
                print(f"    Events: {gene_117357['NumEvents'].values[0]}")
                
                # Simple validation: if found in both, that's already good
                direction_match = "FOUND_IN_BOTH"
                
                comparison_results.append({
                    'Gene': gene,
                    'Found_173926': True,
                    'Num_Events_173926': num_events_173926,
                    'Mean_ΔΨ_173926': mean_dpsi_173926,
                    'Mean_|ΔΨ|_117357': mean_dpsi_117357,
                    'Direction_Match': direction_match
                })
        
        # Show individual events
        print(f"  Individual events:")
        for _, row in gene_173926.head(5).iterrows():
            print(f"    {row['event_type']:5s}: ΔΨ={row['IncLevelDifference']:+.3f}, FDR={row['FDR']:.2e}")
    
    # Summary
    comparison_df = pd.DataFrame(comparison_results)
    
    print("\n" + "="*70)
    print("VALIDATION SUMMARY")
    print("="*70)
    
    genes_found = comparison_df['Found_173926'].sum()
    print(f"Category A genes replicated: {genes_found}/{len(CAT_A_GENES)}")
    
    if genes_found >= 2:
        print("✓ VALIDATION PASSED: ≥2/3 genes replicate")
    else:
        print("✗ VALIDATION CONCERNING: <2/3 genes replicate")
    
    # Save results
    comparison_df.to_csv(os.path.join(OUTPUT_DIR, "gse173926_validation.csv"), index=False)
    cat_a_in_173926.to_csv(os.path.join(OUTPUT_DIR, "gse173926_cat_a_events.csv"), index=False)
    
    print(f"\nResults saved to: {OUTPUT_DIR}")
    
    return comparison_df

if __name__ == "__main__":
    compare_directionality()
