#!/usr/bin/env python3
"""
Strict Core Event Analysis
===========================
Find genes where the EXACT SAME EXON shows differential splicing
across all 3 tissues (not just any event in the gene).
"""

import pandas as pd
import os
from collections import defaultdict

# Configuration
RESULTS_DIR = "/Volumes/Untitled/NeuroSplice/results/rmats_gse117357"
OUTPUT_DIR = "/Volumes/Untitled/NeuroSplice/results/strict_core_events"
TISSUES = ["hippocampus", "prefrontal", "striatum"]
EVENT_TYPES = ["SE", "A5SS", "A3SS", "MXE", "RI"]
FDR_THRESHOLD = 0.05
DPSI_THRESHOLD = 0.1

os.makedirs(OUTPUT_DIR, exist_ok=True)

def create_event_signature(row, event_type):
    """Create unique signature for an exon/junction based on coordinates."""
    
    if event_type == "SE":
        # Skipped exon: use exon coordinates
        sig = f"{row['chr']}:{row['exonStart_0base']}-{row['exonEnd']}"
    elif event_type == "A5SS":
        # Alt 5' splice site
        sig = f"{row['chr']}:{row['longExonStart_0base']}-{row['longExonEnd']}-{row['shortES']}-{row['shortEE']}"
    elif event_type == "A3SS":
        # Alt 3' splice site
        sig = f"{row['chr']}:{row['longExonStart_0base']}-{row['longExonEnd']}-{row['shortES']}-{row['shortEE']}"
    elif event_type == "MXE":
        # Mutually exclusive exons
        sig = f"{row['chr']}:{row['1stExonStart_0base']}-{row['1stExonEnd']}-{row['2ndExonStart_0base']}-{row['2ndExonEnd']}"
    elif event_type == "RI":
        # Retained intron
        sig = f"{row['chr']}:{row['riExonStart_0base']}-{row['riExonEnd']}"
    else:
        sig = None
    
    return sig

def load_events_with_signatures(tissue, event_type):
    """Load events and add coordinate-based signatures."""
    filepath = os.path.join(RESULTS_DIR, tissue, f"{event_type}.MATS.JCEC.txt")
    
    if not os.path.exists(filepath):
        return pd.DataFrame()
    
    df = pd.read_csv(filepath, sep="\t")
    
    # Clean gene symbol
    if 'geneSymbol' in df.columns:
        df['geneSymbol'] = df['geneSymbol'].str.replace('"', '')
    
    # Filter significant events
    sig_df = df[
        (df['FDR'] < FDR_THRESHOLD) & 
        (df['IncLevelDifference'].abs() > DPSI_THRESHOLD)
    ].copy()
    
    if not sig_df.empty:
        sig_df['tissue'] = tissue
        sig_df['event_type'] = event_type
        
        # Create event signatures
        sig_df['event_signature'] = sig_df.apply(
            lambda row: create_event_signature(row, event_type), axis=1
        )
    
    return sig_df

def find_strict_core_events():
    """Find events that occur in all 3 tissues with same exon coordinates."""
    
    print("="*60)
    print("Strict Core Event Analysis")
    print("="*60)
    print("Finding genes with IDENTICAL exon events across all tissues\n")
    
    # Load all significant events from all tissues
    all_events = []
    for tissue in TISSUES:
        for event_type in EVENT_TYPES:
            events = load_events_with_signatures(tissue, event_type)
            if not events.empty:
                all_events.append(events)
    
    if not all_events:
        print("No significant events found!")
        return pd.DataFrame()
    
    combined = pd.concat(all_events, ignore_index=True)
    print(f"Total significant events: {len(combined)}")
    
    # Group by gene + event_type + event_signature
    grouped = combined.groupby(['geneSymbol', 'event_type', 'event_signature'])
    
    strict_core = []
    for (gene, etype, sig), group in grouped:
        tissues_present = set(group['tissue'].unique())
        
        # Check if present in all 3 tissues
        if len(tissues_present) == 3:
            # Get ΔΨ values for each tissue
            tissue_dpsi = {}
            for tissue in TISSUES:
                tissue_rows = group[group['tissue'] == tissue]
                if not tissue_rows.empty:
                    tissue_dpsi[tissue] = tissue_rows['IncLevelDifference'].values[0]
            
            strict_core.append({
                'Gene': gene,
                'EventType': etype,
                'EventSignature': sig,
                'Hippo_ΔΨ': tissue_dpsi.get('hippocampus', None),
                'PFC_ΔΨ': tissue_dpsi.get('prefrontal', None),
                'Stri_ΔΨ': tissue_dpsi.get('striatum', None)
            })
    
    strict_core_df = pd.DataFrame(strict_core)
    
    print(f"\nStrict core events (same exon in all 3 tissues): {len(strict_core_df)}")
    
    if not strict_core_df.empty:
        # Count unique genes
        unique_genes = strict_core_df['Gene'].unique()
        print(f"Unique genes with strict core events: {len(unique_genes)}")
        
        print("\nGenes with strict core events:")
        for gene in sorted(unique_genes):
            gene_events = strict_core_df[strict_core_df['Gene'] == gene]
            print(f"  {gene:15s}: {len(gene_events)} event(s)")
            
            for _, row in gene_events.iterrows():
                h_dir = "↑WT" if row['Hippo_ΔΨ'] > 0 else "↓KO"
                p_dir = "↑WT" if row['PFC_ΔΨ'] > 0 else "↓KO"
                s_dir = "↑WT" if row['Stri_ΔΨ'] > 0 else "↓KO"
                
                # Check consistency
                all_same_sign = (
                    (row['Hippo_ΔΨ'] > 0 and row['PFC_ΔΨ'] > 0 and row['Stri_ΔΨ'] > 0) or
                    (row['Hippo_ΔΨ'] < 0 and row['PFC_ΔΨ'] < 0 and row['Stri_ΔΨ'] < 0)
                )
                consistency = "✓ Consistent" if all_same_sign else "✗ Mixed"
                
                print(f"    {row['EventType']:5s}: H={row['Hippo_ΔΨ']:+.3f} ({h_dir}), "
                      f"P={row['PFC_ΔΨ']:+.3f} ({p_dir}), "
                      f"S={row['Stri_ΔΨ']:+.3f} ({s_dir}) - {consistency}")
        
        # Save results
        strict_core_df.to_csv(os.path.join(OUTPUT_DIR, "strict_core_events.csv"), index=False)
        
        # Create summary by gene
        gene_summary = strict_core_df.groupby('Gene').agg({
            'EventType': lambda x: ','.join(x),
            'Hippo_ΔΨ': ['mean', 'min', 'max'],
            'PFC_ΔΨ': ['mean', 'min', 'max'],
            'Stri_ΔΨ': ['mean', 'min', 'max']
        }).reset_index()
        
        gene_summary.to_csv(os.path.join(OUTPUT_DIR, "strict_core_genes_summary.csv"), index=False)
        
        # Compare with original "core 6"
        original_core = ['Arid5a', 'Bcl2l11', 'Gm10419', 'Lrp8', 'Myo9b', 'Pts']
        print("\n" + "="*60)
        print("Comparison with Original 'Core 6'")
        print("="*60)
        
        for gene in original_core:
            if gene in unique_genes:
                print(f"  {gene:15s}: ✓ Has strict core event(s)")
            else:
                print(f"  {gene:15s}: ✗ No strict core events (different exons per tissue)")
        
        # New genes not in original 6
        new_strict_genes = set(unique_genes) - set(original_core)
        if new_strict_genes:
            print(f"\nAdditional genes with strict core events (not in original 6):")
            for gene in sorted(new_strict_genes):
                gene_events = strict_core_df[strict_core_df['Gene'] == gene]
                print(f"  {gene:15s}: {len(gene_events)} event(s)")
    else:
        print("\nNo strict core events found!")
        print("(No genes have the exact same exon affected in all 3 tissues)")
    
    return strict_core_df

def main():
    strict_df = find_strict_core_events()
    
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    
    if not strict_df.empty:
        print(f"Total strict core events: {len(strict_df)}")
        print(f"Unique genes: {len(strict_df['Gene'].unique())}")
        print(f"\nOutput saved to: {OUTPUT_DIR}")
    else:
        print("No genes meet the strict core event criterion")
        print("(Same exon differential splicing in all 3 tissues)")

if __name__ == "__main__":
    main()
