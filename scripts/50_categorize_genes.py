#!/usr/bin/env python3
"""
Comprehensive Gene Categorization
==================================
Stratify ALL differential splicing genes into 3 categories:

Category 1: True Core (Same Exon)
   - Exact same genomic coordinates in all 3 tissues
   
Category 2: Multi-Event Core (Same Gene, Different Exons)  
   - Gene affected in all 3 tissues but different splice sites
   
Category 3: Tissue-Specific
   - Affected in 1-2 tissues only
   - High bar: |ΔΨ| > 0.15 for 2-tissue, > 0.20 for 1-tissue
"""

import pandas as pd
import os
from collections import defaultdict

# Configuration
RESULTS_DIR = "/Volumes/Untitled/NeuroSplice/results/rmats_gse117357"
OUTPUT_DIR = "/Volumes/Untitled/NeuroSplice/results/gene_categorization"
TISSUES = ["hippocampus", "prefrontal", "striatum"]
EVENT_TYPES = ["SE", "A5SS", "A3SS", "MXE", "RI"]

# Thresholds
FDR_THRESHOLD = 0.05
DPSI_BASE = 0.1        # Base threshold for all
DPSI_2TISSUE = 0.15    # Higher bar for 2-tissue genes
DPSI_1TISSUE = 0.20    # Even higher bar for 1-tissue genes

os.makedirs(OUTPUT_DIR, exist_ok=True)

def create_event_signature(row, event_type):
    """Create unique signature based on genomic coordinates."""
    try:
        if event_type == "SE":
            # Skipped exon
            sig = f"{row['chr']}:{row['strand']}:{row['exonStart_0base']}-{row['exonEnd']}"
        elif event_type == "A5SS":
            # Alt 5' splice site - use all relevant coordinates
            sig = f"{row['chr']}:{row['strand']}:{row['longExonStart_0base']}-{row['longExonEnd']}-{row['shortES']}-{row['shortEE']}"
        elif event_type == "A3SS":
            # Alt 3' splice site
            sig = f"{row['chr']}:{row['strand']}:{row['longExonStart_0base']}-{row['longExonEnd']}-{row['shortES']}-{row['shortEE']}"
        elif event_type == "MXE":
            # Mutually exclusive exons
            sig = f"{row['chr']}:{row['strand']}:{row['1stExonStart_0base']}-{row['1stExonEnd']}-{row['2ndExonStart_0base']}-{row['2ndExonEnd']}"
        elif event_type == "RI":
            # Retained intron
            sig = f"{row['chr']}:{row['strand']}:{row['riExonStart_0base']}-{row['riExonEnd']}"
        else:
            sig = None
        return sig
    except:
        return None

def load_all_significant_events():
    """Load all significant events from all tissues."""
    all_events = []
    
    for tissue in TISSUES:
        for event_type in EVENT_TYPES:
            filepath = os.path.join(RESULTS_DIR, tissue, f"{event_type}.MATS.JCEC.txt")
            
            if not os.path.exists(filepath):
                continue
            
            df = pd.read_csv(filepath, sep="\t")
            
            # Clean gene symbol
            if 'geneSymbol' in df.columns:
                df['geneSymbol'] = df['geneSymbol'].str.replace('"', '')
            
            # Filter significant
            sig_df = df[
                (df['FDR'] < FDR_THRESHOLD) & 
                (df['IncLevelDifference'].abs() > DPSI_BASE)
            ].copy()
            
            if not sig_df.empty:
                sig_df['tissue'] = tissue
                sig_df['event_type'] = event_type
                sig_df['event_signature'] = sig_df.apply(
                    lambda row: create_event_signature(row, event_type), axis=1
                )
                all_events.append(sig_df)
    
    if all_events:
        return pd.concat(all_events, ignore_index=True)
    return pd.DataFrame()

def categorize_genes(events_df):
    """Categorize genes into the 3 categories."""
    
    # Group by gene and event signature
    gene_event_groups = events_df.groupby(['geneSymbol', 'event_signature', 'event_type'])
    
    # Track genes and their event patterns
    gene_data = defaultdict(lambda: {
        'tissues': set(),
        'event_signatures': [],
        'events': []
    })
    
    for (gene, sig, etype), group in gene_event_groups:
        if sig is None:
            continue
            
        tissues = set(group['tissue'].unique())
        
        gene_data[gene]['event_signatures'].append({
            'signature': sig,
            'event_type': etype,
            'tissues': tissues,
            'num_tissues': len(tissues),
            'events': group.to_dict('records')
        })
        gene_data[gene]['tissues'].update(tissues)
    
    # Categorize
    category1 = []  # True Core (same exon in all 3 tissues)
    category2 = []  # Multi-Event Core (different exons in all 3 tissues)
    category3 = []  # Tissue-Specific
    
    for gene, data in gene_data.items():
        num_tissues = len(data['tissues'])
        
        # Check if any event signature appears in all 3 tissues
        has_true_core = False
        true_core_events = []
        
        for sig_data in data['event_signatures']:
            if sig_data['num_tissues'] == 3:
                has_true_core = True
                true_core_events.append(sig_data)
        
        if has_true_core:
            # Category 1: True Core
            for sig_data in true_core_events:
                # Get ΔΨ for each tissue
                tissue_dpsi = {}
                for event in sig_data['events']:
                    tissue_dpsi[event['tissue']] = event['IncLevelDifference']
                
                # Check direction consistency
                dpsi_values = list(tissue_dpsi.values())
                consistent = all(x > 0 for x in dpsi_values) or all(x < 0 for x in dpsi_values)
                
                category1.append({
                    'Gene': gene,
                    'EventType': sig_data['event_type'],
                    'EventSignature': sig_data['signature'],
                    'Hippo_ΔΨ': tissue_dpsi.get('hippocampus'),
                    'PFC_ΔΨ': tissue_dpsi.get('prefrontal'),
                    'Stri_ΔΨ': tissue_dpsi.get('striatum'),
                    'Consistent_Direction': consistent,
                    'Mean_|ΔΨ|': sum(abs(x) for x in dpsi_values) / len(dpsi_values)
                })
        
        elif num_tissues == 3:
            # Category 2: Multi-Event Core (gene in all 3 but different exons)
            # Get representative ΔΨ values (max absolute for each tissue)
            tissue_max_dpsi = defaultdict(float)
            event_types = set()
            
            for sig_data in data['event_signatures']:
                event_types.add(sig_data['event_type'])
                for event in sig_data['events']:
                    tissue = event['tissue']
                    dpsi = event['IncLevelDifference']
                    if abs(dpsi) > abs(tissue_max_dpsi[tissue]):
                        tissue_max_dpsi[tissue] = dpsi
            
            category2.append({
                'Gene': gene,
                'EventTypes': ','.join(sorted(event_types)),
                'NumEvents': len(data['event_signatures']),
                'Hippo_Max|ΔΨ|': abs(tissue_max_dpsi.get('hippocampus', 0)),
                'PFC_Max|ΔΨ|': abs(tissue_max_dpsi.get('prefrontal', 0)),
                'Stri_Max|ΔΨ|': abs(tissue_max_dpsi.get('striatum', 0)),
                'Mean_|ΔΨ|': sum(abs(tissue_max_dpsi[t]) for t in TISSUES) / 3
            })
        
        else:
            # Category 3: Tissue-Specific (1-2 tissues)
            # Apply higher thresholds
            tissues_list = sorted(list(data['tissues']))
            
            # Get max |ΔΨ| across all events for this gene
            max_dpsi = 0
            event_types = set()
            tissue_dpsi = defaultdict(float)
            
            for sig_data in data['event_signatures']:
                event_types.add(sig_data['event_type'])
                for event in sig_data['events']:
                    dpsi = abs(event['IncLevelDifference'])
                    if dpsi > max_dpsi:
                        max_dpsi = dpsi
                    tissue = event['tissue']
                    if dpsi > abs(tissue_dpsi[tissue]):
                        tissue_dpsi[tissue] = event['IncLevelDifference']
            
            # Apply threshold
            threshold = DPSI_2TISSUE if num_tissues == 2 else DPSI_1TISSUE
            
            if max_dpsi >= threshold:
                category3.append({
                    'Gene': gene,
                    'NumTissues': num_tissues,
                    'Tissues': ','.join(tissues_list),
                    'EventTypes': ','.join(sorted(event_types)),
                    'NumEvents': len(data['event_signatures']),
                    'Max_|ΔΨ|': max_dpsi,
                    'Hippo_ΔΨ': tissue_dpsi.get('hippocampus', 0),
                    'PFC_ΔΨ': tissue_dpsi.get('prefrontal', 0),
                    'Stri_ΔΨ': tissue_dpsi.get('striatum', 0),
                    'Passes_Threshold': f">={threshold}"
                })
    
    return (
        pd.DataFrame(category1),
        pd.DataFrame(category2),
        pd.DataFrame(category3)
    )

def main():
    print("="*70)
    print("Comprehensive Gene Categorization")
    print("="*70)
    
    # Load all events
    print("\nLoading all significant events...")
    events_df = load_all_significant_events()
    print(f"Total significant events: {len(events_df)}")
    print(f"Unique genes: {len(events_df['geneSymbol'].unique())}")
    
    # Categorize
    print("\nCategorizing genes...")
    cat1_df, cat2_df, cat3_df = categorize_genes(events_df)
    
    # Print results
    print("\n" + "="*70)
    print("CATEGORY 1: TRUE CORE (Same Exon in All 3 Tissues)")
    print("="*70)
    print(f"Total events: {len(cat1_df)}")
    print(f"Unique genes: {len(cat1_df['Gene'].unique()) if not cat1_df.empty else 0}")
    
    if not cat1_df.empty:
        print("\nGenes:")
        for gene in sorted(cat1_df['Gene'].unique()):
            gene_events = cat1_df[cat1_df['Gene'] == gene]
            print(f"\n  {gene}:")
            for _, row in gene_events.iterrows():
                consistency = "✓ Consistent" if row['Consistent_Direction'] else "✗ Mixed"
                print(f"    {row['EventType']:5s}: H={row['Hippo_ΔΨ']:+.3f}, "
                      f"P={row['PFC_ΔΨ']:+.3f}, S={row['Stri_ΔΨ']:+.3f} "
                      f"(Mean |ΔΨ|={row['Mean_|ΔΨ|']:.3f}) - {consistency}")
        
        cat1_df.to_csv(os.path.join(OUTPUT_DIR, "category1_true_core.csv"), index=False)
    
    print("\n" + "="*70)
    print("CATEGORY 2: MULTI-EVENT CORE (Different Exons in All 3 Tissues)")
    print("="*70)
    print(f"Total genes: {len(cat2_df)}")
    
    if not cat2_df.empty:
        # Sort by mean effect size
        cat2_df_sorted = cat2_df.sort_values('Mean_|ΔΨ|', ascending=False)
        
        print("\nTop 10 genes by mean |ΔΨ|:")
        for _, row in cat2_df_sorted.head(10).iterrows():
            print(f"  {row['Gene']:15s}: {row['NumEvents']} event(s), "
                  f"Mean |ΔΨ|={row['Mean_|ΔΨ|']:.3f}, "
                  f"Types: {row['EventTypes']}")
        
        cat2_df.to_csv(os.path.join(OUTPUT_DIR, "category2_multi_event_core.csv"), index=False)
    
    print("\n" + "="*70)
    print(f"CATEGORY 3: TISSUE-SPECIFIC (|ΔΨ| > {DPSI_2TISSUE} for 2-tissue, > {DPSI_1TISSUE} for 1-tissue)")
    print("="*70)
    print(f"Total genes: {len(cat3_df)}")
    
    if not cat3_df.empty:
        # Split by tissue count
        cat3_2tissue = cat3_df[cat3_df['NumTissues'] == 2]
        cat3_1tissue = cat3_df[cat3_df['NumTissues'] == 1]
        
        print(f"\n2-tissue genes: {len(cat3_2tissue)}")
        if not cat3_2tissue.empty:
            cat3_2tissue_sorted = cat3_2tissue.sort_values('Max_|ΔΨ|', ascending=False)
            print("  Top 5:")
            for _, row in cat3_2tissue_sorted.head(5).iterrows():
                print(f"    {row['Gene']:15s}: {row['Tissues']:30s}, Max |ΔΨ|={row['Max_|ΔΨ|']:.3f}")
        
        print(f"\n1-tissue genes: {len(cat3_1tissue)}")
        if not cat3_1tissue.empty:
            cat3_1tissue_sorted = cat3_1tissue.sort_values('Max_|ΔΨ|', ascending=False)
            print("  Top 5:")
            for _, row in cat3_1tissue_sorted.head(5).iterrows():
                print(f"    {row['Gene']:15s}: {row['Tissues']:30s}, Max |ΔΨ|={row['Max_|ΔΨ|']:.3f}")
        
        cat3_df.to_csv(os.path.join(OUTPUT_DIR, "category3_tissue_specific.csv"), index=False)
    
    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"Category 1 (True Core):           {len(cat1_df)} events, {len(cat1_df['Gene'].unique()) if not cat1_df.empty else 0} genes")
    print(f"Category 2 (Multi-Event Core):    {len(cat2_df)} genes")
    print(f"Category 3 (Tissue-Specific):     {len(cat3_df)} genes")
    total_genes = (
        (len(cat1_df['Gene'].unique()) if not cat1_df.empty else 0) +
        len(cat2_df) +
        len(cat3_df)
    )
    print(f"Total high-confidence genes:      {total_genes}")
    
    print(f"\nOutput directory: {OUTPUT_DIR}")

if __name__ == "__main__":
    main()
