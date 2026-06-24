#!/usr/bin/env python3
"""
Core Gene Quality Verification
===============================
Verify the 6 core genes identified across all 3 tissues:
1. Gene type (protein-coding vs lncRNA)
2. ΔΨ direction consistency across tissues
3. Junction novelty (annotated vs novel)
"""

import pandas as pd
import os
import gzip

# Configuration
RESULTS_DIR = "/Volumes/Untitled/NeuroSplice/results/rmats_gse117357"
GTF_FILE = "/Volumes/Untitled/NeuroSplice/reference/gencode.vM25.annotation.gtf"
OUTPUT_DIR = "/Volumes/Untitled/NeuroSplice/results/core_gene_verification"
TISSUES = ["hippocampus", "prefrontal", "striatum"]
EVENT_TYPES = ["SE", "A5SS", "A3SS", "MXE", "RI"]

CORE_GENES = ['Arid5a', 'Bcl2l11', 'Gm10419', 'Lrp8', 'Myo9b', 'Pts']

os.makedirs(OUTPUT_DIR, exist_ok=True)

def get_gene_biotype(gene_symbol, gtf_file):
    """Extract gene biotype from GTF file."""
    biotypes = {}
    
    print(f"  Searching GTF for gene types...")
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            if fields[2] == 'gene':
                attrs = fields[8]
                
                # Extract gene_name and gene_type
                gene_name = None
                gene_type = None
                
                for attr in attrs.split(';'):
                    attr = attr.strip()
                    if attr.startswith('gene_name'):
                        gene_name = attr.split('"')[1]
                    elif attr.startswith('gene_type'):
                        gene_type = attr.split('"')[1]
                
                if gene_name and gene_type:
                    biotypes[gene_name] = gene_type
    
    return biotypes

def load_gene_events(gene, tissue, event_type):
    """Load all events for a specific gene from rMATS output."""
    filepath = os.path.join(RESULTS_DIR, tissue, f"{event_type}.MATS.JCEC.txt")
    
    if not os.path.exists(filepath):
        return pd.DataFrame()
    
    df = pd.read_csv(filepath, sep="\t")
    
    # Clean gene symbol
    if 'geneSymbol' in df.columns:
        df['geneSymbol'] = df['geneSymbol'].str.replace('"', '')
    
    # Filter for this gene
    gene_events = df[df['geneSymbol'] == gene].copy()
    
    if not gene_events.empty:
        gene_events['tissue'] = tissue
        gene_events['event_type'] = event_type
    
    return gene_events

def check_junction_novelty(tissue, event_type):
    """Check which events use novel junctions."""
    
    # Check if novel junction file exists
    novel_file = os.path.join(RESULTS_DIR, tissue, f"fromGTF.novelJunction.{event_type}.txt")
    regular_file = os.path.join(RESULTS_DIR, tissue, f"fromGTF.{event_type}.txt")
    
    novel_ids = set()
    regular_ids = set()
    
    if os.path.exists(novel_file):
        df_novel = pd.read_csv(novel_file, sep="\t")
        novel_ids = set(df_novel['ID'].tolist()) if 'ID' in df_novel.columns else set()
    
    if os.path.exists(regular_file):
        df_regular = pd.read_csv(regular_file, sep="\t")
        regular_ids = set(df_regular['ID'].tolist()) if 'ID' in df_regular.columns else set()
    
    return novel_ids, regular_ids

def analyze_core_genes():
    """Comprehensive analysis of core genes."""
    
    print("="*60)
    print("Core Gene Quality Verification")
    print("="*60)
    
    # 1. Get gene biotypes
    print("\n1. Checking Gene Types...")
    biotypes = get_gene_biotype(None, GTF_FILE)
    
    gene_info = []
    for gene in CORE_GENES:
        biotype = biotypes.get(gene, 'UNKNOWN')
        is_protein_coding = biotype == 'protein_coding'
        
        gene_info.append({
            'Gene': gene,
            'Biotype': biotype,
            'Protein_Coding': is_protein_coding
        })
        
        status = "✓ Protein-coding" if is_protein_coding else f"✗ {biotype}"
        print(f"  {gene:12s}: {status}")
    
    gene_df = pd.DataFrame(gene_info)
    
    # 2. Analyze ΔΨ patterns across tissues
    print("\n2. Analyzing ΔΨ Patterns Across Tissues...")
    
    all_events = []
    for gene in CORE_GENES:
        gene_events = []
        
        for tissue in TISSUES:
            for event_type in EVENT_TYPES:
                events = load_gene_events(gene, tissue, event_type)
                if not events.empty:
                    gene_events.append(events)
        
        if gene_events:
            combined = pd.concat(gene_events, ignore_index=True)
            
            # Get significant events only
            sig_events = combined[
                (combined['FDR'] < 0.05) & 
                (combined['IncLevelDifference'].abs() > 0.1)
            ]
            
            all_events.append(sig_events)
            
            if not sig_events.empty:
                print(f"\n  {gene}:")
                for _, row in sig_events.iterrows():
                    dpsi = row['IncLevelDifference']
                    direction = "↑ More in WT" if dpsi > 0 else "↓ More in KO"
                    print(f"    {row['tissue']:15s} {row['event_type']:5s}: ΔΨ={dpsi:+.3f} ({direction})")
    
    if all_events:
        all_events_df = pd.concat(all_events, ignore_index=True)
    
    # 3. Check junction novelty
    print("\n3. Checking Junction Novelty...")
    
    junction_summary = []
    for gene in CORE_GENES:
        for tissue in TISSUES:
            for event_type in EVENT_TYPES:
                events = load_gene_events(gene, tissue, event_type)
                
                if not events.empty:
                    # Check if these event IDs appear in novel junction files
                    novel_ids, regular_ids = check_junction_novelty(tissue, event_type)
                    
                    for _, row in events.iterrows():
                        event_id = row['ID']
                        is_novel = event_id in novel_ids
                        is_annotated = event_id in regular_ids
                        
                        status = "Novel" if is_novel else ("Annotated" if is_annotated else "Unknown")
                        
                        if (row['FDR'] < 0.05) and (abs(row['IncLevelDifference']) > 0.1):
                            junction_summary.append({
                                'Gene': gene,
                                'Tissue': tissue,
                                'EventType': event_type,
                                'EventID': event_id,
                                'Status': status,
                                'ΔΨ': row['IncLevelDifference']
                            })
    
    junction_df = pd.DataFrame(junction_summary)
    
    if not junction_df.empty:
        print("\n  Junction Status Summary:")
        status_counts = junction_df.groupby(['Gene', 'Status']).size().reset_index(name='Count')
        for gene in CORE_GENES:
            gene_data = status_counts[status_counts['Gene'] == gene]
            if not gene_data.empty:
                print(f"    {gene}:")
                for _, row in gene_data.iterrows():
                    print(f"      {row['Status']:12s}: {row['Count']} events")
    
    # 4. Create summary report
    print("\n" + "="*60)
    print("SUMMARY REPORT")
    print("="*60)
    
    summary_file = os.path.join(OUTPUT_DIR, "core_gene_verification.txt")
    with open(summary_file, 'w') as f:
        f.write("Core Gene Quality Verification Report\n")
        f.write("="*60 + "\n\n")
        
        f.write("1. GENE TYPES\n")
        f.write("-"*60 + "\n")
        for _, row in gene_df.iterrows():
            status = "PASS" if row['Protein_Coding'] else "FAIL"
            f.write(f"{row['Gene']:12s}: {row['Biotype']:20s} [{status}]\n")
        
        protein_coding_count = gene_df['Protein_Coding'].sum()
        f.write(f"\nProtein-coding: {protein_coding_count}/6 genes\n")
        
        f.write("\n2. ΔΨ PATTERNS\n")
        f.write("-"*60 + "\n")
        
        if all_events:
            for gene in CORE_GENES:
                gene_events = all_events_df[all_events_df['geneSymbol'] == gene]
                if not gene_events.empty:
                    f.write(f"\n{gene}:\n")
                    for tissue in TISSUES:
                        tissue_events = gene_events[gene_events['tissue'] == tissue]
                        if not tissue_events.empty:
                            dpsi_values = tissue_events['IncLevelDifference'].tolist()
                            dpsi_str = ", ".join([f"{x:+.3f}" for x in dpsi_values])
                            f.write(f"  {tissue:15s}: {dpsi_str}\n")
                    
                    # Check consistency
                    all_dpsi = gene_events['IncLevelDifference'].tolist()
                    consistent = all(x > 0 for x in all_dpsi) or all(x < 0 for x in all_dpsi)
                    consistency_str = "Consistent direction" if consistent else "Mixed directions"
                    f.write(f"  → {consistency_str}\n")
        
        f.write("\n3. JUNCTION STATUS\n")
        f.write("-"*60 + "\n")
        if not junction_df.empty:
            for gene in CORE_GENES:
                gene_junc = junction_df[junction_df['Gene'] == gene]
                if not gene_junc.empty:
                    f.write(f"\n{gene}:\n")
                    status_counts = gene_junc['Status'].value_counts()
                    for status, count in status_counts.items():
                        f.write(f"  {status:12s}: {count} events\n")
    
    print(f"\n  Full report saved to: {summary_file}")
    
    # Save detailed tables
    gene_df.to_csv(os.path.join(OUTPUT_DIR, "gene_types.csv"), index=False)
    if not junction_df.empty:
        junction_df.to_csv(os.path.join(OUTPUT_DIR, "junction_status.csv"), index=False)
    if all_events:
        all_events_df.to_csv(os.path.join(OUTPUT_DIR, "all_core_gene_events.csv"), index=False)
    
    print("\n  Additional files:")
    print("    - gene_types.csv")
    print("    - junction_status.csv")
    print("    - all_core_gene_events.csv")
    
    return gene_df, junction_df

if __name__ == "__main__":
    analyze_core_genes()
