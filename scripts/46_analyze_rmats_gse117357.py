#!/usr/bin/env python3
"""
GSE117357 Differential Splicing Analysis
=========================================
Comprehensive analysis of rMATS results:
1. Filter for FDR < 0.05 and |ΔΨ| > 0.1
2. Identify genes with differential splicing
3. Compare splicing patterns across tissues
4. Functional enrichment analysis (using gprofiler-official)
"""

import pandas as pd
import numpy as np
import os
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Configuration
RESULTS_DIR = "/Volumes/Untitled/NeuroSplice/results/rmats_gse117357"
OUTPUT_DIR = "/Volumes/Untitled/NeuroSplice/results/rmats_gse117357_analysis"
TISSUES = ["hippocampus", "prefrontal", "striatum"]
EVENT_TYPES = ["SE", "A5SS", "A3SS", "MXE", "RI"]
FDR_THRESHOLD = 0.05
DELTA_PSI_THRESHOLD = 0.1

os.makedirs(OUTPUT_DIR, exist_ok=True)

def load_rmats_results(tissue, event_type):
    """Load rMATS JCEC results for a tissue and event type."""
    filepath = os.path.join(RESULTS_DIR, tissue, f"{event_type}.MATS.JCEC.txt")
    if not os.path.exists(filepath):
        print(f"Warning: {filepath} not found")
        return pd.DataFrame()
    
    df = pd.read_csv(filepath, sep="\t")
    df['tissue'] = tissue
    df['event_type'] = event_type
    
    # Clean gene symbol (remove quotes)
    if 'geneSymbol' in df.columns:
        df['geneSymbol'] = df['geneSymbol'].str.replace('"', '')
    if 'GeneID' in df.columns:
        df['GeneID'] = df['GeneID'].str.replace('"', '')
    
    return df

def filter_significant_events(df, fdr_threshold=FDR_THRESHOLD, dpsi_threshold=DELTA_PSI_THRESHOLD):
    """Filter for significant differential splicing events."""
    if df.empty:
        return df
    
    # Filter by FDR and |ΔΨ|
    significant = df[
        (df['FDR'] < fdr_threshold) & 
        (df['IncLevelDifference'].abs() > dpsi_threshold)
    ].copy()
    
    return significant

def analyze_tissue(tissue):
    """Analyze all event types for a tissue."""
    print(f"\n{'='*60}")
    print(f"Analyzing {tissue.upper()}")
    print('='*60)
    
    all_events = []
    summary = {}
    
    for event_type in EVENT_TYPES:
        df = load_rmats_results(tissue, event_type)
        if df.empty:
            continue
            
        total = len(df)
        sig = filter_significant_events(df)
        
        summary[event_type] = {
            'total': total,
            'significant': len(sig),
            'up_in_wt': len(sig[sig['IncLevelDifference'] > 0]),
            'up_in_ko': len(sig[sig['IncLevelDifference'] < 0])
        }
        
        print(f"  {event_type}: {len(sig)}/{total} significant (FDR<{FDR_THRESHOLD}, |ΔΨ|>{DELTA_PSI_THRESHOLD})")
        
        if not sig.empty:
            all_events.append(sig)
    
    if all_events:
        combined = pd.concat(all_events, ignore_index=True)
        return combined, summary
    return pd.DataFrame(), summary

def identify_genes(sig_events):
    """Identify unique genes with differential splicing."""
    if sig_events.empty:
        return pd.DataFrame()
    
    # Group by gene and count events
    gene_summary = sig_events.groupby('geneSymbol').agg({
        'event_type': lambda x: ','.join(sorted(set(x))),
        'FDR': 'min',
        'IncLevelDifference': ['mean', 'min', 'max', 'count']
    }).reset_index()
    
    # Flatten column names
    gene_summary.columns = ['Gene', 'EventTypes', 'MinFDR', 'MeanDPSI', 'MinDPSI', 'MaxDPSI', 'NumEvents']
    gene_summary = gene_summary.sort_values('NumEvents', ascending=False)
    
    return gene_summary

def compare_tissues(tissue_genes):
    """Compare differential splicing genes across tissues."""
    print("\n" + "="*60)
    print("CROSS-TISSUE COMPARISON")
    print("="*60)
    
    # Get gene sets for each tissue
    gene_sets = {tissue: set(genes['Gene'].tolist()) for tissue, genes in tissue_genes.items() if not genes.empty}
    
    if len(gene_sets) < 2:
        print("Not enough tissues with significant genes for comparison")
        return {}
    
    # Calculate overlaps
    all_tissues = list(gene_sets.keys())
    
    results = {
        'tissue_counts': {t: len(g) for t, g in gene_sets.items()},
        'overlaps': {}
    }
    
    # Pairwise overlaps
    for i, t1 in enumerate(all_tissues):
        for t2 in all_tissues[i+1:]:
            overlap = gene_sets[t1] & gene_sets[t2]
            results['overlaps'][f"{t1}_&_{t2}"] = len(overlap)
            print(f"  {t1} ∩ {t2}: {len(overlap)} genes")
    
    # Three-way overlap (if 3 tissues)
    if len(all_tissues) == 3:
        core = gene_sets[all_tissues[0]] & gene_sets[all_tissues[1]] & gene_sets[all_tissues[2]]
        results['core_genes'] = list(core)
        print(f"\n  Core genes (all 3 tissues): {len(core)}")
        if len(core) <= 30:
            print(f"  Genes: {', '.join(sorted(core))}")
    
    # Tissue-specific genes
    print("\n  Tissue-specific genes:")
    for tissue in all_tissues:
        others = set().union(*[g for t, g in gene_sets.items() if t != tissue])
        specific = gene_sets[tissue] - others
        results[f'{tissue}_specific'] = list(specific)
        print(f"    {tissue}: {len(specific)} unique genes")
    
    return results

def run_enrichment(genes, tissue_name, output_dir):
    """Run functional enrichment using gprofiler."""
    if len(genes) == 0:
        print(f"  No genes for enrichment in {tissue_name}")
        return pd.DataFrame()
    
    print(f"\n  Running enrichment for {tissue_name} ({len(genes)} genes)...")
    
    try:
        from gprofiler import GProfiler
        gp = GProfiler(return_dataframe=True)
        
        results = gp.profile(
            organism='mmusculus',
            query=genes,
            sources=['GO:BP', 'GO:MF', 'GO:CC', 'KEGG', 'REAC'],
            user_threshold=0.05,
            significance_threshold_method='fdr'
        )
        
        if not results.empty:
            # Save results
            output_file = os.path.join(output_dir, f"{tissue_name}_enrichment.csv")
            results.to_csv(output_file, index=False)
            
            # Print top results
            print(f"  Top enriched terms for {tissue_name}:")
            for _, row in results.head(10).iterrows():
                print(f"    {row['source']}: {row['name']} (p={row['p_value']:.2e}, genes={row['intersection_size']})")
            
            return results
        else:
            print(f"  No significant enrichment found for {tissue_name}")
            return pd.DataFrame()
            
    except ImportError:
        print("  gprofiler not installed. Installing...")
        os.system("pip install gprofiler-official")
        return run_enrichment(genes, tissue_name, output_dir)
    except Exception as e:
        print(f"  Enrichment error: {e}")
        return pd.DataFrame()

def create_summary_plots(tissue_data, comparison, output_dir):
    """Create summary visualizations."""
    
    # 1. Event type distribution by tissue
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Bar plot of significant events
    event_counts = []
    for tissue, (events, summary) in tissue_data.items():
        for event_type, counts in summary.items():
            event_counts.append({
                'Tissue': tissue.capitalize(),
                'Event Type': event_type,
                'Count': counts['significant']
            })
    
    if event_counts:
        df_counts = pd.DataFrame(event_counts)
        df_pivot = df_counts.pivot(index='Tissue', columns='Event Type', values='Count')
        df_pivot.plot(kind='bar', ax=axes[0], colormap='Set2')
        axes[0].set_title('Significant Splicing Events by Tissue\n(FDR<0.05, |ΔΨ|>0.1)')
        axes[0].set_ylabel('Number of Events')
        axes[0].legend(title='Event Type')
        axes[0].tick_params(axis='x', rotation=45)
    
    # 2. Venn-like comparison (simple bar for now)
    if 'tissue_counts' in comparison:
        tissues = list(comparison['tissue_counts'].keys())
        counts = list(comparison['tissue_counts'].values())
        
        colors = ['#ff6b6b', '#4ecdc4', '#45b7d1']
        axes[1].bar(tissues, counts, color=colors[:len(tissues)])
        axes[1].set_title('Genes with Differential Splicing by Tissue')
        axes[1].set_ylabel('Number of Genes')
        
        # Add core genes annotation if available
        if 'core_genes' in comparison:
            axes[1].axhline(y=len(comparison['core_genes']), color='gray', 
                          linestyle='--', label=f"Core genes: {len(comparison['core_genes'])}")
            axes[1].legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'splicing_summary.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\n  Saved summary plot: {output_dir}/splicing_summary.png")

def main():
    print("="*60)
    print("GSE117357 Differential Splicing Analysis")
    print("="*60)
    print(f"FDR threshold: {FDR_THRESHOLD}")
    print(f"|ΔΨ| threshold: {DELTA_PSI_THRESHOLD}")
    
    tissue_data = {}
    tissue_genes = {}
    all_enrichment = {}
    
    # 1. Analyze each tissue
    for tissue in TISSUES:
        events, summary = analyze_tissue(tissue)
        tissue_data[tissue] = (events, summary)
        
        # 2. Identify genes
        if not events.empty:
            genes = identify_genes(events)
            tissue_genes[tissue] = genes
            
            # Save gene list
            output_file = os.path.join(OUTPUT_DIR, f"{tissue}_significant_genes.csv")
            genes.to_csv(output_file, index=False)
            print(f"\n  Saved {len(genes)} genes to {output_file}")
            
            # Also save full significant events
            events_file = os.path.join(OUTPUT_DIR, f"{tissue}_significant_events.csv")
            events.to_csv(events_file, index=False)
    
    # 3. Compare across tissues
    comparison = compare_tissues(tissue_genes)
    
    # Save comparison results
    if comparison:
        with open(os.path.join(OUTPUT_DIR, "tissue_comparison.txt"), 'w') as f:
            f.write("Tissue Comparison Results\n")
            f.write("="*40 + "\n\n")
            
            f.write("Gene counts per tissue:\n")
            for tissue, count in comparison.get('tissue_counts', {}).items():
                f.write(f"  {tissue}: {count}\n")
            
            f.write("\nOverlaps:\n")
            for overlap, count in comparison.get('overlaps', {}).items():
                f.write(f"  {overlap}: {count}\n")
            
            if 'core_genes' in comparison:
                f.write(f"\nCore genes (all 3 tissues): {len(comparison['core_genes'])}\n")
                f.write(f"  {', '.join(sorted(comparison['core_genes']))}\n")
    
    # 4. Functional enrichment
    print("\n" + "="*60)
    print("FUNCTIONAL ENRICHMENT ANALYSIS")
    print("="*60)
    
    for tissue, genes in tissue_genes.items():
        if not genes.empty:
            gene_list = genes['Gene'].tolist()
            enrichment = run_enrichment(gene_list, tissue, OUTPUT_DIR)
            all_enrichment[tissue] = enrichment
    
    # Core genes enrichment
    if 'core_genes' in comparison and len(comparison['core_genes']) > 0:
        core_enrichment = run_enrichment(list(comparison['core_genes']), "core_genes", OUTPUT_DIR)
        all_enrichment['core'] = core_enrichment
    
    # 5. Create summary visualizations
    print("\n" + "="*60)
    print("CREATING VISUALIZATIONS")
    print("="*60)
    create_summary_plots(tissue_data, comparison, OUTPUT_DIR)
    
    # Final summary
    print("\n" + "="*60)
    print("ANALYSIS COMPLETE")
    print("="*60)
    print(f"Results saved to: {OUTPUT_DIR}")
    print("\nFiles generated:")
    for f in os.listdir(OUTPUT_DIR):
        print(f"  - {f}")

if __name__ == "__main__":
    main()
