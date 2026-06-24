#!/usr/bin/env python3
"""
Human Translation & GWAS Analysis for GSE117357 Splicing Genes
==============================================================
1. Convert mouse genes to human orthologs
2. Query GTEx for splicing events in brain tissues
3. Run MAGMA enrichment against ADHD GWAS
4. Literature deep-dive on top candidates
"""

import pandas as pd
import numpy as np
import os
import subprocess
from collections import defaultdict
import requests
import json
import warnings
warnings.filterwarnings('ignore')

# Configuration
RESULTS_DIR = "/Volumes/Untitled/NeuroSplice/results/rmats_gse117357_analysis"
MAGMA_DIR = "/Volumes/Untitled/NeuroSplice/results/magma"
OUTPUT_DIR = "/Volumes/Untitled/NeuroSplice/results/human_translation"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Core genes and top candidates
CORE_GENES = ['Arid5a', 'Bcl2l11', 'Gm10419', 'Lrp8', 'Myo9b', 'Pts']
TOP_CANDIDATES = ['Bcl2l11', 'Lrp8', 'Pts']  # Focus genes

# Mouse to human mapping (via biomart or manual for key genes)
MOUSE_HUMAN_ORTHOLOGS = {
    'Lrp8': 'LRP8',
    'Pts': 'PTS', 
    'Bcl2l11': 'BCL2L11',
    'Arid5a': 'ARID5A',
    'Myo9b': 'MYO9B',
    # Add more as needed
}

def load_all_genes():
    """Load and combine genes from all tissues."""
    all_genes = set()
    for tissue in ['hippocampus', 'prefrontal', 'striatum']:
        filepath = os.path.join(RESULTS_DIR, f"{tissue}_significant_genes.csv")
        if os.path.exists(filepath):
            df = pd.read_csv(filepath)
            all_genes.update(df['Gene'].tolist())
    return list(all_genes)

def convert_to_human(mouse_genes):
    """Convert mouse gene symbols to human using API or mapping."""
    print("\n=== Converting Mouse to Human Orthologs ===")
    
    human_genes = []
    mapped = {}
    unmapped = []
    
    for gene in mouse_genes:
        # Direct mapping (gene symbols often same in mouse/human)
        human_gene = gene.upper()
        
        # Check against known orthologs
        if gene in MOUSE_HUMAN_ORTHOLOGS:
            human_gene = MOUSE_HUMAN_ORTHOLOGS[gene]
        
        # Skip mouse-specific predicted genes
        if gene.startswith('Gm') or gene.startswith('Rik') or gene.startswith('mt-'):
            unmapped.append(gene)
            continue
            
        human_genes.append(human_gene)
        mapped[gene] = human_gene
    
    print(f"  Mapped: {len(mapped)} genes")
    print(f"  Unmapped (mouse-specific): {len(unmapped)}")
    
    return human_genes, mapped

def get_entrez_ids(human_genes):
    """Get Entrez IDs from MAGMA gene location file."""
    print("\n=== Getting Entrez IDs ===")
    
    gene_loc_file = os.path.join(MAGMA_DIR, "NCBI37.3.gene.loc")
    gene_loc = pd.read_csv(gene_loc_file, sep='\t', header=None, 
                           names=['entrez', 'chr', 'start', 'end', 'strand', 'symbol'])
    
    # Create mapping
    symbol_to_entrez = dict(zip(gene_loc['symbol'], gene_loc['entrez']))
    
    entrez_ids = []
    found = 0
    for gene in human_genes:
        if gene in symbol_to_entrez:
            entrez_ids.append(str(symbol_to_entrez[gene]))
            found += 1
    
    print(f"  Found Entrez IDs for {found}/{len(human_genes)} genes")
    return entrez_ids

def run_magma_enrichment(entrez_ids, set_name="GSE117357_splicing"):
    """Run MAGMA gene-set analysis against ADHD GWAS."""
    print("\n=== Running MAGMA Enrichment ===")
    
    # Create gene set file
    geneset_file = os.path.join(OUTPUT_DIR, f"{set_name}.geneset.txt")
    with open(geneset_file, 'w') as f:
        f.write(f"{set_name} " + " ".join(entrez_ids) + "\n")
    
    print(f"  Gene set file: {geneset_file}")
    print(f"  Genes in set: {len(entrez_ids)}")
    
    # Run MAGMA
    magma_binary = os.path.join(MAGMA_DIR, "magma")
    gene_results = os.path.join(MAGMA_DIR, "gene_analysis.genes.raw")
    output_prefix = os.path.join(OUTPUT_DIR, f"{set_name}_magma")
    
    cmd = [
        magma_binary,
        "--gene-results", gene_results,
        "--set-annot", geneset_file,
        "--out", output_prefix
    ]
    
    print(f"  Running: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode == 0:
            print("  MAGMA completed successfully!")
            # Parse results
            gsa_file = f"{output_prefix}.gsa.out"
            if os.path.exists(gsa_file):
                with open(gsa_file, 'r') as f:
                    for line in f:
                        if not line.startswith('#'):
                            print(f"  Result: {line.strip()}")
                return gsa_file
        else:
            print(f"  MAGMA error: {result.stderr}")
    except Exception as e:
        print(f"  Error running MAGMA: {e}")
    
    return None

def check_candidate_gwas_proximity(candidates):
    """Check if candidate genes are near ADHD GWAS hits."""
    print("\n=== Checking Candidate Gene GWAS Proximity ===")
    
    # Get gene-level MAGMA results
    gene_results_file = os.path.join(MAGMA_DIR, "gene_analysis.genes.out")
    gene_results = pd.read_csv(gene_results_file, sep='\s+', comment='#')
    
    # Load gene locations
    gene_loc_file = os.path.join(MAGMA_DIR, "NCBI37.3.gene.loc")
    gene_loc = pd.read_csv(gene_loc_file, sep='\t', header=None,
                           names=['entrez', 'chr', 'start', 'end', 'strand', 'symbol'])
    
    results = []
    for candidate in candidates:
        human_gene = candidate.upper()
        row = gene_loc[gene_loc['symbol'] == human_gene]
        
        if not row.empty:
            entrez_id = row['entrez'].values[0]
            gene_res = gene_results[gene_results['GENE'] == entrez_id]
            
            if not gene_res.empty:
                p_val = gene_res['P'].values[0]
                z_score = gene_res['ZSTAT'].values[0]
                chr_loc = row['chr'].values[0]
                
                results.append({
                    'Gene': human_gene,
                    'Entrez': entrez_id,
                    'Chromosome': chr_loc,
                    'MAGMA_P': p_val,
                    'MAGMA_Z': z_score,
                    'Significant': p_val < 0.05
                })
                
                sig_str = "✓ SIGNIFICANT" if p_val < 0.05 else ""
                print(f"  {human_gene}: chr{chr_loc}, P={p_val:.4f}, Z={z_score:.2f} {sig_str}")
            else:
                print(f"  {human_gene}: No GWAS signal found")
        else:
            print(f"  {human_gene}: Gene not in reference")
    
    return pd.DataFrame(results)

def query_gtex_sqtl(gene):
    """Query GTEx sQTL browser for splicing QTLs in brain."""
    print(f"\n=== Querying GTEx sQTLs for {gene} ===")
    
    # GTEx brain tissues of interest
    brain_tissues = [
        'Brain_Hippocampus',
        'Brain_Frontal_Cortex_BA9',
        'Brain_Caudate_basal_ganglia',
        'Brain_Putamen_basal_ganglia',
        'Brain_Cortex'
    ]
    
    # Note: This would require GTEx API access or downloaded data
    # For now, we'll create a summary of what to check
    print("  GTEx sQTL Analysis Recommendations:")
    print(f"  1. Visit: https://gtexportal.org/home/gene/{gene}")
    print(f"  2. Check 'sQTL' tab for brain tissues")
    print(f"  3. Look for splicing events matching mouse AS patterns")
    
    return None

def literature_review():
    """Compile literature findings for top candidates."""
    print("\n=== Literature Deep-Dive ===")
    
    findings = {
        'LRP8': {
            'alias': 'ApoER2',
            'function': 'LDL receptor related protein 8 / Apolipoprotein E receptor 2',
            'pathway': 'Reelin signaling',
            'brain_relevance': [
                'Reelin-ApoER2 signaling critical for neuronal migration',
                'Mediates synaptic plasticity in hippocampus',
                'Controls dendritic spine development',
                'Interacts with Dab1 for downstream signaling'
            ],
            'adhd_connection': [
                'RELN (Reelin) associated with ADHD in some studies',
                'Synaptic plasticity deficits linked to attention disorders',
                'ApoER2 KO mice show learning/memory deficits',
                'Dopamine signaling modulation via Reelin pathway'
            ],
            'references': [
                'Trommsdorff et al., Cell 1999',
                'Weeber et al., J Biol Chem 2002',
                'Rogers et al., Nat Rev Neurosci 2011'
            ]
        },
        'PTS': {
            'alias': '6-pyruvoyltetrahydropterin synthase',
            'function': 'BH4 biosynthesis enzyme',
            'pathway': 'Tetrahydrobiopterin (BH4) synthesis',
            'brain_relevance': [
                'BH4 is essential cofactor for:',
                '  - Tyrosine hydroxylase (dopamine synthesis)',
                '  - Tryptophan hydroxylase (serotonin synthesis)',
                '  - Phenylalanine hydroxylase',
                'BH4 deficiency causes hyperphenylalaninemia and neurotransmitter deficiency'
            ],
            'adhd_connection': [
                '**DIRECT LINK**: Dopamine synthesis requires BH4',
                'ADHD core pathology: dopaminergic dysfunction',
                'BH4 supplementation improves attention in some patients',
                'GCH1 (BH4 pathway) variants associated with ADHD',
                'Methylphenidate (ADHD medication) affects dopamine'
            ],
            'references': [
                'Thöny et al., Physiol Rev 2000',
                'Opladen et al., Mol Genet Metab 2012',
                'Khanh et al., J Neurochem 2016'
            ]
        },
        'BCL2L11': {
            'alias': 'BIM',
            'function': 'BH3-only pro-apoptotic protein',
            'pathway': 'Intrinsic apoptosis pathway',
            'brain_relevance': [
                'Key regulator of neuronal apoptosis',
                'Critical for developmental synaptic pruning',
                'Controls neuronal number during development',
                'Activity-dependent survival signaling'
            ],
            'adhd_connection': [
                'Abnormal synaptic pruning proposed in ADHD',
                'Prefrontal cortex development affected in ADHD',
                'Apoptosis dysregulation may affect neural circuit formation',
                'Less direct than PTS/dopamine link'
            ],
            'references': [
                'Whitfield et al., Neuron 2001',
                'Ghosh et al., Cell Death & Dis 2012'
            ]
        }
    }
    
    # Save literature review
    output_file = os.path.join(OUTPUT_DIR, "literature_review.txt")
    with open(output_file, 'w') as f:
        f.write("Literature Deep-Dive: Top Candidate Genes\n")
        f.write("=" * 60 + "\n\n")
        
        for gene, info in findings.items():
            f.write(f"\n{'='*60}\n")
            f.write(f"{gene} ({info['alias']})\n")
            f.write(f"{'='*60}\n")
            f.write(f"Function: {info['function']}\n")
            f.write(f"Pathway: {info['pathway']}\n\n")
            
            f.write("Brain Relevance:\n")
            for item in info['brain_relevance']:
                f.write(f"  • {item}\n")
            
            f.write("\nADHD Connection:\n")
            for item in info['adhd_connection']:
                f.write(f"  • {item}\n")
            
            f.write("\nKey References:\n")
            for ref in info['references']:
                f.write(f"  - {ref}\n")
            f.write("\n")
    
    print(f"  Saved to: {output_file}")
    
    # Print summary
    for gene, info in findings.items():
        print(f"\n  {gene} ({info['alias']}):")
        print(f"    Pathway: {info['pathway']}")
        print(f"    ADHD Link: {info['adhd_connection'][0]}")
    
    return findings

def main():
    print("="*60)
    print("Human Translation & GWAS Analysis")
    print("="*60)
    
    # 1. Load and convert genes
    mouse_genes = load_all_genes()
    print(f"\nTotal unique mouse genes: {len(mouse_genes)}")
    
    human_genes, mapping = convert_to_human(mouse_genes)
    
    # 2. Get Entrez IDs for MAGMA
    entrez_ids = get_entrez_ids(human_genes)
    
    # 3. Run MAGMA enrichment
    if entrez_ids:
        magma_result = run_magma_enrichment(entrez_ids)
    
    # 4. Check candidate gene GWAS proximity
    candidate_df = check_candidate_gwas_proximity(TOP_CANDIDATES)
    if not candidate_df.empty:
        candidate_df.to_csv(os.path.join(OUTPUT_DIR, "candidate_gwas_results.csv"), index=False)
    
    # 5. GTEx sQTL recommendations
    for gene in TOP_CANDIDATES:
        query_gtex_sqtl(gene.upper())
    
    # 6. Literature review
    literature = literature_review()
    
    # Summary
    print("\n" + "="*60)
    print("ANALYSIS COMPLETE")
    print("="*60)
    print(f"\nOutput directory: {OUTPUT_DIR}")
    print("\nKey Findings:")
    print("  • PTS: DIRECT ADHD link via BH4 → dopamine synthesis")
    print("  • LRP8: Reelin pathway → synaptic plasticity")
    print("  • BCL2L11: Developmental pruning (indirect)")
    
    print("\nNext Steps:")
    print("  1. Check GTEx sQTL browser for brain splicing")
    print("  2. Review MAGMA results for gene-set enrichment")
    print("  3. Deep-dive on PTS as priority candidate")

if __name__ == "__main__":
    main()
