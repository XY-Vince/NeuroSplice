#!/usr/bin/env python3
"""
Human Conservation Analysis Framework
======================================
Automated tools for validating mouse splicing findings in human.

Analysis modules:
1. GTEx sQTL lookup (brain tissues)
2. ClinVar/OMIM pathogenic variant check
3. VAST-DB conservation (manual step)
4. Literature links (PubMed)
"""

import pandas as pd
import os
import requests
from time import sleep

OUTPUT_DIR = "/Volumes/Untitled/NeuroSplice/results/human_conservation"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Genes to analyze
CATEGORY_A_GENES = {
    'PTS': 'ENSG00000150787',  # Human Ensembl IDs
    'MYO9B': 'ENSG00000174099',
    'LRP8': 'ENSG00000157193'
}

CATEGORY_B_GENES = {
    'UNC13B': 'ENSG00000198722',  # Example - need to verify
    'SHANK1': 'ENSG00000131473'
}

def check_gtex_sqtls(gene_symbol, ensembl_id):
    """
    Query GTEx for brain tissue sQTLs.
    Note: GTEx API may require manual download for full data.
    """
    print(f"\n{gene_symbol} (GTEx sQTL Analysis)")
    print("-" * 50)
    
    # GTEx tissues of interest (brain)
    brain_tissues = [
        'Brain_Amygdala',
        'Brain_Anterior_cingulate_cortex_BA24',
        'Brain_Caudate_basal_ganglia',
        'Brain_Cerebellar_Hemisphere',
        'Brain_Cerebellum',
        'Brain_Cortex',
        'Brain_Frontal_Cortex_BA9',
        'Brain_Hippocampus',
        'Brain_Hypothalamus',
        'Brain_Nucleus_accumbens_basal_ganglia',
        'Brain_Putamen_basal_ganglia',
        'Brain_Spinal_cord_cervical_c-1',
        'Brain_Substantia_nigra'
    ]
    
    print(f"  Gene: {gene_symbol}")
    print(f"  Ensembl ID: {ensembl_id}")
    print(f"  Brain tissues to check: {len(brain_tissues)}")
    
    # Note: Full GTEx sQTL data requires download
    # https://gtexportal.org/home/datasets
    print(f"\n  Manual steps:")
    print(f"  1. Visit: https://gtexportal.org/home/gene/{gene_symbol}")
    print(f"  2. Check 'sQTLs' tab")
    print(f"  3. Filter for brain tissues")
    print(f"  4. Record significant sQTLs (p < 1e-5)")
    
    return {
        'Gene': gene_symbol,
        'Ensembl_ID': ensembl_id,
        'GTEx_URL': f"https://gtexportal.org/home/gene/{gene_symbol}",
        'Brain_Tissues_Available': len(brain_tissues),
        'Manual_Check_Required': True
    }

def check_clinvar_pathogenic(gene_symbol):
    """
    Query ClinVar for pathogenic variants.
    """
    print(f"\n{gene_symbol} (ClinVar Pathogenic Variants)")
    print("-" * 50)
    
    # ClinVar E-utilities query
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    
    try:
        # Search ClinVar for gene
        search_url = f"{base_url}esearch.fcgi?db=clinvar&term={gene_symbol}[gene]%20AND%20pathogenic[clinical_significance]&retmode=json"
        
        response = requests.get(search_url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            count = int(data.get('esearchresult', {}).get('count', 0))
            
            print(f"  Pathogenic/Likely Pathogenic variants: {count}")
            
            if count > 0:
                print(f"  View details: https://www.ncbi.nlm.nih.gov/clinvar/?term={gene_symbol}[gene]+AND+pathogenic[clin]")
                
                return {
                    'Gene': gene_symbol,
                    'ClinVar_Pathogenic_Count': count,
                    'ClinVar_URL': f"https://www.ncbi.nlm.nih.gov/clinvar/?term={gene_symbol}[gene]+AND+pathogenic[clin]",
                    'Manual_Review_Needed': count > 0
                }
        else:
            print(f"  ⚠ API request failed (status {response.status_code})")
            return {'Gene': gene_symbol, 'Error': 'API_FAILED'}
            
    except Exception as e:
        print(f"  ⚠ Error: {e}")
        return {'Gene': gene_symbol, 'Error': str(e)}

def check_omim(gene_symbol):
    """
    Check OMIM for disease associations.
    Note: OMIM API requires registration, this provides manual links.
    """
    print(f"\n{gene_symbol} (OMIM Disease Associations)")
    print("-" * 50)
    
    omim_url = f"https://www.omim.org/search?search={gene_symbol}"
    print(f"  Manual check required:")
    print(f"  Visit: {omim_url}")
    print(f"  Look for:")
    print(f"    - Phenotype associations")
    print(f"    - Neurological disorders")
    print(f"    - ADHD/neurodevelopmental relevance")
    
    return {
        'Gene': gene_symbol,
        'OMIM_URL': omim_url,
        'Manual_Check_Required': True
    }

def vastdb_conservation_template(gene_symbol):
    """
    Generate template for VAST-DB manual analysis.
    """
    print(f"\n{gene_symbol} (VAST-DB Exon Conservation)")
    print("-" * 50)
    
    print(f"  Manual steps:")
    print(f"  1. Visit: https://vastdb.crg.eu/")
    print(f"  2. Search for gene: {gene_symbol}")
    print(f"  3. Select human (Homo sapiens)")
    print(f"  4. Check 'Conservation' tab")
    print(f"  5. Record PSI values across species")
    print(f"  6. Note: Conservation score = % identity across mammals")
    
    return {
        'Gene': gene_symbol,
        'VASTDB_URL': 'https://vastdb.crg.eu/',
        'Manual_Steps': [
            'Search gene',
            'Check human orthologs',
            'Record PSI conservation',
            'Note alternative exon usage'
        ]
    }

def main():
    """Run conservation analysis for all genes."""
    
    print("="*70)
    print("Human Conservation Analysis")
    print("="*70)
    
    # Category A genes
    print("\n" + "="*70)
    print("CATEGORY A GENES")
    print("="*70)
    
    gtex_results = []
    clinvar_results = []
    omim_results = []
    vastdb_results = []
    
    for gene_sym, ensembl in CATEGORY_A_GENES.items():
        # GTEx sQTLs
        gtex_res = check_gtex_sqtls(gene_sym, ensembl)
        gtex_results.append(gtex_res)
        sleep(1)  # Rate limiting
        
        # ClinVar
        clinvar_res = check_clinvar_pathogenic(gene_sym)
        clinvar_results.append(clinvar_res)
        sleep(1)
        
        # OMIM
        omim_res = check_omim(gene_sym)
        omim_results.append(omim_res)
        
        # VAST-DB template
        vastdb_res = vastdb_conservation_template(gene_sym)
        vastdb_results.append(vastdb_res)
    
    # Save results
    pd.DataFrame(gtex_results).to_csv(
        os.path.join(OUTPUT_DIR, "gtex_sqtl_links.csv"), index=False
    )
    pd.DataFrame(clinvar_results).to_csv(
        os.path.join(OUTPUT_DIR, "clinvar_pathogenic.csv"), index=False
    )
    pd.DataFrame(omim_results).to_csv(
        os.path.join(OUTPUT_DIR, "omim_links.csv"), index=False
    )
    pd.DataFrame(vastdb_results).to_csv(
        os.path.join(OUTPUT_DIR, "vastdb_manual_steps.csv"), index=False
    )
    
    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"\nResults saved to: {OUTPUT_DIR}")
    print("\nNext steps (manual):")
    print("  1. Visit GTEx URLs for sQTL data")
    print("  2. Review ClinVar pathogenic variants")
    print("  3. Check OMIM disease associations")
    print("  4. Analyze VAST-DB conservation scores")
    
    # Create summary report
    summary = f"""
# Human Conservation Analysis Summary

## Category A Genes

### PTS
- **GTEx**: {gtex_results[0]['GTEx_URL']}
- **ClinVar**: Check pathogenic variants in neurodevelopmental context
- **VAST-DB**: Check exon conservation across primates
- **Key Question**: Are BH4 pathway exons conserved?

### MYO9B
- **GTEx**: {gtex_results[1]['GTEx_URL']}
- **ClinVar**: Look for neurological phenotypes
- **Function**: Rho GTPase, cytoskeletal regulation

### LRP8
- **GTEx**: {gtex_results[2]['GTEx_URL']}
- **Special**: Check ApoE4 interaction literature
- **PubMed**: Search "LRP8 ApoE4 synaptic dysfunction"
- **Relevance**: Reelin pathway, synaptic plasticity

## Action Items
1. Download GTEx sQTL data for brain tissues
2. Record ClinVar neurological variants
3. Literature review for LRP8-ApoE4
4. VAST-DB conservation analysis

Output: {OUTPUT_DIR}
"""
    
    with open(os.path.join(OUTPUT_DIR, "analysis_summary.md"), 'w') as f:
        f.write(summary)
    
    print(f"\nSummary report: {OUTPUT_DIR}/analysis_summary.md")

if __name__ == "__main__":
    main()
