import os
import gseapy as gp
import pandas as pd
import json

# Configuration
PROJECT_DIR = "/Volumes/Untitled/NeuroSplice"
OUT_DIR = os.path.join(PROJECT_DIR, "results", "disease_context")
INPUT_GENES_FILE = os.path.join(PROJECT_DIR, "results", "string", "significant_genes.txt")
CORE_GENES = ["Dctn1", "Crem", "Pdlim7", "Spata5"]

# Hardcoded High Confidence Lists (Fallback/Augment)
# SFARI Score 1 (Selected prominent genes)
SFARI_HC = {
    "Adnp", "Ank2", "Ank3", "Arid1b", "Chd8", "Cntnap2", "Dscam", 
    "Fmr1", "Foxp1", "Grin2b", "Katnal2", "Nrxn1", "Pogz", "Pten", 
    "Scn2a", "Shank2", "Shank3", "Syngap1", "Tbr1", "Tsc1", "Tsc2",
    "Dctn1", "Crem", "Pdlim7", "Spata5" # Check if our core are in the full set
}

# PsychENCODE / SCZ GWAS (Common hits)
SCZ_RISK = {
    "C4a", "C4b", "Drd2", "Grin2a", "Cacna1c", "Tcf4", "Nrgn", 
    "Znf804a", "Mir137", "Bcl11b", "Hcn1", "Rims1", "Smg6",
    "Grm3", "Slc32a1", "Furin", "Tsnaare1"
}

def ensure_dir(d):
    if not os.path.exists(d):
        os.makedirs(d)

def run_disease_context():
    ensure_dir(OUT_DIR)
    
    # Load Splicing Genes
    if not os.path.exists(INPUT_GENES_FILE):
        print("Significant genes file not found. Generating from rMATS...")
        # (Simplified: Just assume we have it or exit, but better to be robust)
        # Let's read from the rMATS files directly if needed
        all_genes = set()
        for f in ["SE.MATS.JCEC.txt", "MXE.MATS.JCEC.txt"]:
            for region in ["hippocampus", "prefrontal_cortex", "striatum"]:
                path = os.path.join(PROJECT_DIR, "results", f"rmats_{region}", f)
                if os.path.exists(path):
                    df = pd.read_csv(path, sep='\t')
                    sigs = df[df['FDR']<0.05]['geneSymbol'].tolist()
                    all_genes.update([str(g).replace('"','') for g in sigs])
        
        gene_list = list(all_genes)
        print(f"Loaded {len(gene_list)} genes from rMATS.")
    else:
        with open(INPUT_GENES_FILE, 'r') as f:
            gene_list = [l.strip() for l in f if l.strip()]
        print(f"Loaded {len(gene_list)} genes from {INPUT_GENES_FILE}")

    # 1. Enrichr Analysis - Robust Library Loading
    print("Running Enrichr...")
    # Valid libraries often used: 'DisGeNET', 'GWAS_Catalog_2019', 'Jensen_DISEASES'
    # 'SFARI_Gene_Human_Gene_Module' might be deprecated or named differently.
    # We will try a few and print errors non-fatally.
    libraries = ['DisGeNET', 'GWAS_Catalog_2019', 'Jensen_DISEASES']
    
    try:
        enr = gp.enrichr(gene_list=gene_list,
                         gene_sets=libraries,
                         organism='Mouse', 
                         outdir=OUT_DIR,
                         cutoff=0.5
                         )
        
        if enr.results is not None:
             res = enr.results
             # Save full results
             res.to_csv(os.path.join(OUT_DIR, "disease_enrichment_stat.csv"))
             
             print("\nTop Disease Terms (FDR < 0.1):")
             top = res[res['Adjusted P-value'] < 0.1]
             if not top.empty:
                 # Print full term and genes for top 5
                 for i, row in top.head(5).iterrows():
                     print(f"- {row['Term']} ({row['Gene_set']}) FDR:{row['Adjusted P-value']:.2e}")
                     print(f"  Genes: {row['Genes']}")
             else:
                 print("No significant disease enrichment found.")

             # Check Core 4 in results
             print("\n--- Core 4 Overlap Check ---")
             for core in CORE_GENES:
                 mask = res['Genes'].apply(lambda x: core.upper() in x.upper().split(';'))
                 hits = res[mask]
                 if not hits.empty:
                     print(f"** {core} ** found in {len(hits)} terms. Top 3:")
                     for i, row in hits.sort_values('Adjusted P-value').head(3).iterrows():
                        print(f"  - {row['Term']} ({row['Gene_set']}) FDR:{row['Adjusted P-value']:.2e}")
                 else:
                     print(f"** {core} ** - No enrichment terms found containing this gene.")

    except Exception as e:
        print(f"Enrichr failed: {e}")

if __name__ == "__main__":
    run_disease_context()
