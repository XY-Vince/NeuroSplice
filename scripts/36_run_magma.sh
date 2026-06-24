#!/bin/bash

# MAGMA Analysis Script for ADHD Splicing Genes

# PREREQUISITES:
# 1. Unblock MAGMA binary: Run `results/magma/magma --help` manually and approve security prompt.
# 2. Download Reference Data: 
#    - Download `g1000_eur.zip` from https://ctg.cncr.nl/software/magma
#    - Unzip into `results/magma/` (should contain g1000_eur.bed/.bim/.fam)

echo "--- Step 1: Annotation ---"
# Maps SNPs to Genes using NCBI37.3 coordinates
# Uses the SNP location file derived from ADHD 2022 SumStats
results/magma/magma --annotate window=0 --snp-loc results/magma/snp_loc.txt --gene-loc results/magma/NCBI37.3.gene.loc --out results/magma/annotation

echo "--- Step 2: Gene Analysis ---"
# Calculates Gene-Level P-values from SNP P-values (Correcting for LD)
# NOTE: This step REQUIRES the g1000_eur reference files. 
# If missing, it will fail or run in 'snp-wise-mean' mode (invalid for GWAS).
results/magma/magma --bfile results/magma/g1000_eur --pval results/magma/pval_file.txt ncol=3 --gene-annot results/magma/annotation.genes.annot --out results/magma/gene_analysis

echo "--- Step 3: Gene Set Analysis ---"
# Tests enrichment of our 886 splicing genes against the ADHD background
results/magma/magma --gene-results results/magma/gene_analysis.genes.raw --set-annot results/magma/geneset.txt --out results/magma/set_analysis

echo "Done. Check results/magma/set_analysis.gsa.out for P-value."
