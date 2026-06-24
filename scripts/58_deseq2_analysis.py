#!/usr/bin/env python3
# ============================================================
# DEPRECATED — DO NOT USE FOR CANONICAL ANALYSIS
# ============================================================
# This script uses PyDESeq2 (Python port) which differs from the
# canonical DESeq2 v1.38.3 (R/Bioconductor) specified in NeuroSplice.md.
# It also reads from a local featureCounts output generated with
# the incorrect GENCODE vM25 annotation.
#
# USE INSTEAD: scripts/hpc/58_deseq2_hpc.R
#   (R/DESeq2 v1.38.3, reads HPC featureCounts output with Ensembl 102)
# ============================================================
"""
DESeq2 Analysis for GSE117357 (Adgrl3 KO)
=========================================
Performs differential gene expression analysis across 3 brain regions.
Uses PyDESeq2 to fit negative binomial GLMs and compute log2 fold changes.

Key Decisions:
- Design Formula: `~Genotype` (Sex covariate dropped because 100% of samples are male).
- Pre-filtering: Retain genes with ≥ 10 counts in at least N samples (N = smallest group size, i.e., 10).
- Significance Thresholds: padj < 0.05 AND |log2FoldChange| > 0.5.

Outputs:
- Full DESeq2 stats matrix (baseMean, log2FoldChange, padj) -> needed for downstream background definition.
- Filtered significant DEGs.
"""

import os
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

print("Loading Data...")
# Paths
PROJECT = os.path.expanduser("~/Desktop/WXY/NeuroSplice")
COUNTS_FILE = os.path.join(PROJECT, "results/featurecounts/gene_counts.txt")
METADATA_FILE = os.path.join(PROJECT, "data/metadata/GSE117357_complete_mapping.csv")
OUT_DIR = os.path.join(PROJECT, "results/deseq2")
os.makedirs(OUT_DIR, exist_ok=True)

# Determine CPUs available, defaulting to 8 if not found or > 8
n_cpus = min(os.cpu_count() or 8, 8)

# 1. Load Count Matrix
counts_df = pd.read_csv(COUNTS_FILE, sep="\t", comment="#", index_col=0)

# Clean up count matrix columns (strip BAM paths to just SRR_ID)
counts_df = counts_df.drop(columns=["Chr", "Start", "End", "Strand", "Length"])
counts_df.columns = [col.split("/")[-1].replace(".Aligned.sortedByCoord.out.bam", "") for col in counts_df.columns]

# pydeseq2 expects samples as rows and genes as columns
counts_df = counts_df.T

# 2. Load Metadata
metadata = pd.read_csv(METADATA_FILE)
metadata = metadata.set_index("SRR_ID")

tissues = ["Hippocampus", "Prefrontal-Cortex", "Striatum"]

# 3. Run DESeq2 per tissue
for tissue in tissues:
    print(f"\n========================================================")
    print(f"Running DESeq2 for: {tissue}")
    print(f"========================================================")
    
    # Subset metadata and counts
    tissue_meta = metadata[metadata["Tissue"] == tissue]
    
    if len(tissue_meta) == 0:
        print(f"WARNING: No samples found in metadata for tissue '{tissue}'. Skipping.")
        continue
        
    common = tissue_meta.index.intersection(counts_df.index)
    if len(common) < len(tissue_meta):
        missing = tissue_meta.index.difference(counts_df.index)
        print(f"WARNING: {len(missing)} sample(s) missing from count matrix: {list(missing)}")
        
    if len(common) == 0:
        print(f"ERROR: No matching samples in count matrix for {tissue}. Skipping.")
        continue
        
    tissue_meta = tissue_meta.loc[common]
    tissue_counts = counts_df.loc[common]
    
    # Dynamically compute sample group sizes
    n_wt = sum(tissue_meta["Genotype"] == "WT")
    n_ko = sum(tissue_meta["Genotype"] == "KO")
    
    if n_wt == 0 or n_ko == 0:
        print(f"ERROR: {tissue} is missing one genotype group after intersection (WT={n_wt}, KO={n_ko}). Skipping.")
        continue
        
    min_group_size = min(n_wt, n_ko)
    
    # Strict Pre-filtering: keep genes with >= 10 counts in at least `min_group_size` samples
    mask = (tissue_counts >= 10).sum(axis=0) >= min_group_size
    genes_to_keep = tissue_counts.columns[mask]
    tissue_counts = tissue_counts[genes_to_keep]
    
    print(f"Samples: {len(tissue_meta)} ({n_wt} WT, {n_ko} KO)")
    print(f"Genes after strict pre-filtering: {len(genes_to_keep)}")
    
    # Create DESeqDataSet object
    dds = DeseqDataSet(
        counts=tissue_counts,
        metadata=tissue_meta,
        design_factors="Genotype",
        refit_cooks=True,
        n_cpus=n_cpus
    )
    
    # Fit size factors and dispersion estimates
    print("Fitting model...")
    dds.deseq2()
    
    # Statistical testing: KO vs WT (WT is the reference)
    print("Testing KO vs WT...")
    stat_res = DeseqStats(dds, contrast=["Genotype", "KO", "WT"])
    
    res_df = stat_res.results_df
    
    # Save raw full results
    out_file = os.path.join(OUT_DIR, f"{tissue}_DESeq2_full.csv")
    res_df.to_csv(out_file)
    print(f"Full results saved: {out_file}")
    
    # Extract significant DEGs
    # Thresholds: padj < 0.05 and |log2FoldChange| > 0.5
    # Note: DESeq2 assigns NaN padj to outlier/low-count genes; pandas drops NaNs automatically here.
    sig_res = res_df[(res_df["padj"] < 0.05) & (res_df["log2FoldChange"].abs() > 0.5)]
    sig_file = os.path.join(OUT_DIR, f"{tissue}_DESeq2_significant.csv")
    sig_res.to_csv(sig_file)
    print(f"Significant DEGs found: {len(sig_res)}")
    print(f"Significant results saved: {sig_file}")

print("\nDone! Processing loop complete. Check output above for any warnings.")
