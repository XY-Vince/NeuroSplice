#!/usr/bin/env python3
# ============================================================
# Canonical Local DESeq2 Execution (via PyDESeq2)
# ============================================================
# Generates canonical DEG results with full reproducibility,
# strict metadata auditing, and QC logging.
# ============================================================

import os
import sys
import json
import logging
import argparse
import hashlib
from datetime import datetime
import pandas as pd
import numpy as np
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

def sha256_file(path):
    """Calculate SHA256 hash of a file for reproducibility manifest."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def setup_directories(project_dir):
    """Create structured output directories."""
    base_out = os.path.join(project_dir, "results/deseq2_canonical")
    dirs = {
        "full_results": os.path.join(base_out, "full_results"),
        "qc": os.path.join(base_out, "qc"),
        "logs": os.path.join(base_out, "logs")
    }
    for d in dirs.values():
        os.makedirs(d, exist_ok=True)
    return base_out, dirs

def parse_args():
    parser = argparse.ArgumentParser(description="Run canonical PyDESeq2 analysis")
    parser.add_argument("--project-dir", required=True, help="Path to NeuroSplice project directory")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility")
    parser.add_argument("--expected-n-per-genotype", type=int, default=10, help="Expected N samples per WT/KO group per tissue")
    parser.add_argument("--lfc-threshold", type=float, default=0.5, help="Absolute log2FoldChange threshold for summary reporting")
    return parser.parse_args()

def main():
    args = parse_args()
    project_dir = args.project_dir
    
    try:
        base_out, dirs = setup_directories(project_dir)
    except OSError as e:
        print(f"FATAL: Cannot create output directories: {e}", file=sys.stderr)
        sys.exit(1)

    # ---------------------------------------------------------
    # 1. Setup Logging
    # ---------------------------------------------------------
    log_file = os.path.join(dirs["logs"], "deseq2_canonical.log")
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        force=True
    )
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    console.setFormatter(formatter)
    logging.getLogger("").addHandler(console)

    logging.info("Starting Canonical DESeq2 Pipeline")
    
    # Reproducibility Seed
    np.random.seed(args.seed)

    # ---------------------------------------------------------
    # 2. Load and Audit Metadata
    # ---------------------------------------------------------
    meta_path = os.path.join(project_dir, "data/metadata/GSE117357_complete_mapping.csv")
    if not os.path.isfile(meta_path):
        logging.critical(f"Metadata file not found: {meta_path}")
        sys.exit(1)
        
    logging.info(f"Loading metadata from: {meta_path}")
    meta_df = pd.read_csv(meta_path)
    
    required_meta_cols = ["SRR_ID", "Tissue", "Genotype"]
    missing_meta_cols = [c for c in required_meta_cols if c not in meta_df.columns]
    if missing_meta_cols:
        logging.critical(f"Missing required metadata columns: {missing_meta_cols}")
        sys.exit(1)
        
    if meta_df[["SRR_ID", "Tissue", "Genotype"]].isna().any().any():
        logging.critical("Metadata contains NA values in required columns.")
        sys.exit(1)
        
    meta_df["SRR_ID"] = meta_df["SRR_ID"].astype(str).str.strip()
    meta_df["Tissue"] = meta_df["Tissue"].astype(str).str.strip()
    meta_df["Genotype"] = meta_df["Genotype"].astype(str).str.strip()
    
    for col in ["SRR_ID", "Tissue", "Genotype"]:
        if (meta_df[col].str.lower() == "nan").any():
            logging.critical(f"Metadata column '{col}' contains literal 'nan' string values.")
            sys.exit(1)
            
    meta_df = meta_df.set_index("SRR_ID")

    # Global genotype label validation
    observed_genotypes = set(meta_df["Genotype"].unique())
    expected_genotypes = {"WT", "KO"}
    if observed_genotypes != expected_genotypes:
        logging.critical(f"Genotype labels mismatch. Expected {expected_genotypes}, observed {observed_genotypes}")
        sys.exit(1)

    # ---------------------------------------------------------
    # 3. Load and Audit Counts
    # ---------------------------------------------------------
    counts_path = os.path.join(project_dir, "results/featurecounts_canonical/gene_counts.txt")
    if not os.path.isfile(counts_path):
        logging.critical(f"Counts file not found: {counts_path}")
        sys.exit(1)
        
    logging.info(f"Loading counts from: {counts_path}")
    counts_df = pd.read_csv(counts_path, sep="\t", comment="#", index_col=0)
    
    required_annot_cols = ["Chr", "Start", "End", "Strand", "Length"]
    missing_annot_cols = [c for c in required_annot_cols if c not in counts_df.columns]
    if missing_annot_cols:
        logging.critical(f"Missing featureCounts annotation columns: {missing_annot_cols}")
        sys.exit(1)
        
    gene_annotations = counts_df[required_annot_cols].copy()
    counts_df = counts_df.drop(columns=required_annot_cols)

    # Clean columns ('/path/to/SRR...bam' -> 'SRR...')
    counts_df.columns = [col.split("/")[-1].split(".")[0] for col in counts_df.columns]
    counts_df = counts_df.T

    # ---------------------------------------------------------
    # 4. Critical Validation Gates
    # ---------------------------------------------------------
    # Duplicates check
    if counts_df.index.duplicated().any():
        logging.critical("Duplicate sample IDs in counts matrix.")
        sys.exit(1)
    if counts_df.columns.duplicated().any():
        logging.critical("Duplicate gene IDs in counts matrix.")
        sys.exit(1)
    if meta_df.index.duplicated().any():
        logging.critical("Duplicate SRR_ID values in metadata.")
        sys.exit(1)

    # Numeric count validation (NA/inf safe)
    counts_df = counts_df.apply(pd.to_numeric, errors="coerce")
    values = counts_df.to_numpy(dtype=float)
    if not np.isfinite(values).all():
        logging.critical("Counts matrix contains NA, inf, or non-finite values.")
        sys.exit(1)
    if (values < 0).any():
        logging.critical("Counts matrix contains negative values.")
        sys.exit(1)
    if not np.all(np.equal(values, np.floor(values))):
        logging.critical("Counts matrix contains non-integer values.")
        sys.exit(1)
    counts_df = counts_df.astype(np.int64)

    # Sample Intersection (Fail Loudly)
    missing_in_meta = sorted(set(counts_df.index) - set(meta_df.index))
    missing_in_counts = sorted(set(meta_df.index) - set(counts_df.index))
    
    with open(os.path.join(dirs["qc"], "sample_matching_audit.csv"), "w") as f:
        f.write("Status,Sample_ID\n")
        for s in missing_in_meta: f.write(f"MissingInMeta,{s}\n")
        for s in missing_in_counts: f.write(f"MissingInCounts,{s}\n")
        for s in counts_df.index.intersection(meta_df.index): f.write(f"Matched,{s}\n")

    if missing_in_meta or missing_in_counts:
        logging.critical(f"Sample mismatch. Missing in meta: {missing_in_meta}. Missing in counts: {missing_in_counts}")
        sys.exit(1)

    # Tissue label validation
    observed_tissues = set(meta_df["Tissue"].dropna().unique())
    expected_tissues = {"Hippocampus", "Prefrontal-Cortex", "Striatum"}
    if observed_tissues != expected_tissues:
        logging.critical(f"Tissue labels mismatch. Expected {expected_tissues}, observed {observed_tissues}")
        sys.exit(1)

    # Confounding Audit
    with open(os.path.join(dirs["qc"], "metadata_balance_audit.txt"), "w") as f:
        f.write("=== Metadata Confounding Audit ===\n\n")
        for tissue in sorted(expected_tissues):
            t_meta = meta_df[meta_df["Tissue"] == tissue]
            f.write(f"=== {tissue} ===\n")
            for col in t_meta.columns:
                if col not in ["Tissue", "Genotype"]:
                    try:
                        f.write(f"--- Genotype vs {col} ---\n")
                        f.write(pd.crosstab(t_meta["Genotype"], t_meta[col]).to_string() + "\n\n")
                    except Exception as e:
                        f.write(f"  [Could not generate crosstab: {e}]\n\n")

    # ---------------------------------------------------------
    # 5. Tissue-Level Analysis
    # ---------------------------------------------------------
    summary_stats = []

    for tissue in sorted(expected_tissues):
        safe_tissue = tissue.lower().replace("-", "_")
        logging.info(f"--- Processing {tissue} ---")
        
        tissue_meta = meta_df[meta_df["Tissue"] == tissue].copy()
        tissue_counts = counts_df.loc[tissue_meta.index]
        
        # Explicit index alignment assertion
        if not tissue_counts.index.equals(tissue_meta.index):
            logging.critical(f"{tissue}: counts and metadata indices are misaligned.")
            sys.exit(1)

        # Enforce exactly expected N
        for genotype in ["WT", "KO"]:
            n = (tissue_meta["Genotype"] == genotype).sum()
            if n != args.expected_n_per_genotype:
                logging.critical(f"{tissue}: expected {args.expected_n_per_genotype} {genotype}, found {n}")
                sys.exit(1)

        # Encode categorical explicitly
        tissue_meta["Genotype"] = pd.Categorical(tissue_meta["Genotype"], categories=["WT", "KO"], ordered=True)

        logging.info("Fitting DESeq2 Model...")
        dds = DeseqDataSet(counts=tissue_counts, metadata=tissue_meta, design_factors="Genotype")
        dds.deseq2()
        
        logging.info("Calculating Stats (Contrast: KO vs WT)...")
        # Positive log2FC = Higher in KO
        stat_res = DeseqStats(dds, contrast=("Genotype", "KO", "WT"))
        import io, contextlib
        _buf = io.StringIO()
        with contextlib.redirect_stdout(_buf):
            stat_res.summary()  # Required: populates results_df
        logging.info(f"{tissue} summary:\n{_buf.getvalue().strip()}")
        res_df = stat_res.results_df.copy()
        
        # Validate gene annotation indices
        missing_gene_annot = res_df.index.difference(gene_annotations.index)
        if len(missing_gene_annot) > 0:
            logging.warning(f"{tissue}: {len(missing_gene_annot)} DESeq2 genes missing annotations.")
        
        # Merge annotation
        res_annotated = res_df.join(gene_annotations)
        
        # Save results
        out_csv = os.path.join(dirs["full_results"], f"{safe_tissue}_DESeq2_full.csv")
        out_csv_annot = os.path.join(dirs["full_results"], f"{safe_tissue}_DESeq2_annotated.csv")
        res_df.to_csv(out_csv)
        res_annotated.to_csv(out_csv_annot)
        
        # Low Count QC
        all_zero = (tissue_counts.sum(axis=0) == 0).sum()
        expressed_ge_10 = (tissue_counts.sum(axis=0) >= 10).sum()
        expressed_ge_10_in_3 = ((tissue_counts >= 10).sum(axis=0) >= 3).sum()
        
        # Stat QC
        n_input = len(res_df)
        n_valid_padj = res_df["padj"].notna().sum()
        n_sig_05 = (res_df["padj"] < 0.05).sum()
        n_sig_10 = (res_df["padj"] < 0.10).sum()
        lfc_thresh = args.lfc_threshold
        n_sig_05_lfc = ((res_df["padj"] < 0.05) & (res_df["log2FoldChange"].abs() > lfc_thresh)).sum()
        
        summary_stats.append({
            "Tissue": tissue,
            "n_samples": len(tissue_meta),
            "WT": int((tissue_meta["Genotype"]=="WT").sum()),
            "KO": int((tissue_meta["Genotype"]=="KO").sum()),
            "genes_all_zero": all_zero,
            "genes_ge_10": expressed_ge_10,
            "genes_ge_10_in_3": expressed_ge_10_in_3,
            "genes_input": n_input,
            "genes_valid_padj": n_valid_padj,
            "FDR<0.05": n_sig_05,
            "FDR<0.1": n_sig_10,
            f"FDR<0.05_and_absLFC>{lfc_thresh}": n_sig_05_lfc
        })

    # Save Summaries
    pd.DataFrame(summary_stats).to_csv(os.path.join(dirs["qc"], "deg_count_summary.csv"), index=False)

    # ---------------------------------------------------------
    # 6. Run Manifest
    # ---------------------------------------------------------
    try:
        import pydeseq2
        pydeseq2_version = pydeseq2.__version__
    except Exception:
        pydeseq2_version = "unknown"

    manifest = {
        "script": "58_deseq2_local_canonical.py",
        "command": " ".join(sys.argv),
        "counts_path": counts_path,
        "metadata_path": meta_path,
        "counts_sha256": sha256_file(counts_path),
        "metadata_sha256": sha256_file(meta_path),
        "counts_mtime": datetime.fromtimestamp(os.path.getmtime(counts_path)).isoformat(),
        "metadata_mtime": datetime.fromtimestamp(os.path.getmtime(meta_path)).isoformat(),
        "contrast": "Genotype KO vs WT",
        "interpretation": "Positive log2FC means higher expression in KO. Negative log2FC means lower expression in KO.",
        "python_version": sys.version.split()[0],
        "pandas_version": pd.__version__,
        "numpy_version": np.__version__,
        "pydeseq2_version": pydeseq2_version,
        "run_datetime": datetime.now().isoformat()
    }
    with open(os.path.join(base_out, "run_manifest.json"), "w") as f:
        json.dump(manifest, f, indent=4)

    logging.info("Canonical DESeq2 pipeline completed successfully.")

if __name__ == "__main__":
    main()
