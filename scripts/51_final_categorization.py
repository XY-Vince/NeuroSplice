#!/usr/bin/env python3
# ============================================================
# Canonical Candidate Categorization (Script 51)
# ============================================================
# Goal: Replaces the outdated categorization script to consume
# the rigorous QC outputs from script 59.
#
# Categories:
#   A1: Same event recurrent at FDR < 0.05 in all 3 tissues
#   A2: Same gene affected in all 3 tissues at distinct event loci
#   B:  Gene restricted to 1-2 tissues with high |ΔPSI|
#   C:  Directional heterogeneity across tissues
#   D:  Technical/low-annotation flag (e.g. Gm10419)
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

# ============================================================
# Constants & Defaults
# ============================================================
DEFAULT_DPSI_B1 = 0.30
DEFAULT_DPSI_B2 = 0.20
TECHNICAL_CONTROLS = ["Gm10419"]

REQUIRED_COLUMNS = [
    "geneSymbol", "region", "event_type", "ID", "GeneID", "chr", "strand",
    "FDR", "PValue", "IncLevelDifference", "sig_tier", "JC_significant",
    "JCEC_FDR_significant", "direction_consistent", "sashimi_eligible",
    "direction"
]

# ============================================================
# Utilities
# ============================================================

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def parse_args():
    parser = argparse.ArgumentParser(description="Canonical Candidate Categorization")
    parser.add_argument("--project-dir", required=True,
                        help="Path to NeuroSplice project root")
    parser.add_argument("--dpsi-b1", type=float, default=DEFAULT_DPSI_B1,
                        help="Min |ΔPSI| for Category B (1 tissue)")
    parser.add_argument("--dpsi-b2", type=float, default=DEFAULT_DPSI_B2,
                        help="Min |ΔPSI| for Category B (2 tissues)")
    return parser.parse_args()

def setup_directories(project_dir):
    qc_dir = os.path.join(project_dir, "results/rmats_gse117357_canonical/qc")
    out_dir = os.path.join(project_dir, "results/rmats_gse117357_canonical/final_categorization")
    log_dir = os.path.join(project_dir, "results/rmats_gse117357_canonical/logs")
    
    if not os.path.isdir(qc_dir):
        print(f"FATAL: QC directory not found at {qc_dir}", file=sys.stderr)
        sys.exit(1)
        
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)
    return qc_dir, out_dir, log_dir

def validate_columns(df, path):
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        logging.critical(f"Missing required columns in {path}: {missing}")
        sys.exit(1)

def construct_event_signature(r):
    """
    Constructs a coordinate-based signature independent of rMATS ID.
    SE: chr:strand:exonStart_0base-exonEnd:upstreamES-upstreamEE:downstreamES-downstreamEE
    RI: chr:strand:riExonStart_0base-riExonEnd:upstreamES-upstreamEE:downstreamES-downstreamEE
    A3SS/A5SS: chr:strand:longExonStart_0base-longExonEnd:shortES-shortEE
    MXE: chr:strand:1stExonStart_0base-1stExonEnd:2ndExonStart_0base-2ndExonEnd
    """
    etype = r.get("event_type")
    c = r.get("chr")
    s = r.get("strand")
    
    try:
        if etype == "SE":
            return f"{c}:{s}:{r['exonStart_0base']}-{r['exonEnd']}:{r['upstreamES']}-{r['upstreamEE']}:{r['downstreamES']}-{r['downstreamEE']}"
        elif etype == "RI":
            return f"{c}:{s}:{r['riExonStart_0base']}-{r['riExonEnd']}:{r['upstreamES']}-{r['upstreamEE']}:{r['downstreamES']}-{r['downstreamEE']}"
        elif etype in ("A3SS", "A5SS"):
            return f"{c}:{s}:{r['longExonStart_0base']}-{r['longExonEnd']}:{r['shortES']}-{r['shortEE']}"
        elif etype == "MXE":
            return f"{c}:{s}:{r['1stExonStart_0base']}-{r['1stExonEnd']}:{r['2ndExonStart_0base']}-{r['2ndExonEnd']}"
        else:
            return np.nan
    except KeyError:
        return np.nan

# ============================================================
# Core Logic
# ============================================================

def categorize_genes(df, dpsi_b1, dpsi_b2):
    """
    Assigns each gene to a category (A1, A2, B, C, D) and builds a gene summary.
    """
    # Only consider significant events for deciding categories A, B, C
    # Technical/supplemental logic happens globally
    sig_df = df[df["FDR_significant"] == True].copy()
    
    gene_summaries = []
    gene_to_cat = {}
    
    for gene, group in df.groupby("geneSymbol"):
        total_events = len(group)
        sig_group = sig_df[sig_df["geneSymbol"] == gene]
        
        tech_flag = gene in TECHNICAL_CONTROLS
        n_strict  = (group["strict_effect"] == True).sum()
        n_sashimi = (group["sashimi_eligible"] == True).sum()
        max_abs_dpsi = group["IncLevelDifference"].abs().max() if not group["IncLevelDifference"].isna().all() else 0.0
        
        # Directions among significant events
        sig_dirs = set(sig_group["direction"].dropna()) - {"unknown", "no_difference"}
        if len(sig_dirs) == 2:
            dir_pattern = "mixed"
        elif "WT_higher" in sig_dirs:
            dir_pattern = "consistent_WT_higher"
        elif "KO_higher" in sig_dirs:
            dir_pattern = "consistent_KO_higher"
        else:
            dir_pattern = "none"
            
        regions = set(sig_group["region"].unique())
        n_regions = len(regions)
        
        category = "None"
        
        if tech_flag or (n_strict == 0 and n_sashimi == 0):
            category = "D"
            event_pan_tissue = False
            gene_pan_tissue = False
        else:
            # Check A1/A2/C
            event_pan_tissue = False
            if not sig_group.empty:
                for sig_str, event_grp in sig_group.groupby("event_signature"):
                    if len(set(event_grp["region"])) == 3:
                        event_pan_tissue = True
                        break
                        
            gene_pan_tissue = (n_regions == 3)
            
            if dir_pattern == "mixed":
                category = "C"
            elif event_pan_tissue:
                category = "A1"
            elif gene_pan_tissue:
                category = "A2"
            else:
                # B check
                if n_regions == 1 and max_abs_dpsi >= dpsi_b1:
                    category = "B"
                elif n_regions == 2 and max_abs_dpsi >= dpsi_b2:
                    category = "B"
                else:
                    category = "B_moderate"
                    
        gene_to_cat[gene] = category
        
        gene_summaries.append({
            "geneSymbol": gene,
            "final_category": category,
            "event_pan_tissue": event_pan_tissue,
            "gene_pan_tissue": gene_pan_tissue,
            "n_total_events": total_events,
            "n_fdr_sig_events": len(sig_group),
            "n_strict_events": n_strict,
            "n_sashimi_eligible": n_sashimi,
            "max_abs_dpsi": max_abs_dpsi,
            "DirectionPattern": dir_pattern,
            "sig_regions": "|".join(sorted(regions)) if regions else "none"
        })
        
    return pd.DataFrame(gene_summaries), gene_to_cat


# ============================================================
# Main
# ============================================================

def main():
    args = parse_args()
    
    qc_dir, out_dir, log_dir = setup_directories(args.project_dir)
    
    log_file = os.path.join(log_dir, "final_categorization.log")
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        force=True
    )
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    console.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    logging.getLogger("").addHandler(console)
    
    logging.info("=== Starting Canonical Candidate Categorization ===")
    
    in_csv = os.path.join(qc_dir, "rmats_candidate_events_ranked.csv")
    if not os.path.isfile(in_csv):
        logging.critical(f"Input file not found: {in_csv}")
        sys.exit(1)
        
    df = pd.read_csv(in_csv)
    
    # We implicitly rely on fields created by script 59.
    # Script 59 produced "JC_significant" and "JCEC_FDR_significant" and "sig_tier".
    # Let's derive the unified boolean flags explicitly to be completely sure.
    df["FDR_significant"] = df["FDR"] < 0.05
    df["strict_effect"] = df["sig_tier"] == "Strict_effect"
    df["near_miss"] = df["sig_tier"] == "Near-miss"
    
    # If the user script 59 had these, great. Let's validate.
    validate_columns(df, in_csv)
    
    logging.info(f"Loaded {len(df)} events from {in_csv}")
    
    # 1. Add event signatures
    df["event_signature"] = df.apply(construct_event_signature, axis=1)
    missing_sigs = df["event_signature"].isna().sum()
    if missing_sigs > 0:
        logging.warning(f"Failed to construct coordinate signature for {missing_sigs} events.")
        
    # 2. Add technical flags
    df["technical_flag"] = df["geneSymbol"].isin(TECHNICAL_CONTROLS)
    
    # 3. Gene-level categorization
    gene_summary_df, gene_cat_map = categorize_genes(df, args.dpsi_b1, args.dpsi_b2)
    
    # 4. Map category back to event-level
    df["final_category"] = df["geneSymbol"].map(gene_cat_map)
    
    # ============================================================
    # Export Tables
    # ============================================================
    
    # 1. Full Events Table
    out_all = os.path.join(out_dir, "final_category_events.csv")
    df.to_csv(out_all, index=False)
    logging.info(f"Wrote full events table: {out_all}")
    
    # 2. Gene Summary
    out_genes = os.path.join(out_dir, "final_category_gene_summary.csv")
    gene_summary_df.to_csv(out_genes, index=False)
    logging.info(f"Wrote gene summary table: {out_genes}")
    
    # 3. Main Figure Candidates
    main_df = df[~df["technical_flag"] & (df["strict_effect"] | df["near_miss"]) & df["sashimi_eligible"]].copy()
    out_main = os.path.join(out_dir, "main_figure_candidates.csv")
    main_df.to_csv(out_main, index=False)
    logging.info(f"Wrote main figure candidates: {out_main} ({len(main_df)} events)")
    
    # 4. Supplemental Candidates
    supp_df = df[~df["technical_flag"] & df["FDR_significant"] & (~df["sashimi_eligible"] | ~df["strict_effect"])].copy()
    out_supp = os.path.join(out_dir, "supplemental_candidates.csv")
    supp_df.to_csv(out_supp, index=False)
    logging.info(f"Wrote supplemental candidates: {out_supp} ({len(supp_df)} events)")
    
    # 5. Technical / Flagged Candidates
    tech_df = df[df["technical_flag"]].copy()
    out_tech = os.path.join(out_dir, "technical_flagged_candidates.csv")
    tech_df.to_csv(out_tech, index=False)
    logging.info(f"Wrote technical flagged candidates: {out_tech} ({len(tech_df)} events)")
    
    # ============================================================
    # Manifest
    # ============================================================
    manifest = {
        "script": "51_final_categorization.py",
        "command": " ".join(sys.argv),
        "inputs": {
            "rmats_candidate_events_ranked": in_csv,
            "sha256": sha256_file(in_csv)
        },
        "thresholds": {
            "dpsi_b1": args.dpsi_b1,
            "dpsi_b2": args.dpsi_b2,
            "fdr_significant": "< 0.05",
            "strict_effect": "FDR < 0.05 and |ΔPSI| >= 0.05"
        },
        "technical_controls": TECHNICAL_CONTROLS,
        "python_version": sys.version.split()[0],
        "pandas_version": pd.__version__,
        "numpy_version": np.__version__,
        "run_datetime": datetime.now().isoformat()
    }
    manifest_path = os.path.join(out_dir, "final_categorization_manifest.json")
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=4)
    logging.info(f"Run manifest written: {manifest_path}")
    logging.info("=== Canonical Candidate Categorization Complete ===")

if __name__ == "__main__":
    main()
