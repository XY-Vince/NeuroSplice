#!/usr/bin/env python3
# ============================================================
# Canonical rMATS QC, Thresholded Summary & Candidate Extraction
# Script: 59_rmats_qc_and_candidates.py
# ============================================================
# Outputs:
#   qc/rmats_group_directionality_audit.txt
#   qc/rmats_thresholded_summary.csv
#   qc/rmats_candidate_events_all.csv
#   qc/rmats_candidate_events_ranked.csv
#   qc/sashimi_eligibility_candidates.csv
#   qc/rmats_qc_candidates_manifest.json
# ============================================================

import os
import sys
import json
import logging
import argparse
import hashlib
from datetime import datetime

import numpy as np
import pandas as pd

# ============================================================
# Constants
# ============================================================

CANDIDATE_GENES = [
    "Pts", "Neil2", "Tox3", "Unc13b", "Lrp8",
    "Myo9b", "Arid5a", "Bcl2l11", "Gm10419",
]

EVENT_TYPES = ["SE", "A3SS", "A5SS", "MXE", "RI"]

# rMATS output dir name → (wt_bamlist_filename, ko_bamlist_filename)
REGIONS = {
    "hippocampus": ("wt_hippocampus_bams.txt", "ko_hippocampus_bams.txt"),
    "prefrontal":  ("wt_prefrontal_bams.txt",  "ko_prefrontal_bams.txt"),
    "striatum":    ("wt_striatum_bams.txt",     "ko_striatum_bams.txt"),
}

DPSI_THRESHOLDS = [0.05, 0.10, 0.20]
FDR_STRICT   = 0.05
FDR_NEAR     = 0.10
DEFAULT_MIN_JUNCTION_READS = 10
DEFAULT_MIN_VALID_PSI      = 8
N_PER_GROUP                = 10


# ============================================================
# Utilities
# ============================================================

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def safe_nanmean(xs):
    arr = np.array(xs, dtype=float)
    return float(np.nanmean(arr)) if np.isfinite(arr).any() else np.nan


def parse_args():
    parser = argparse.ArgumentParser(description="Canonical rMATS QC and Candidate Extraction")
    parser.add_argument("--project-dir", required=True,
                        help="Path to NeuroSplice project root")
    parser.add_argument("--min-junction-reads", type=int, default=DEFAULT_MIN_JUNCTION_READS,
                        help="Min total junction reads per group for sashimi eligibility")
    parser.add_argument("--min-valid-psi", type=int, default=DEFAULT_MIN_VALID_PSI,
                        help="Min samples with valid (non-NA) PSI per group")
    return parser.parse_args()


def setup_directories(project_dir):
    rmats_dir = os.path.join(project_dir, "results/rmats_gse117357_canonical")
    if not os.path.isdir(rmats_dir):
        print(f"FATAL: rMATS output directory not found: {rmats_dir}", file=sys.stderr)
        sys.exit(1)
    qc_dir  = os.path.join(rmats_dir, "qc")
    log_dir = os.path.join(rmats_dir, "logs")
    os.makedirs(qc_dir,  exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)
    return rmats_dir, qc_dir, log_dir


def load_rmats_table(path, event_type):
    """Load one rMATS MATS.JC/JCEC file with full edge-case handling."""
    if not os.path.isfile(path):
        logging.warning(f"File not found, skipping: {path}")
        return None
    df = pd.read_csv(path, sep="\t")
    if "geneSymbol" not in df.columns:
        logging.critical(f"{path}: missing geneSymbol column")
        sys.exit(1)
    # rMATS produces a duplicate 'ID' column — pandas renames it 'ID.1'. Drop it.
    if "ID.1" in df.columns:
        df = df.drop(columns=["ID.1"])
    # Strip literal double-quotes wrapping GeneID / geneSymbol values
    for col in ("GeneID", "geneSymbol"):
        if col in df.columns:
            df[col] = df[col].astype(str).str.strip('"').str.strip()
    # Coerce numeric stats
    for col in ("PValue", "FDR", "IncLevelDifference"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    df["event_type"] = event_type
    return df


def parse_comma_floats(value, na_token="NA"):
    """Parse a comma-separated PSI string to a list of floats (np.nan for NA)."""
    if pd.isna(value):
        return [np.nan]
    parts = str(value).split(",")
    result = []
    for p in parts:
        p = p.strip()
        if p in (na_token, "nan", ""):
            result.append(np.nan)
        else:
            try:
                result.append(float(p))
            except ValueError:
                result.append(np.nan)
    return result


def parse_comma_ints(value):
    """Parse a comma-separated read-count string to a list of ints."""
    if pd.isna(value):
        return [0]
    parts = str(value).split(",")
    result = []
    for p in parts:
        p = p.strip()
        try:
            result.append(int(p))
        except ValueError:
            result.append(0)
    return result


def annotate_tier(fdr, pvalue, dpsi_abs):
    """Return 'Strict_effect', 'FDR_significant', 'Near-miss', or 'Background'."""
    if pd.isna(fdr) and pd.isna(pvalue):
        return "Background"
    fdr_v   = fdr   if not pd.isna(fdr)   else 1.0
    pval_v  = pvalue if not pd.isna(pvalue) else 1.0
    dpsi_v  = dpsi_abs if not pd.isna(dpsi_abs) else 0.0

    if fdr_v < FDR_STRICT and dpsi_v >= 0.05:
        return "Strict_effect"
    if fdr_v < FDR_STRICT:
        return "FDR_significant"
    if fdr_v < FDR_NEAR or (pval_v < 0.05 and dpsi_v >= 0.05):
        return "Near-miss"
    return "Background"


# ============================================================
# Step 1 — Directionality Audit
# ============================================================

def run_directionality_audit(project_dir, meta_df, qc_dir):
    bam_dir = os.path.join(project_dir, "data/rmats_input")
    lines = [
        "=== rMATS Group Directionality Audit ===\n",
        "\n",
        "rMATS convention:\n",
        "  IncLevelDifference = mean(IncLevel1) - mean(IncLevel2)\n",
        "  Sample1 = b1 file = WT\n",
        "  Sample2 = b2 file = KO\n",
        "  Positive ΔPSI = higher inclusion in WT\n",
        "  Negative ΔPSI = higher inclusion in KO\n",
        "  Event-specific biological interpretation depends on event type.\n",
        "\n",
    ]
    all_ok = True

    for region, (wt_fname, ko_fname) in REGIONS.items():
        wt_path = os.path.join(bam_dir, wt_fname)
        ko_path = os.path.join(bam_dir, ko_fname)

        for fpath in (wt_path, ko_path):
            if not os.path.isfile(fpath):
                logging.critical(f"BAM list file not found: {fpath}")
                sys.exit(1)

        # BAM lists are comma-separated on a single line
        with open(wt_path) as f:
            wt_bams = [b.strip() for b in f.read().strip().split(",") if b.strip()]
        with open(ko_path) as f:
            ko_bams = [b.strip() for b in f.read().strip().split(",") if b.strip()]

        wt_srrs = [os.path.basename(b).split(".")[0] for b in wt_bams]
        ko_srrs = [os.path.basename(b).split(".")[0] for b in ko_bams]

        lines.append(f"=== {region} ===\n")
        lines.append(f"  b1 (Sample1): {wt_fname}  [{len(wt_srrs)} BAMs]\n")
        lines.append(f"  b2 (Sample2): {ko_fname}  [{len(ko_srrs)} BAMs]\n")
        lines.append(f"  Sample1 SRRs: {', '.join(wt_srrs)}\n")
        lines.append(f"  Sample2 SRRs: {', '.join(ko_srrs)}\n")

        region_ok = True

        # Assert Sample1 = WT
        missing_wt = [s for s in wt_srrs if s not in meta_df.index]
        if missing_wt:
            logging.critical(f"{region}: b1 SRRs not in metadata: {missing_wt}")
            sys.exit(1)
        wt_genos = meta_df.loc[wt_srrs, "Genotype"].unique().tolist()
        if wt_genos != ["WT"]:
            lines.append(f"  ASSERTION FAILED: b1 expected WT, got {wt_genos}\n")
            region_ok = False

        # Assert Sample2 = KO
        missing_ko = [s for s in ko_srrs if s not in meta_df.index]
        if missing_ko:
            logging.critical(f"{region}: b2 SRRs not in metadata: {missing_ko}")
            sys.exit(1)
        ko_genos = meta_df.loc[ko_srrs, "Genotype"].unique().tolist()
        if ko_genos != ["KO"]:
            lines.append(f"  ASSERTION FAILED: b2 expected KO, got {ko_genos}\n")
            region_ok = False

        if region_ok:
            lines.append(f"  ASSERTION PASSED: b1=WT, b2=KO confirmed.\n\n")
            logging.info(f"{region}: directionality confirmed (b1=WT, b2=KO)")
        else:
            lines.append("\n")
            all_ok = False

    out_path = os.path.join(qc_dir, "rmats_group_directionality_audit.txt")
    with open(out_path, "w") as f:
        f.writelines(lines)
    logging.info(f"Directionality audit written: {out_path}")

    if not all_ok:
        logging.critical("Directionality assertion failed. Aborting.")
        sys.exit(1)


# ============================================================
# Step 2 — Thresholded Event Summary
# ============================================================

def build_thresholded_summary(rmats_dir, qc_dir):
    rows = []
    for region in REGIONS:
        for event_type in EVENT_TYPES:
            for model in ("JC", "JCEC"):
                fpath = os.path.join(rmats_dir, region, f"{event_type}.MATS.{model}.txt")
                df = load_rmats_table(fpath, event_type)
                if df is None:
                    continue

                total      = len(df)
                valid_fdr  = int(df["FDR"].notna().sum())
                valid_dpsi = int(df["IncLevelDifference"].notna().sum())
                valid_both = int((df["FDR"].notna() & df["IncLevelDifference"].notna()).sum())
                fdr        = df["FDR"]
                dpsi       = df["IncLevelDifference"].abs()
                sig_mask   = fdr < FDR_STRICT

                # Directionality within significant events
                incd = df["IncLevelDifference"]
                s1_higher = int((sig_mask & (incd > 0)).sum())
                s2_higher = int((sig_mask & (incd < 0)).sum())

                row = {
                    "region":             region,
                    "model":              model,
                    "event_type":         event_type,
                    "total_events":       total,
                    "valid_fdr_events":   valid_fdr,
                    "valid_dpsi_events":  valid_dpsi,
                    "valid_fdr_and_dpsi_events": valid_both,
                    "sig_fdr05":          int(sig_mask.sum()),
                    "sig_fraction_fdr05": round(sig_mask.sum() / total, 6) if total else 0,
                    "sample1_WT_higher":  s1_higher,
                    "sample2_KO_higher":  s2_higher,
                }
                for thresh in DPSI_THRESHOLDS:
                    key  = f"sig_fdr05_dpsi{int(thresh*100):02d}"
                    fkey = f"sig_fraction_dpsi{int(thresh*100):02d}"
                    n    = int((sig_mask & (dpsi >= thresh)).sum())
                    row[key]  = n
                    row[fkey] = round(n / total, 6) if total else 0

                rows.append(row)

    out = pd.DataFrame(rows)
    out_path = os.path.join(qc_dir, "rmats_thresholded_summary.csv")
    out.to_csv(out_path, index=False)
    logging.info(f"Thresholded summary written ({len(rows)} rows): {out_path}")
    return out


# ============================================================
# Step 3 — Candidate Extraction
# ============================================================

def extract_candidates(rmats_dir, qc_dir):
    jc_records   = []
    jcec_records = []

    for region in REGIONS:
        for event_type in EVENT_TYPES:
            jc_path   = os.path.join(rmats_dir, region, f"{event_type}.MATS.JC.txt")
            jcec_path = os.path.join(rmats_dir, region, f"{event_type}.MATS.JCEC.txt")
            jc_df   = load_rmats_table(jc_path,   event_type)
            jcec_df = load_rmats_table(jcec_path, event_type)

            for df, model, records in [
                (jc_df,   "JC",   jc_records),
                (jcec_df, "JCEC", jcec_records),
            ]:
                if df is None:
                    continue
                mask = df["geneSymbol"].isin(CANDIDATE_GENES)
                subset = df[mask].copy()
                if subset.empty:
                    continue
                subset["region"] = region
                subset["model"]  = model
                subset["sig_tier"] = subset.apply(
                    lambda r: annotate_tier(
                        r["FDR"], r["PValue"],
                        abs(r["IncLevelDifference"]) if not pd.isna(r["IncLevelDifference"]) else np.nan
                    ), axis=1
                )
                records.append(subset)

    if not jc_records:
        logging.warning("No candidate gene events found in any JC file.")
        return pd.DataFrame(), pd.DataFrame()

    jc_all   = pd.concat(jc_records,   ignore_index=True)
    jcec_all = pd.concat(jcec_records, ignore_index=True) if jcec_records else pd.DataFrame()

    # Save the raw all-candidate table
    jc_all.to_csv(os.path.join(qc_dir, "rmats_candidate_events_all.csv"), index=False)
    logging.info(
        f"Candidate events (JC): {len(jc_all)} rows, "
        f"{jc_all['geneSymbol'].nunique()} genes, "
        f"{jc_all['region'].nunique()} regions"
    )

    # Merge JCEC significance flag onto JC using (region, event_type, ID)
    if not jcec_all.empty:
        key = ["region", "event_type", "ID"]
        if jc_all.duplicated(key).any() or jcec_all.duplicated(key).any():
            logging.critical("Duplicate (region, event_type, ID) keys found before merge.")
            sys.exit(1)
            
        jcec_cols = jcec_all[
            ["region", "event_type", "ID", "FDR", "PValue", "IncLevelDifference", "sig_tier"]
        ].copy()
        jcec_cols.columns = [
            "region", "event_type", "ID", 
            "JCEC_FDR", "JCEC_PValue", "JCEC_IncLevelDifference", "JCEC_sig_tier"
        ]
        
        merged = jc_all.merge(jcec_cols, on=["region", "event_type", "ID"], how="left")
        merged["JCEC_FDR_significant"] = merged["JCEC_sig_tier"].isin(["Strict_effect", "FDR_significant"])
        merged["JCEC_strict_effect"] = merged["JCEC_sig_tier"] == "Strict_effect"
        merged["JC_significant"] = merged["sig_tier"].isin(["Strict_effect", "FDR_significant"])
        
        both_have_dpsi = (
            merged["IncLevelDifference"].notna() & 
            merged["JCEC_IncLevelDifference"].notna()
        )
        merged["direction_consistent"] = (
            merged["JC_significant"] &
            merged["JCEC_FDR_significant"] &
            both_have_dpsi &
            (np.sign(merged["IncLevelDifference"]) == np.sign(merged["JCEC_IncLevelDifference"]))
        )
    else:
        merged = jc_all.copy()
        merged["JC_significant"]          = merged["sig_tier"].isin(["Strict_effect", "FDR_significant"])
        merged["JCEC_FDR"]                = np.nan
        merged["JCEC_PValue"]             = np.nan
        merged["JCEC_IncLevelDifference"] = np.nan
        merged["JCEC_sig_tier"]           = np.nan
        merged["JCEC_FDR_significant"]    = False
        merged["JCEC_strict_effect"]      = False
        merged["direction_consistent"]    = np.nan

    # Add direction column
    merged["direction"] = np.where(
        merged["IncLevelDifference"].isna(), "unknown",
        np.where(merged["IncLevelDifference"] > 0, "WT_higher",
            np.where(merged["IncLevelDifference"] < 0, "KO_higher", "no_difference"))
    )

    # Rank: FDR ASC → |IncLevelDifference| DESC
    merged["_abs_dpsi"] = merged["IncLevelDifference"].abs()
    
    tier_map = {"Strict_effect": 0, "FDR_significant": 1, "Near-miss": 2, "Background": 3}
    merged["sig_tier_rank"] = merged["sig_tier"].map(tier_map)
    
    ranked = merged.sort_values(
        ["sig_tier_rank", "FDR", "_abs_dpsi", "geneSymbol", "region"],
        ascending=[True, True, False, True, True]
    ).drop(columns=["_abs_dpsi", "sig_tier_rank"]).reset_index(drop=True)

    # We will save the ranked table AFTER adding sashimi eligibility columns.
    return ranked, jcec_all


# ============================================================
# Step 4 — Sashimi Eligibility
# ============================================================

def compute_sashimi_eligibility(ranked_df, min_reads, min_valid_psi, qc_dir):
    if ranked_df.empty:
        logging.warning("No ranked candidates — skipping sashimi eligibility.")
        return pd.DataFrame()

    rows = []
    for _, r in ranked_df.iterrows():
        ijc1 = parse_comma_ints(r.get("IJC_SAMPLE_1"))
        sjc1 = parse_comma_ints(r.get("SJC_SAMPLE_1"))
        ijc2 = parse_comma_ints(r.get("IJC_SAMPLE_2"))
        sjc2 = parse_comma_ints(r.get("SJC_SAMPLE_2"))
        psi1 = parse_comma_floats(r.get("IncLevel1"))
        psi2 = parse_comma_floats(r.get("IncLevel2"))

        totals1 = [i + s for i, s in zip(ijc1, sjc1)]
        totals2 = [i + s for i, s in zip(ijc2, sjc2)]
        
        if len(totals1) != N_PER_GROUP or len(totals2) != N_PER_GROUP:
            logging.warning(f"{r.get('event_type')} {r.get('ID')}: sample count mismatch ({len(totals1)} vs {len(totals2)})")

        median_g1 = np.median(totals1) if totals1 else 0
        median_g2 = np.median(totals2) if totals2 else 0
        
        n_ok_1 = sum(1 for t in totals1 if t >= min_reads)
        n_ok_2 = sum(1 for t in totals2 if t >= min_reads)

        valid1 = sum(1 for x in psi1 if not np.isnan(x))
        valid2 = sum(1 for x in psi2 if not np.isnan(x))
        
        mean_psi1 = safe_nanmean(psi1)
        mean_psi2 = safe_nanmean(psi2)

        reads_ok = (median_g1 >= min_reads and median_g2 >= min_reads and 
                    n_ok_1 >= min_valid_psi and n_ok_2 >= min_valid_psi)
        psi_ok   = valid1 >= min_valid_psi and valid2 >= min_valid_psi
        eligible = reads_ok and psi_ok

        rows.append({
            "gene":              r.get("geneSymbol"),
            "region":            r.get("region"),
            "event_type":        r.get("event_type"),
            "event_id":          r.get("ID"),
            "GeneID":            r.get("GeneID"),
            "chr":               r.get("chr"),
            "strand":            r.get("strand"),
            "sig_tier":          r.get("sig_tier"),
            "FDR":               r.get("FDR"),
            "PValue":            r.get("PValue"),
            "IncLevelDifference": r.get("IncLevelDifference"),
            "direction":         r.get("direction"),
            "mean_WT_PSI":       round(mean_psi1, 4) if not np.isnan(mean_psi1) else np.nan,
            "mean_KO_PSI":       round(mean_psi2, 4) if not np.isnan(mean_psi2) else np.nan,
            "mean_IJC_WT":       round(float(np.mean(ijc1)), 2),
            "mean_SJC_WT":       round(float(np.mean(sjc1)), 2),
            "mean_IJC_KO":       round(float(np.mean(ijc2)), 2),
            "mean_SJC_KO":       round(float(np.mean(sjc2)), 2),
            "median_total_reads_WT": median_g1,
            "median_total_reads_KO": median_g2,
            "n_samples_reads_ok_WT": n_ok_1,
            "n_samples_reads_ok_KO": n_ok_2,
            "valid_psi_WT":      valid1,
            "valid_psi_KO":      valid2,
            "JC_significant":    r.get("JC_significant"),
            "JCEC_FDR_significant": r.get("JCEC_FDR_significant"),
            "JCEC_strict_effect": r.get("JCEC_strict_effect"),
            "direction_consistent": r.get("direction_consistent"),
            "reads_ok":          reads_ok,
            "psi_ok":            psi_ok,
            "sashimi_eligible":  eligible,
        })

    out = pd.DataFrame(rows)
    out_path = os.path.join(qc_dir, "sashimi_eligibility_candidates.csv")
    out.to_csv(out_path, index=False)
    n_eligible = int(out["sashimi_eligible"].sum())
    logging.info(
        f"Sashimi eligibility written: {n_eligible}/{len(out)} events eligible "
        f"(min_reads>={min_reads}, valid_psi>={min_valid_psi})"
    )

    # Add the key sashimi fields back into ranked_df for a single unified table
    ranked_df["mean_WT_PSI"] = out["mean_WT_PSI"]
    ranked_df["mean_KO_PSI"] = out["mean_KO_PSI"]
    ranked_df["median_total_reads_WT"] = out["median_total_reads_WT"]
    ranked_df["median_total_reads_KO"] = out["median_total_reads_KO"]
    ranked_df["valid_psi_WT"] = out["valid_psi_WT"]
    ranked_df["valid_psi_KO"] = out["valid_psi_KO"]
    ranked_df["sashimi_eligible"] = out["sashimi_eligible"]

    ranked_path = os.path.join(qc_dir, "rmats_candidate_events_ranked.csv")
    ranked_df.to_csv(ranked_path, index=False)
    logging.info(f"Ranked candidate table written (with sashimi fields): {len(ranked_df)} events")

    return out


# ============================================================
# Main
# ============================================================

def main():
    args = parse_args()
    project_dir = args.project_dir

    rmats_dir, qc_dir, log_dir = setup_directories(project_dir)

    log_file = os.path.join(log_dir, "rmats_qc_candidates.log")
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        force=True,
    )
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    console.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    logging.getLogger("").addHandler(console)

    logging.info("=== Starting Canonical rMATS QC & Candidate Extraction ===")

    # Load metadata for directionality verification
    meta_path = os.path.join(project_dir, "data/metadata/GSE117357_complete_mapping.csv")
    if not os.path.isfile(meta_path):
        logging.critical(f"Metadata not found: {meta_path}")
        sys.exit(1)
    meta_df = pd.read_csv(meta_path)
    meta_df["SRR_ID"]   = meta_df["SRR_ID"].astype(str).str.strip()
    meta_df["Genotype"] = meta_df["Genotype"].astype(str).str.strip()
    meta_df = meta_df.set_index("SRR_ID")

    # Step 1: Directionality audit
    run_directionality_audit(project_dir, meta_df, qc_dir)

    # Step 2: Thresholded event summary
    build_thresholded_summary(rmats_dir, qc_dir)

    # Step 3: Candidate extraction
    ranked_df, jcec_all = extract_candidates(rmats_dir, qc_dir)

    # Step 4: Sashimi eligibility
    compute_sashimi_eligibility(
        ranked_df,
        min_reads=args.min_junction_reads,
        min_valid_psi=args.min_valid_psi,
        qc_dir=qc_dir,
    )

    # Run manifest
    manifest = {
        "script":           "59_rmats_qc_and_candidates.py",
        "command":          " ".join(sys.argv),
        "project_dir":      project_dir,
        "rmats_dir":        rmats_dir,
        "metadata_path":    meta_path,
        "metadata_sha256":  sha256_file(meta_path),
        "candidate_genes":  CANDIDATE_GENES,
        "regions":          list(REGIONS.keys()),
        "event_types":      EVENT_TYPES,
        "fdr_strict":       FDR_STRICT,
        "fdr_near":         FDR_NEAR,
        "dpsi_thresholds":  DPSI_THRESHOLDS,
        "sashimi_min_junction_reads": args.min_junction_reads,
        "sashimi_min_valid_psi":      args.min_valid_psi,
        "directionality":   "IncLevelDifference = mean(IncLevel1[WT]) - mean(IncLevel2[KO]); Positive = higher inclusion in WT",
        "python_version":   sys.version.split()[0],
        "pandas_version":   pd.__version__,
        "numpy_version":    np.__version__,
        "run_datetime":     datetime.now().isoformat(),
    }
    manifest_path = os.path.join(qc_dir, "rmats_qc_candidates_manifest.json")
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=4)
    logging.info(f"Run manifest written: {manifest_path}")
    logging.info("=== Canonical rMATS QC & Candidate Extraction Complete ===")


if __name__ == "__main__":
    main()
