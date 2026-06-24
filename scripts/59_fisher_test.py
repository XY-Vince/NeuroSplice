#!/usr/bin/env python3
"""
Phase 3: Human Digital Phenotype Intersection (Fisher's Exact Test)
===================================================================
Tests whether Adgrl3 KO mouse AS/DEG gene sets overlap significantly
with the 37 ADHD/psychiatric-associated genes from:

    Liu & Borsari et al., Cell 2025 (DOI: 10.1016/j.cell.2024.11.012)
    "Digital phenotyping from wearables using AI characterizes
     psychiatric disorders and identifies genetic associations"

Statistical Design:
- 6 independent Fisher's Exact Tests (3 regions × 2 molecular layers)
- Background: all expressed genes (baseMean >= 10 from DESeq2)
- Bonferroni correction: significance threshold p < 0.0083 (0.05/6)

Gene Mapping:
- Human→Mouse ortholog mapping by standard capitalization conversion
  (Human: NOVA1 → Mouse: Nova1)
- Ensembl→Symbol mapping extracted from rMATS output files
"""

import os

import pandas as pd
import numpy as np
import sys
from scipy.stats import fisher_exact

# ============================================================
# Configuration
# ============================================================
PROJECT = os.path.expanduser("~/Desktop/WXY/NeuroSplice")
DESEQ2_DIR = os.path.join(PROJECT, "results/deseq2_hpc")
RMATS_DIR = os.path.join(PROJECT, "results/rmats_gse117357_hpc")
OUT_DIR = os.path.join(PROJECT, "results/fisher_test")
os.makedirs(OUT_DIR, exist_ok=True)

REGIONS = {
    "Hippocampus": "hippocampus",
    "Prefrontal-Cortex": "prefrontal",
    "Striatum": "striatum",
}

EVENT_TYPES = ["SE", "A3SS", "A5SS", "MXE", "RI"]

N_TESTS = 6  # 3 regions × 2 layers
BONFERRONI_THRESHOLD = 0.05 / N_TESTS  # 0.0083

# ============================================================
# 37 Human Genes (Liu & Borsari, Cell 2025, Table 1)
# ============================================================
HUMAN_GENES_37 = [
    "TMIGD3", "ADORA3", "RAP1A", "CHI3L2",
    "IRF2", "CASP3", "PRIMPOL",
    "FGFR2",
    "NOVA1",
    "KCNH5", "RHOJ",
    "CLEC10A", "DLG4",
    "RHBDL3", "RHOT1", "C17orf75", "ZNF207", "PSMD11",
    "LRRC37B", "CDK5R1", "MYO1D",
    "AMY1C",
    "GCLC", "CILK1", "ELOVL5", "FBXO9", "GCM1",
    "MAD1L1", "ELFN1", "PSMG3", "MAFK",
    "MYH6", "CMTM5", "IL25", "BCL2L2",
    "WWOX",
]
# Note: BCL2L2-PABPN1 (readthrough) excluded — no clean mouse ortholog.
# Note: C17orf75 is human-specific nomenclature; mouse ortholog is Cgref1.

# Human → Mouse ortholog mapping
# Standard conversion: capitalize first letter only (e.g., NOVA1 → Nova1)
# With manual overrides for known exceptions
MANUAL_OVERRIDES = {
    "C17orf75": "Cgref1",   # Human C17orf75 = Mouse Cgref1
    "CHI3L2": "Chil1",      # Chitinase-3-like protein family divergence
}


def human_to_mouse(gene: str) -> str:
    """Convert human gene symbol to mouse ortholog symbol."""
    if gene in MANUAL_OVERRIDES:
        return MANUAL_OVERRIDES[gene]
    # Standard: NOVA1 → Nova1, DLG4 → Dlg4
    return gene[0].upper() + gene[1:].lower()


def build_ensembl_to_symbol_map() -> dict:
    """
    Build Ensembl ID → gene symbol mapping from rMATS output files.
    rMATS files contain both GeneID (Ensembl) and geneSymbol columns.
    """
    mapping = {}
    for region_rmats in REGIONS.values():
        for event in EVENT_TYPES:
            fpath = os.path.join(RMATS_DIR, region_rmats, f"{event}.MATS.JCEC.txt")
            if not os.path.exists(fpath):
                continue
            df = pd.read_csv(fpath, sep="\t", usecols=["GeneID", "geneSymbol"])
            df["GeneID"] = df["GeneID"].str.strip('"')
            df["symbol"] = df["geneSymbol"].str.strip('"')
            df["GeneID_base"] = df["GeneID"].str.split(".").str[0]
            mapping.update(dict(zip(df["GeneID"], df["symbol"])))
            mapping.update(dict(zip(df["GeneID_base"], df["symbol"])))
    return mapping


def get_as_genes(region_rmats: str) -> set:
    """Extract significant AS genes (FDR < 0.05) from rMATS for a region."""
    sig_genes = set()
    for event in EVENT_TYPES:
        fpath = os.path.join(RMATS_DIR, region_rmats, f"{event}.MATS.JCEC.txt")
        if not os.path.exists(fpath):
            continue
        df = pd.read_csv(fpath, sep="\t")
        sig = df[(df["FDR"] < 0.05) & (df["IncLevelDifference"].abs() > 0.1)]
        symbols = sig["geneSymbol"].str.strip('"').unique()
        sig_genes.update(symbols)
    return sig_genes


def get_deg_genes(region_name: str, ens_to_sym: dict) -> tuple:
    """
    Extract DEG foreground and background gene sets.
    Returns (sig_gene_symbols, background_gene_symbols).
    """
    full_path = os.path.join(DESEQ2_DIR, f"{region_name}_DESeq2_full.csv")
    if not os.path.exists(full_path):
        print(f"  WARNING: DESeq2 file not found, skipping: {full_path}")
        return set(), set()
    
    df = pd.read_csv(full_path, index_col=0)

    # Map Ensembl IDs to symbols
    df["symbol"] = df.index.map(lambda x: ens_to_sym.get(x, ens_to_sym.get(x.split(".")[0], None)))

    # Background: baseMean >= 10, with a valid symbol
    bg = df[(df["baseMean"] >= 10) & (df["symbol"].notna())]
    bg_genes = set(bg["symbol"].unique())

    # Foreground: padj < 0.05, |log2FC| > 0.5
    # NaN padj values are excluded by this comparison (intentional)
    sig = bg[(bg["padj"] < 0.05) & (bg["log2FoldChange"].abs() > 0.5)]
    sig_genes = set(sig["symbol"].unique())

    return sig_genes, bg_genes


def run_fisher(foreground: set, background: set, target_genes: set, label: str) -> dict:
    """
    Run a one-sided Fisher's Exact Test (greater).

    Contingency table:
                    In Target   Not In Target
    In Foreground      a             b
    Not In Foreground  c             d
    """
    # Only consider target genes that are in the background
    target_in_bg = target_genes & background

    a = len(foreground & target_in_bg)        # In both
    b = len(foreground - target_in_bg)          # In foreground only
    c = len(target_in_bg - foreground)          # In target only
    d = len(background - foreground - target_in_bg)  # In neither

    table = [[a, b], [c, d]]
    odds_ratio, p_value = fisher_exact(table, alternative="greater")

    overlapping = sorted(foreground & target_in_bg)

    return {
        "Test": label,
        "Foreground_Size": len(foreground),
        "Background_Size": len(background),
        "Target_In_Background": len(target_in_bg),
        "Overlap": a,
        "Odds_Ratio": round(odds_ratio, 3),
        "P_Value": p_value,
        "Significant": "YES" if p_value < BONFERRONI_THRESHOLD else "no",
        "Overlapping_Genes": ", ".join(overlapping) if overlapping else "—",
    }


def main():
    print("=" * 70)
    print("Phase 3: Human Digital Phenotype Intersection (Fisher's Exact Test)")
    print("=" * 70)
    print(f"Reference: Liu & Borsari, Cell 2025")
    print(f"Human genes: {len(HUMAN_GENES_37)}")
    print(f"Bonferroni threshold: p < {BONFERRONI_THRESHOLD:.4f}")

    if not os.path.exists(RMATS_DIR):
        sys.exit(f"ERROR: RMATS directory not found: {RMATS_DIR}")
    if not os.path.exists(DESEQ2_DIR):
        sys.exit(f"ERROR: DESeq2 directory not found: {DESEQ2_DIR}")

    # Convert human genes to mouse orthologs
    mouse_targets = set()
    mapping_log = []
    for hg in HUMAN_GENES_37:
        mg = human_to_mouse(hg)
        mouse_targets.add(mg)
        mapping_log.append({"Human": hg, "Mouse": mg})

    print(f"Mouse orthologs: {len(mouse_targets)}")

    # Save mapping
    pd.DataFrame(mapping_log).to_csv(
        os.path.join(OUT_DIR, "human_to_mouse_mapping.csv"), index=False
    )

    # Build Ensembl → Symbol map
    print("\nBuilding Ensembl → Symbol mapping from rMATS files...")
    ens_to_sym = build_ensembl_to_symbol_map()
    print(f"  Mapped {len(ens_to_sym)} Ensembl IDs to symbols")

    # Run 6 Fisher's Exact Tests
    results = []

    for region_name, region_rmats in REGIONS.items():
        print(f"\n--- {region_name} ---")

        # Get DEG foreground + background
        deg_fg, deg_bg = get_deg_genes(region_name, ens_to_sym)
        print(f"  DEG: {len(deg_fg)} sig genes / {len(deg_bg)} background")

        # Get AS foreground (use DEG background for consistency)
        as_fg = get_as_genes(region_rmats)
        print(f"  AS:  {len(as_fg)} sig genes")

        # Fisher test 1: DEG
        r1 = run_fisher(deg_fg, deg_bg, mouse_targets, f"{region_name}_DEG")
        results.append(r1)

        # Fisher test 2: AS
        r2 = run_fisher(as_fg, deg_bg, mouse_targets, f"{region_name}_AS")
        results.append(r2)

    # Format and save results
    results_df = pd.DataFrame(results)

    print(f"\n{'=' * 70}")
    print("RESULTS")
    print(f"{'=' * 70}")
    for _, row in results_df.iterrows():
        sig_marker = "🟢" if row["Significant"] == "YES" else "⚪"
        print(f"  {sig_marker} {row['Test']:30s}  "
              f"overlap={row['Overlap']}/{row['Target_In_Background']}  "
              f"OR={row['Odds_Ratio']:>7}  "
              f"p={row['P_Value']:.4e}  "
              f"{row['Significant']}")
        if row["Overlap"] > 0:
            print(f"     → Genes: {row['Overlapping_Genes']}")

    # Save
    out_file = os.path.join(OUT_DIR, "fisher_results.csv")
    results_df.to_csv(out_file, index=False)
    print(f"\nResults saved: {out_file}")

    # Interpretation
    n_sig = sum(results_df["Significant"] == "YES")
    print(f"\n{'=' * 70}")
    print("INTERPRETATION")
    print(f"{'=' * 70}")
    if n_sig > 0:
        print(f"✓ {n_sig}/6 tests reached Bonferroni significance (p < {BONFERRONI_THRESHOLD:.4f})")
        print("  → Evidence for overlap between Adgrl3 KO and human ADHD digital phenotype")
    else:
        print(f"✗ 0/6 tests reached Bonferroni significance (p < {BONFERRONI_THRESHOLD:.4f})")
        print("  → No significant overlap detected.")
        print("  → This is not unexpected: the 37-gene list was derived from human GWAS,")
        print("    while our data reflects downstream transcriptomic effects of a single")
        print("    gene knockout. The pathway-level convergence may be more meaningful")
        print("    than direct gene-level overlap.")


if __name__ == "__main__":
    main()
