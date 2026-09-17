#!/usr/bin/env python3
"""
Generate correct rMATS BAM input lists from the verified metadata.

Source of truth: data/metadata/GSE117357_complete_mapping.csv
BAM location:   /Volumes/Untitled/Archive/Projects/NeuroSplice/aligned/bam/
Output:          data/rmats_input/  (6 files: WT/KO × 3 brain regions)

rMATS format: comma-separated BAM paths on a single line per group.
"""

import csv
import os
from collections import defaultdict

METADATA = "data/metadata/GSE117357_complete_mapping.csv"
BAM_ROOT = "/Volumes/Untitled/Archive/Projects/NeuroSplice/aligned/bam"
BAM_SUFFIX = "_Aligned.sortedByCoord.out.bam"
OUTPUT_DIR = "data/rmats_input"

# Tissue name normalization for filenames
TISSUE_MAP = {
    "Hippocampus": "hippocampus",
    "Prefrontal-Cortex": "prefrontal",
    "Striatum": "striatum",
}

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Read metadata
    groups = defaultdict(list)
    with open(METADATA) as f:
        reader = csv.DictReader(f)
        for row in reader:
            genotype = row["Genotype"].lower()  # wt or ko
            tissue = TISSUE_MAP[row["Tissue"]]
            srr = row["SRR_ID"]
            bam_path = os.path.join(BAM_ROOT, f"{srr}{BAM_SUFFIX}")
            groups[(genotype, tissue)].append((srr, bam_path))

    # Write output files and audit log
    print(f"{'File':<35} {'Samples':>7}  SRR range")
    print("-" * 75)

    for (genotype, tissue), samples in sorted(groups.items()):
        filename = f"{genotype}_{tissue}_bams.txt"
        filepath = os.path.join(OUTPUT_DIR, filename)

        # Sort by SRR ID for consistency
        samples.sort(key=lambda x: x[0])
        bam_paths = [s[1] for s in samples]

        # rMATS format: single line, comma-separated
        with open(filepath, "w") as f:
            f.write(",".join(bam_paths))

        srr_ids = [s[0] for s in samples]
        print(f"{filename:<35} {len(samples):>7}  {srr_ids[0]}–{srr_ids[-1]}")

    print(f"\n✅ 6 files written to {OUTPUT_DIR}/")
    print(f"   Source: {METADATA}")
    print(f"   BAM root: {BAM_ROOT}")

if __name__ == "__main__":
    main()
