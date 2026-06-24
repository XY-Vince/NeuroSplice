#!/bin/bash
# Generate rMATS BAM lists for GSE117357
# rMATS expects comma-separated BAM paths on a SINGLE line

ALIGNED_DIR="/Volumes/Untitled/NeuroSplice/data/aligned_star_reprocess"
OUT_DIR="/Volumes/Untitled/NeuroSplice/data/rmats_input"

mkdir -p "$OUT_DIR"

# Hippocampus samples
echo "Creating Hippocampus BAM lists..."
# KO: SRR7540620-SRR7540629
echo "${ALIGNED_DIR}/SRR7540620.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540621.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540622.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540623.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540624.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540625.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540626.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540627.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540628.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540629.Aligned.sortedByCoord.out.bam" > "${OUT_DIR}/ko_hippocampus_bams.txt"

# WT: SRR7540650-SRR7540659
echo "${ALIGNED_DIR}/SRR7540650.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540651.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540652.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540653.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540654.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540655.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540656.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540657.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540658.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540659.Aligned.sortedByCoord.out.bam" > "${OUT_DIR}/wt_hippocampus_bams.txt"

# Prefrontal Cortex samples
echo "Creating Prefrontal-Cortex BAM lists..."
# KO: SRR7540630-SRR7540639
echo "${ALIGNED_DIR}/SRR7540630.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540631.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540632.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540633.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540634.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540635.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540636.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540637.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540638.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540639.Aligned.sortedByCoord.out.bam" > "${OUT_DIR}/ko_prefrontal_bams.txt"

# WT: SRR7540660-SRR7540669
echo "${ALIGNED_DIR}/SRR7540660.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540661.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540662.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540663.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540664.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540665.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540666.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540667.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540668.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540669.Aligned.sortedByCoord.out.bam" > "${OUT_DIR}/wt_prefrontal_bams.txt"

# Striatum samples
echo "Creating Striatum BAM lists..."
# KO: SRR7540640-SRR7540649
echo "${ALIGNED_DIR}/SRR7540640.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540641.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540642.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540643.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540644.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540645.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540646.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540647.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540648.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540649.Aligned.sortedByCoord.out.bam" > "${OUT_DIR}/ko_striatum_bams.txt"

# WT: SRR7540670-SRR7540679
echo "${ALIGNED_DIR}/SRR7540670.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540671.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540672.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540673.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540674.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540675.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540676.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540677.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540678.Aligned.sortedByCoord.out.bam,${ALIGNED_DIR}/SRR7540679.Aligned.sortedByCoord.out.bam" > "${OUT_DIR}/wt_striatum_bams.txt"

echo "BAM lists created in $OUT_DIR"
echo "Verifying format (should be one comma-separated line each):"
for f in "$OUT_DIR"/*.txt; do
    echo "$(basename $f): $(wc -l < $f) lines"
done
