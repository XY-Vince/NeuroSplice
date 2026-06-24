#!/bin/bash
# ============================================================
# DEPRECATED — DO NOT USE FOR CANONICAL ANALYSIS
# ============================================================
# This script runs featureCounts locally against BAMs in
# data/aligned_star_reprocess/ using GENCODE vM25 annotation.
# These BAMs and annotation are NOT compatible with the canonical
# Ensembl 102-based rMATS and DESeq2 analyses.
#
# USE INSTEAD: scripts/hpc/56_featurecounts_hpc.sh
#   (targets HPC BAMs + Ensembl 102 GTF, SLURM-ready)
# ============================================================
# featureCounts — Generate gene-level count matrix
# ============================================================
# Dataset:  GSE117357 (Adgrl3 KO, 60 samples)
# Library:  SINGLE-END, reverse-stranded
# Aligner:  STAR v2.7.10b
# ============================================================

set -euo pipefail

# Activate conda
source /opt/miniconda3/bin/activate combio

# Paths
PROJECT="/Volumes/Untitled/NeuroSplice"
BAM_DIR="${PROJECT}/data/aligned_star_reprocess"
GTF="${PROJECT}/reference/gencode.vM25.annotation.gtf"
OUT_DIR="${PROJECT}/results/featurecounts"
THREADS=8

mkdir -p "${OUT_DIR}"

# Collect all 60 BAM files into a bash array (robust to spaces)
shopt -s nullglob
BAMS=("${BAM_DIR}"/*.bam)
N_BAMS=${#BAMS[@]}

if (( N_BAMS == 0 )); then
    echo "ERROR: No BAM files found in ${BAM_DIR}"
    exit 1
fi

echo "============================================================"
echo "featureCounts — Gene Count Matrix"
echo "============================================================"
echo "BAM directory:  ${BAM_DIR}"
echo "GTF annotation: ${GTF}"
echo "Output:         ${OUT_DIR}/gene_counts.txt"
echo "Samples:        ${N_BAMS}"
echo "Strand:         reverse (-s 2)"
echo "Library:        single-end (no -p flag)"
echo "Threads:        ${THREADS}"
echo "============================================================"
echo ""

# Turn off exit-on-error temporarily to manually catch featureCounts exit code
set +e

# Run featureCounts
# KEY FLAGS:
#   -s 2       : reverse-stranded (confirmed by SRA metadata audit)
#   -t exon    : count reads overlapping exons
#   -g gene_id : summarize at gene level
#   -T ${THREADS} : threads
#   NO -p      : data is SINGLE-END (not paired)
#   -a         : GTF annotation
#   -o         : output file
featureCounts \
    -a "${GTF}" \
    -o "${OUT_DIR}/gene_counts.txt" \
    -s 2 \
    -t exon \
    -g gene_id \
    -T "${THREADS}" \
    --verbose \
    "${BAMS[@]}" \
    2>&1 | tee "${OUT_DIR}/featurecounts.log"

STATUS=${PIPESTATUS[0]} # Catch exit code of featureCounts specifically, not tee

# Re-enable exit-on-error
set -e

echo ""
echo "============================================================"
if [ ${STATUS} -eq 0 ]; then
    echo "✓ featureCounts completed successfully"
    echo ""
    # Quick stats
    NGENES=$(tail -n +3 "${OUT_DIR}/gene_counts.txt" | wc -l | tr -d ' ')
    echo "  Genes counted: ${NGENES}"
    echo "  Output file:   ${OUT_DIR}/gene_counts.txt"
    echo "  Summary file:  ${OUT_DIR}/gene_counts.txt.summary"
    echo ""
    echo "  File size:"
    ls -lh "${OUT_DIR}/gene_counts.txt"
else
    echo "✗ featureCounts failed with exit code ${STATUS}"
fi
echo "============================================================"
echo "Completed: $(date)"
