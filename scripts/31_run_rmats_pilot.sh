#!/bin/bash
set -euo pipefail

# Activate rMATS env
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate combio-rmats

# Paths
PROJECT_DIR="/Volumes/Untitled/NeuroSplice"
STAR_BAM_DIR="${PROJECT_DIR}/aligned/star_pilot"
HISAT2_BAM_DIR="${PROJECT_DIR}/aligned/bam"
RESULTS_DIR="${PROJECT_DIR}/results/pilot_comparison"

mkdir -p "${RESULTS_DIR}/star" "${RESULTS_DIR}/hisat2"

# Define Samples
WT="SRR7540620"
KO="SRR7540650"

# 1. Run rMATS on STAR (Pilot)
echo "--- Running rMATS on STAR Pilot BAMs ---"
# Note: filenames in star_pilot are likely SRR..._Aligned.sortedByCoord.out.bam
echo "${STAR_BAM_DIR}/${WT}_Aligned.sortedByCoord.out.bam" > "${RESULTS_DIR}/b1_star.txt"
echo "${STAR_BAM_DIR}/${KO}_Aligned.sortedByCoord.out.bam" > "${RESULTS_DIR}/b2_star.txt"

rmats.py --b1 "${RESULTS_DIR}/b1_star.txt" \
         --b2 "${RESULTS_DIR}/b2_star.txt" \
         --gtf "${PROJECT_DIR}/reference/gencode.vM25.annotation.gtf" \
         -t single \
         --readLength 100 \
         --od "${RESULTS_DIR}/star" \
         --tmp "${RESULTS_DIR}/star/tmp" \
         --nthread 4

# 2. Run rMATS on HISAT2 (Subset)
echo "--- Running rMATS on HISAT2 Subset BAMs ---"
# Filenames: SRR..._Aligned.sortedByCoord.out.bam
echo "${HISAT2_BAM_DIR}/${WT}_Aligned.sortedByCoord.out.bam" > "${RESULTS_DIR}/b1_hisat.txt"
echo "${HISAT2_BAM_DIR}/${KO}_Aligned.sortedByCoord.out.bam" > "${RESULTS_DIR}/b2_hisat.txt"

rmats.py --b1 "${RESULTS_DIR}/b1_hisat.txt" \
         --b2 "${RESULTS_DIR}/b2_hisat.txt" \
         --gtf "${PROJECT_DIR}/reference/gencode.vM25.annotation.gtf" \
         -t single \
         --readLength 100 \
         --od "${RESULTS_DIR}/hisat2" \
         --tmp "${RESULTS_DIR}/hisat2/tmp" \
         --nthread 4

echo "Pilot Comparison Complete."
