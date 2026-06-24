#!/bin/bash
# 10_rmats.sh - Run rMATS for Alternative Splicing Analysis
set -euo pipefail

# Activate environment (fallback to 'combio-rmats' as 'combio' had conflicts)
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate combio-rmats

# Configuration
PROJECT_DIR="/Volumes/Untitled/NeuroSplice"
B1="${PROJECT_DIR}/data/metadata/wt_bams.txt"
B2="${PROJECT_DIR}/data/metadata/ko_bams.txt"
GTF="${PROJECT_DIR}/reference/gencode.vM25.annotation.gtf"
OUT_DIR="${PROJECT_DIR}/results/rmats"
TMP_DIR="${PROJECT_DIR}/results/rmats/tmp"
THREADS=16

mkdir -p "${OUT_DIR}" "${TMP_DIR}"

echo "[$(date)] Starting rMATS Analysis..."
echo "WT Group: $(awk -F, '{print NF}' $B1) samples"
echo "KO Group: $(awk -F, '{print NF}' $B2) samples"

# Run rMATS
# Note: --readLength 76 based on alignment logs
# Note: -t single based on alignment logs
rmats.py \
    --b1 "${B1}" \
    --b2 "${B2}" \
    --gtf "${GTF}" \
    -t single \
    --readLength 76 \
    --nthread "${THREADS}" \
    --od "${OUT_DIR}" \
    --tmp "${TMP_DIR}" \
    --variable-read-length

if [ $? -eq 0 ]; then
    echo "[$(date)] rMATS completed successfully."
    echo "Results in: ${OUT_DIR}"
else
    echo "[ERROR] rMATS failed."
    exit 1
fi
