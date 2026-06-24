#!/bin/bash
# 06_process_existing.sh - Align samples that have FASTQ but no BAM
set -euo pipefail

# Activate environment
source /opt/miniconda3/bin/activate combio

DATA_DIR="/Volumes/Untitled/NeuroSplice"
RAW_DIR="${DATA_DIR}/data/raw"
BAM_DIR="${DATA_DIR}/aligned/bam"

echo "[$(date)] Scanning for unprocessed FASTQs..."

for fastq in "${RAW_DIR}"/SRR*.fastq.gz; do
    if [[ ! -f "$fastq" ]]; then continue; fi
    
    # Extract SRR ID
    filename=$(basename "$fastq")
    SRR="${filename%%.*}"
    SRR="${SRR%_trimmed}" # Handle possible suffixes
    
    # Check if BAM exists
    if [[ -f "${BAM_DIR}/${SRR}_Aligned.sortedByCoord.out.bam" ]]; then
        echo "[SKIP] ${SRR} already aligned."
    else
        echo "[PROCESS] ${SRR} missing alignment. Triggering..."
        bash scripts/03_preprocess.sh --sample "${SRR}"
    fi
done
