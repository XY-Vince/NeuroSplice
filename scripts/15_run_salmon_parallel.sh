#!/bin/bash
set -euo pipefail

# Activate Salmon Environment
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate salmon_env

# Directories
PROJECT_DIR="/Volumes/Untitled/NeuroSplice"
INDEX_DIR="${PROJECT_DIR}/reference/salmon_index_mm10"
FASTQ_DIR="${PROJECT_DIR}/data/trimmed"
RESULTS_DIR="${PROJECT_DIR}/results/salmon_quant"

mkdir -p "${RESULTS_DIR}"

if [ ! -d "${INDEX_DIR}" ]; then
    echo "Error: Salmon Index not found at ${INDEX_DIR}"
    exit 1
fi

echo "Starting Parallel Quantification..."

ALL_SAMPLES=$(cat ${PROJECT_DIR}/data/metadata/wt_*.txt ${PROJECT_DIR}/data/metadata/ko_*.txt | tr ',' '\n' | awk -F'/' '{print $NF}' | sed 's/_Aligned.*//' | sort | uniq)

# Parallel Config
MAX_JOBS=8
job_count=0

for SRR in ${ALL_SAMPLES}; do
    FASTQ="${FASTQ_DIR}/${SRR}_trimmed.fastq.gz"
    OUT_DIR="${RESULTS_DIR}/${SRR}"
    
    if [ ! -f "${FASTQ}" ]; then
        echo "Warning: FASTQ not found for ${SRR}"
        continue
    fi
    
    if [ -d "${OUT_DIR}" ]; then
        # Check if completed (logs exist and contain 'done'?)
        # For now, skip if dir exists to resume.
        # But if it was interrupted, dir exists but partial.
        # We can check for 'quant.sf'.
        if [ -f "${OUT_DIR}/quant.sf" ]; then
             echo "Skipping ${SRR} (Already exists)"
             continue
        fi
    fi

    echo "Processing ${SRR} [Job ${job_count}]..."
    (
        salmon quant -i "${INDEX_DIR}" -l A \
            -r "${FASTQ}" \
            -p 2 --validateMappings \
            -o "${OUT_DIR}" > "${RESULTS_DIR}/${SRR}.log" 2>&1
    ) &
    
    ((job_count++))
    if (( job_count >= MAX_JOBS )); then
        wait
        job_count=0
    fi
done
wait
echo "Salmon Quantification Complete."
