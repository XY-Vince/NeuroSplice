#!/bin/bash
set -euo pipefail

# Activate Salmon Environment
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate salmon_env

# Directories
PROJECT_DIR="/Volumes/Untitled/NeuroSplice"
INDEX_DIR="${PROJECT_DIR}/reference/salmon_index_mm10"
TRANSCRIPT_FA="${PROJECT_DIR}/reference/gencode.vM25.transcripts.fa.gz"
FASTQ_DIR="${PROJECT_DIR}/data/raw" # Or trimmed? Using fastp cleaned data is better if available.
# Check where trimmed data is. Ideally "data/trimmed" but user script might have put it elsewhere.
# User mentioned "SRR7540620_trimmed.fastq.gz" in "data/trimmed/" in prev logs.
FASTQ_DIR="${PROJECT_DIR}/data/trimmed"
RESULTS_DIR="${PROJECT_DIR}/results/salmon_quant"

mkdir -p "${RESULTS_DIR}"

# 1. Build Index (if not exists)
if [ ! -d "${INDEX_DIR}" ]; then
    echo "Building Salmon Index..."
    salmon index -t "${TRANSCRIPT_FA}" -i "${INDEX_DIR}" -p 16
else
    echo "Salmon Index found at ${INDEX_DIR}"
fi

# 2. Quantify Samples
echo "Starting Quantification..."

# Get list of Samples (SRR)
# We can use metadata/sample_names.txt or just iterate FASTQs
# Let's iterate the manifest files we created
ALL_SAMPLES=$(cat ${PROJECT_DIR}/data/metadata/wt_*.txt ${PROJECT_DIR}/data/metadata/ko_*.txt | tr ',' '\n' | awk -F'/' '{print $NF}' | sed 's/_Aligned.*//' | sort | uniq)

# That gives SRRs. But we need FASTQ paths.
# Assuming single-end (based on fastp commands seen earlier/rMATS config "single")
# FASTQ should be data/trimmed/SRRxxxxxxx_trimmed.fastq.gz

for SRR in ${ALL_SAMPLES}; do
    echo "Processing ${SRR}..."
    FASTQ="${FASTQ_DIR}/${SRR}_trimmed.fastq.gz"
    OUT_DIR="${RESULTS_DIR}/${SRR}"
    
    if [ ! -f "${FASTQ}" ]; then
        echo "Warning: FASTQ not found for ${SRR} at ${FASTQ}"
        continue
    fi
    
    if [ -d "${OUT_DIR}" ]; then
        echo "Skipping ${SRR} (Already exists)"
        continue
    fi

    salmon quant -i "${INDEX_DIR}" -l A \
        -r "${FASTQ}" \
        -p 8 --validateMappings \
        -o "${OUT_DIR}"
        
    echo "Finished ${SRR}"
done

echo "Salmon Quantification Complete."
