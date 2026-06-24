#!/bin/bash
set -euo pipefail

# Activate environment
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate combio-rmats

# Project Directories
PROJECT_DIR="/Volumes/Untitled/NeuroSplice"
METADATA_DIR="${PROJECT_DIR}/data/metadata"
RESULTS_DIR="${PROJECT_DIR}/results"
BAM_DIR="${PROJECT_DIR}/aligned/bam"
GTF="${PROJECT_DIR}/reference/gencode.vM25.annotation.gtf"

# Define Regions
REGIONS=("Hippocampus" "Prefrontal-Cortex" "Striatum")

echo "--- Generating Input Lists ---"

# Create temporary file with SRR and SampleName
paste "${METADATA_DIR}/srr_accessions.txt" "${METADATA_DIR}/sample_names.txt" > "${METADATA_DIR}/sample_map.txt"

for REGION in "${REGIONS[@]}"; do
    echo "Processing ${REGION}..."
    
    # Define Output Files (Lowercase for filenames)
    R_LOWER=$(echo "$REGION" | tr '[:upper:]' '[:lower:]' | sed 's/-/_/g')
    WT_FILE="${METADATA_DIR}/wt_${R_LOWER}.txt"
    KO_FILE="${METADATA_DIR}/ko_${R_LOWER}.txt"
    
    # Extract SRRs for WT and KO for this region
    # logic: grep region -> grep WT/KO -> get SRR -> format as BAM path
    
    grep "${REGION}" "${METADATA_DIR}/sample_map.txt" | grep "WT_" | awk -v bamParams="${BAM_DIR}" '{print bamParams "/" $1 "_Aligned.sortedByCoord.out.bam"}' | tr '\n' ',' | sed 's/,$//' > "${WT_FILE}"
    
    grep "${REGION}" "${METADATA_DIR}/sample_map.txt" | grep "KO_" | awk -v bamParams="${BAM_DIR}" '{print bamParams "/" $1 "_Aligned.sortedByCoord.out.bam"}' | tr '\n' ',' | sed 's/,$//' > "${KO_FILE}"
    
    # Validation
    WT_COUNT=$(awk -F, '{print NF}' "${WT_FILE}")
    KO_COUNT=$(awk -F, '{print NF}' "${KO_FILE}")
    
    echo "  ${REGION}: WT=${WT_COUNT}, KO=${KO_COUNT} samples."
    
    if [[ "$WT_COUNT" -ne 10 || "$KO_COUNT" -ne 10 ]]; then
        echo "Error: Expected 10 samples per group for ${REGION}, found WT=${WT_COUNT} KO=${KO_COUNT}."
        exit 1
    fi
done

echo "--- Launching rMATS Analyses ---"

for REGION in "${REGIONS[@]}"; do
    R_LOWER=$(echo "$REGION" | tr '[:upper:]' '[:lower:]' | sed 's/-/_/g')
    OUT_DIR="${RESULTS_DIR}/rmats_${R_LOWER}"
    WT_FILE="${METADATA_DIR}/wt_${R_LOWER}.txt"
    KO_FILE="${METADATA_DIR}/ko_${R_LOWER}.txt"
    
    echo "Starting rMATS for ${REGION}..."
    mkdir -p "${OUT_DIR}"
    
    # Run rMATS (Detailed logging per region)
    rmats.py \
        --b1 "${WT_FILE}" \
        --b2 "${KO_FILE}" \
        --gtf "${GTF}" \
        -t single \
        --readLength 76 \
        --nthread 16 \
        --od "${OUT_DIR}" \
        --tmp "${OUT_DIR}/tmp" \
        --variable-read-length \
        > "${RESULTS_DIR}/rmats_${R_LOWER}.log" 2>&1
        
    echo "  Finished ${REGION}. Results in ${OUT_DIR}"
done

echo "--- All Regions Processed ---"
