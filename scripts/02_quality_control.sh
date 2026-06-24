#!/bin/bash
# =============================================================================
# GSE117357 Quality Control Script
# FastQC + MultiQC
# =============================================================================

set -euo pipefail

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate combio

# Configuration
PROJECT_DIR="/Users/guozhenghui/Desktop/WXY/ComBio/GSE117357_Analysis"
RAW_DIR="${PROJECT_DIR}/data/raw"
QC_DIR="${PROJECT_DIR}/results/qc"
LOG_DIR="${PROJECT_DIR}/logs"

# Threads
THREADS=8

# =============================================================================
# Step 1: Run FastQC
# =============================================================================
echo "[$(date)] Running FastQC on raw reads..."

mkdir -p "${QC_DIR}/fastqc_raw"

# Run FastQC on all FASTQ files
fastqc \
    --threads ${THREADS} \
    --outdir "${QC_DIR}/fastqc_raw" \
    "${RAW_DIR}"/*.fastq.gz \
    2>&1 | tee "${LOG_DIR}/fastqc.log"

echo "[$(date)] FastQC complete"

# =============================================================================
# Step 2: Aggregate with MultiQC
# =============================================================================
echo "[$(date)] Aggregating QC reports with MultiQC..."

multiqc \
    "${QC_DIR}/fastqc_raw" \
    --outdir "${QC_DIR}" \
    --filename "multiqc_raw_reads" \
    --title "GSE117357 Raw Read Quality" \
    --force \
    2>&1 | tee -a "${LOG_DIR}/multiqc.log"

echo "[$(date)] MultiQC report: ${QC_DIR}/multiqc_raw_reads.html"

# =============================================================================
# Step 3: Quick QC Summary
# =============================================================================
echo "[$(date)] Generating QC summary..."

# Extract key metrics from FastQC data
echo "Sample,Total_Sequences,Sequence_Length,GC_Content,Adapter_Content" > "${QC_DIR}/qc_summary.csv"

for zip in "${QC_DIR}/fastqc_raw"/*.zip; do
    sample=$(basename "${zip}" _fastqc.zip)
    
    # Extract fastqc_data.txt from zip
    unzip -p "${zip}" "*/fastqc_data.txt" 2>/dev/null | \
    awk -v sample="${sample}" '
        /^Total Sequences/ {total=$3}
        /^Sequence length/ {len=$3}
        /^%GC/ {gc=$2}
        END {print sample","total","len","gc",N/A"}
    ' >> "${QC_DIR}/qc_summary.csv"
done

echo "[$(date)] QC summary: ${QC_DIR}/qc_summary.csv"
echo "[$(date)] Quality control complete!"
