#!/bin/bash
# =============================================================================
# GSE117357 Data Acquisition Script (No entrez-direct)
# Uses direct SRA Run Selector export
# =============================================================================

set -euo pipefail

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate combio

# Configuration
PROJECT_DIR="/Volumes/Untitled/NeuroSplice"
DATA_DIR="${PROJECT_DIR}/data"
RAW_DIR="${DATA_DIR}/raw"
META_DIR="${DATA_DIR}/metadata"
LOG_DIR="${PROJECT_DIR}/logs"

# BioProject/SRA identifiers
BIOPROJECT="PRJNA481875"

# Number of parallel downloads
PARALLEL_DOWNLOADS=2

mkdir -p "${RAW_DIR}" "${META_DIR}" "${LOG_DIR}"

# =============================================================================
# Step 1: Fetch SRA Run Info via API
# =============================================================================
echo "[$(date)] Fetching run info from SRA..."

# Use SRA Run Selector API endpoint
curl -s "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_new?acc=${BIOPROJECT}&format=csv" \
    > "${META_DIR}/SRA_runinfo.csv"

# Check if we got data
if [[ ! -s "${META_DIR}/SRA_runinfo.csv" ]]; then
    echo "[ERROR] Failed to fetch run info. Trying alternative method..."
    
    # Alternative: use ENA API
    curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${BIOPROJECT}&result=read_run&fields=run_accession,sample_title,library_layout,fastq_ftp&format=tsv" \
        > "${META_DIR}/ENA_runinfo.tsv"
fi

echo "[$(date)] Run info saved"

# =============================================================================
# Step 2: Extract SRR Accessions
# =============================================================================
echo "[$(date)] Using existing sample accessions..."

# SKIPPING AUTO-GENERATION TO PREVENT CORRUPTION
# Parse the CSV/TSV to get SRR accessions
# if [[ -f "${META_DIR}/SRA_runinfo.csv" ]] && [[ -s "${META_DIR}/SRA_runinfo.csv" ]]; then
#     # SRA format: Run column
#     tail -n +2 "${META_DIR}/SRA_runinfo.csv" | cut -d',' -f1 | grep "^SRR" > "${META_DIR}/srr_accessions.txt" || true
# fi

# if [[ ! -s "${META_DIR}/srr_accessions.txt" ]] && [[ -f "${META_DIR}/ENA_runinfo.tsv" ]]; then
#     # ENA format: run_accession column
#     tail -n +2 "${META_DIR}/ENA_runinfo.tsv" | cut -f1 > "${META_DIR}/srr_accessions.txt"
# fi

# # If still empty, use known GSE117357 accessions
# if [[ ! -s "${META_DIR}/srr_accessions.txt" ]]; then
#     echo "[INFO] Using known accessions for GSE117357..."
#     # GSE117357 has 60 samples: SRR7642247 to SRR7642306
#     for i in $(seq 7642247 7642306); do
#         echo "SRR${i}"
#     done > "${META_DIR}/srr_accessions.txt"
# fi

N_SAMPLES=$(wc -l < "${META_DIR}/srr_accessions.txt" | tr -d ' ')
echo "[$(date)] Found ${N_SAMPLES} samples"
head -5 "${META_DIR}/srr_accessions.txt"

# =============================================================================
# Step 3: Check for Required Tools
# =============================================================================
check_tools() {
    local missing=0
    
    for tool in prefetch fasterq-dump; do
        if ! command -v $tool &> /dev/null; then
            echo "[ERROR] $tool not found. Install SRA Toolkit:"
            echo "  conda install -c bioconda sra-tools"
            echo "  OR: brew install sratoolkit"
            missing=1
        fi
    done
    
    if [[ $missing -eq 1 ]]; then
        echo ""
        echo "After installing, configure vdb-config:"
        echo "  vdb-config --interactive"
        exit 1
    fi
}

check_tools

# =============================================================================
# Step 4: Download FASTQ Files
# =============================================================================
echo "[$(date)] Starting FASTQ downloads..."
echo "[$(date)] This will download ${N_SAMPLES} samples..."

cd "${RAW_DIR}"

download_sample() {
    local SRR=$1
    local LOG="${LOG_DIR}/${SRR}_download.log"
    
    # Skip if already downloaded
    if ls "${SRR}"*.fastq.gz 1>/dev/null 2>&1; then
        echo "[SKIP] ${SRR} already downloaded"
        return 0
    fi
    
    echo "[$(date)] Downloading ${SRR}..." | tee -a "${LOG}"
    
    # Prefetch SRA file (more reliable than direct fastq-dump)
    prefetch "${SRR}" -O . --max-size 50G 2>> "${LOG}" || {
        echo "[ERROR] prefetch failed for ${SRR}" >> "${LOG}"
        return 1
    }
    
    # HANDLED BY 05_sra_stream.sh NOW
    # Convert to FASTQ
    # fasterq-dump "${SRR}" \
    #     --split-files \
    #     --threads 4 \
    #     --progress \
    #     2>> "${LOG}" || {
    #     echo "[ERROR] fasterq-dump failed for ${SRR}" >> "${LOG}"
    #     return 1
    # }
    
    # Compress with pigz (fast) or gzip (fallback)
    # if command -v pigz &> /dev/null; then
    #     pigz -p 4 "${SRR}"*.fastq 2>> "${LOG}"
    # else
    #     gzip "${SRR}"*.fastq 2>> "${LOG}"
    # fi
    
    # Cleanup SRA cache (Let sra_stream.sh handle this)
    # rm -rf "${SRR}"
    
    echo "[$(date)] Completed ${SRR}" | tee -a "${LOG}"
}

export -f download_sample
export LOG_DIR RAW_DIR

# Download samples sequentially for stability (parallel can overwhelm connections)
while read SRR; do
    download_sample "${SRR}"
done < "${META_DIR}/srr_accessions.txt"

echo "[$(date)] All downloads complete!"

# =============================================================================
# Step 5: Verify Downloads
# =============================================================================
echo "[$(date)] Verifying downloads..."

EXPECTED="${N_SAMPLES}"
DOWNLOADED=$(ls -1 "${RAW_DIR}"/*.fastq.gz 2>/dev/null | wc -l | tr -d ' ')

echo "Expected samples: ${EXPECTED}"
echo "Downloaded files: ${DOWNLOADED}"

# List what we have
ls -lh "${RAW_DIR}"/*.fastq.gz 2>/dev/null | head -10

echo ""
echo "[$(date)] Data acquisition complete!"
echo "Raw data location: ${RAW_DIR}"
echo "Metadata location: ${META_DIR}"
