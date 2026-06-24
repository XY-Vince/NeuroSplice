#!/bin/bash
# =============================================================================
# Download mm10 Reference Genome and GENCODE Annotation
# =============================================================================

set -euo pipefail

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate combio

# Configuration
PROJECT_DIR="/Users/guozhenghui/Desktop/WXY/ComBio/GSE117357_Analysis"
REF_DIR="${PROJECT_DIR}/reference"
LOG_DIR="${PROJECT_DIR}/logs"

mkdir -p "${REF_DIR}" "${LOG_DIR}"

echo "[$(date)] Downloading mm10 reference files..."

# =============================================================================
# Step 1: Download GENCODE vM25 Annotation (GTF)
# =============================================================================
echo "[$(date)] Downloading GENCODE vM25 annotation..."

GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz"
GTF_FILE="${REF_DIR}/gencode.vM25.annotation.gtf"

if [[ ! -f "${GTF_FILE}" ]]; then
    curl -L -o "${GTF_FILE}.gz" "${GTF_URL}"
    gunzip "${GTF_FILE}.gz"
    echo "✓ GTF downloaded: ${GTF_FILE}"
else
    echo "[SKIP] GTF already exists"
fi

# =============================================================================
# Step 2: Download mm10 Genome FASTA (Primary Assembly)
# =============================================================================
echo "[$(date)] Downloading mm10 genome FASTA..."

# Using GENCODE primary assembly (matches GTF chromosome naming)
FASTA_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz"
FASTA_FILE="${REF_DIR}/mm10.fa"

if [[ ! -f "${FASTA_FILE}" ]]; then
    curl -L -o "${FASTA_FILE}.gz" "${FASTA_URL}"
    gunzip "${FASTA_FILE}.gz"
    echo "✓ FASTA downloaded: ${FASTA_FILE}"
else
    echo "[SKIP] FASTA already exists"
fi

# =============================================================================
# Step 3: Verify Chromosome Naming Consistency
# =============================================================================
echo "[$(date)] Verifying chromosome naming..."

FASTA_CHR=$(head -1 "${FASTA_FILE}")
GTF_CHR=$(head -100 "${GTF_FILE}" | grep -v "^#" | cut -f1 | head -1)

echo "FASTA first line: ${FASTA_CHR:0:30}"
echo "GTF first chromosome: ${GTF_CHR}"

if [[ "${FASTA_CHR}" == *"chr"* ]] && [[ "${GTF_CHR}" == "chr"* ]]; then
    echo "✓ Both use 'chr' prefix - compatible!"
elif [[ "${FASTA_CHR}" != *"chr"* ]] && [[ "${GTF_CHR}" != "chr"* ]]; then
    echo "✓ Both use numeric chromosomes - compatible!"
else
    echo "✗ WARNING: Chromosome naming mismatch detected!"
fi

# =============================================================================
# Step 4: Create Chromosome Sizes File
# =============================================================================
echo "[$(date)] Creating chromosome sizes file..."

if command -v samtools &> /dev/null; then
    samtools faidx "${FASTA_FILE}"
    cut -f1,2 "${FASTA_FILE}.fai" > "${REF_DIR}/mm10.chrom.sizes"
    echo "✓ Chrom sizes: ${REF_DIR}/mm10.chrom.sizes"
else
    echo "[SKIP] samtools not found - chrom sizes not created"
fi

# =============================================================================
# Step 5: Update Script Paths
# =============================================================================
echo "[$(date)] Updating script paths..."

# Update 03_preprocess.sh
sed -i.bak \
    -e "s|GENOME_DIR=\"/path/to/mm10_STAR_index\"|GENOME_DIR=\"${REF_DIR}/mm10_STAR_index\"|" \
    -e "s|GTF=\"/path/to/gencode.vM25.annotation.gtf\"|GTF=\"${REF_DIR}/gencode.vM25.annotation.gtf\"|" \
    -e "s|FASTA=\"/path/to/mm10.fa\"|FASTA=\"${REF_DIR}/mm10.fa\"|" \
    "${PROJECT_DIR}/scripts/03_preprocess.sh"

# Update 04_sanity_check.sh
sed -i.bak \
    -e "s|GTF=\"/path/to/gencode.vM25.annotation.gtf\"|GTF=\"${REF_DIR}/gencode.vM25.annotation.gtf\"|" \
    "${PROJECT_DIR}/scripts/04_sanity_check.sh"

echo "✓ Script paths updated"

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "============================================"
echo "Reference files downloaded to: ${REF_DIR}"
echo "============================================"
echo "FASTA: ${FASTA_FILE}"
echo "GTF:   ${GTF_FILE}"
echo ""
echo "NEXT STEPS:"
echo "1. Build STAR index: bash scripts/03_preprocess.sh --build-index"
echo "   (This takes ~30-60 minutes and ~30GB RAM)"
echo "2. Run data download: bash scripts/01_download_data.sh"
echo "============================================"
