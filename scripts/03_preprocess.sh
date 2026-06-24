#!/bin/bash
# =============================================================================
# GSE117357 Preprocessing Script
# Adapter Trimming + HISAT2 Alignment
# =============================================================================

set -euo pipefail

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate combio

# Configuration
PROJECT_DIR="/Volumes/Untitled/NeuroSplice"
RAW_DIR="${PROJECT_DIR}/data/raw"
TRIM_DIR="${PROJECT_DIR}/data/trimmed"
ALIGN_DIR="${PROJECT_DIR}/aligned/bam"
QC_DIR="${PROJECT_DIR}/results/qc"
LOG_DIR="${PROJECT_DIR}/logs"
META_DIR="${PROJECT_DIR}/data/metadata"

# Reference paths
GENOME_DIR="${PROJECT_DIR}/reference/mm10_hisat2_index"
GTF="${PROJECT_DIR}/reference/gencode.vM25.annotation.gtf"
FASTA="${PROJECT_DIR}/reference/mm10.fa"
INDEX_PREFIX="${GENOME_DIR}/mm10"

# Parameters
THREADS=16

# =============================================================================
# Pre-flight Check: Chromosome Naming
# =============================================================================
check_chromosome_names() {
    echo "[$(date)] Checking chromosome naming consistency..."
    
    # Check FASTA
    FASTA_CHR=$(head -1 "${FASTA}" | grep -o "^>chr" || echo "")
    
    # Check GTF
    GTF_CHR=$(head -100 "${GTF}" | cut -f1 | grep "^chr" | head -1 || echo "")
    
    if [[ -n "${FASTA_CHR}" ]] && [[ -n "${GTF_CHR}" ]]; then
        echo "✓ Both FASTA and GTF use 'chr' prefix (UCSC style)"
    elif [[ -z "${FASTA_CHR}" ]] && [[ -z "${GTF_CHR}" ]]; then
        echo "✓ Both FASTA and GTF use numeric chromosomes (Ensembl style)"
    else
        echo "✗ WARNING: Chromosome naming mismatch!"
        echo "  FASTA: $(head -1 "${FASTA}" | cut -c1-20)"
        echo "  GTF: $(head -100 "${GTF}" | cut -f1 | sort -u | head -5 | tr '\n' ' ')"
        echo "  This WILL cause alignment failures. Fix before proceeding."
        exit 1
    fi
}

# =============================================================================
# Step 1: Adapter Trimming (fastp)
# =============================================================================
trim_reads() {
    local SRR=$1
    local R1="${RAW_DIR}/${SRR}_1.fastq.gz"
    local R2="${RAW_DIR}/${SRR}_2.fastq.gz"
    local R1_ALT="${RAW_DIR}/${SRR}.fastq.gz"  # Alternative naming for single-end
    local OUT1="${TRIM_DIR}/${SRR}_1_trimmed.fastq.gz"
    local OUT1_SE="${TRIM_DIR}/${SRR}_trimmed.fastq.gz"  # Single-end output
    local OUT2="${TRIM_DIR}/${SRR}_2_trimmed.fastq.gz"
    local LOG="${LOG_DIR}/${SRR}_trim.log"
    local JSON="${TRIM_DIR}/${SRR}_fastp.json"
    local HTML="${TRIM_DIR}/${SRR}_fastp.html"
    
    # Check if paired-end (R2 exists)
    if [[ -f "${R2}" ]]; then
        echo "[$(date)] Trimming ${SRR} (paired-end)..."
        fastp \
            --in1 "${R1}" \
            --in2 "${R2}" \
            --out1 "${OUT1}" \
            --out2 "${OUT2}" \
            --json "${JSON}" \
            --html "${HTML}" \
            --thread 4 \
            --qualified_quality_phred 20 \
            --length_required 36 \
            --detect_adapter_for_pe \
            2>&1 | tee "${LOG}"
    elif [[ -f "${R1}" ]]; then
        # Single-end with _1 suffix
        echo "[$(date)] Trimming ${SRR} (single-end with _1 suffix)..."
        fastp \
            --in1 "${R1}" \
            --out1 "${OUT1}" \
            --json "${JSON}" \
            --html "${HTML}" \
            --thread 4 \
            --qualified_quality_phred 20 \
            --length_required 36 \
            2>&1 | tee "${LOG}"
    elif [[ -f "${R1_ALT}" ]]; then
        # Single-end without suffix (SRR.fastq.gz)
        echo "[$(date)] Trimming ${SRR} (single-end, no suffix)..."
        fastp \
            --in1 "${R1_ALT}" \
            --out1 "${OUT1_SE}" \
            --json "${JSON}" \
            --html "${HTML}" \
            --thread 4 \
            --qualified_quality_phred 20 \
            --length_required 36 \
            2>&1 | tee "${LOG}"
    else
        echo "[ERROR] No input file found for ${SRR}"
        return 1
    fi
}

# =============================================================================
# Step 2: HISAT2 Alignment
# =============================================================================
align_reads() {
    local SRR=$1
    local R1_GZ="${TRIM_DIR}/${SRR}_1_trimmed.fastq.gz"
    local R2_GZ="${TRIM_DIR}/${SRR}_2_trimmed.fastq.gz"
    local R1_SE_GZ="${TRIM_DIR}/${SRR}_trimmed.fastq.gz"  # Single-end alternative
    local PREFIX="${ALIGN_DIR}/${SRR}"
    local LOG="${LOG_DIR}/${SRR}_align.log"
    
    echo "[$(date)] Aligning ${SRR} with HISAT2..."
    
    # Check if paired-end
    if [[ -f "${R2_GZ}" ]] && [[ -f "${R1_GZ}" ]]; then
        echo "[$(date)] Running HISAT2 (paired-end)..."
        hisat2 -p ${THREADS} \
               --dta \
               -x "${INDEX_PREFIX}" \
               -1 "${R1_GZ}" \
               -2 "${R2_GZ}" \
               --summary-file "${LOG}" \
               | samtools sort -@ 4 -o "${PREFIX}_Aligned.sortedByCoord.out.bam" -
        
    elif [[ -f "${R1_GZ}" ]]; then
        # Single-end with _1 suffix
        echo "[$(date)] Running HISAT2 (single-end)..."
        hisat2 -p ${THREADS} \
               --dta \
               -x "${INDEX_PREFIX}" \
               -U "${R1_GZ}" \
               --summary-file "${LOG}" \
               | samtools sort -@ 4 -o "${PREFIX}_Aligned.sortedByCoord.out.bam" -

    elif [[ -f "${R1_SE_GZ}" ]]; then
        # Single-end without _1 suffix
        echo "[$(date)] Running HISAT2 (single-end, no suffix)..."
        hisat2 -p ${THREADS} \
               --dta \
               -x "${INDEX_PREFIX}" \
               -U "${R1_SE_GZ}" \
               --summary-file "${LOG}" \
               | samtools sort -@ 4 -o "${PREFIX}_Aligned.sortedByCoord.out.bam" -

    else
        echo "[ERROR] No trimmed input file found for ${SRR}"
        return 1
    fi
    
    # Index BAM
    samtools index "${PREFIX}_Aligned.sortedByCoord.out.bam"
}

# =============================================================================
# Step 3: Build HISAT2 Index (if needed)
# =============================================================================
build_hisat2_index() {
    echo "[$(date)] Building HISAT2 index..."
    
    mkdir -p "${GENOME_DIR}"
    
    # Build index with splice sites and exons for better accuracy
    echo "Extracting splice sites and exons..."
    hisat2_extract_splice_sites.py "${GTF}" > "${GENOME_DIR}/splicesites.txt"
    hisat2_extract_exons.py "${GTF}" > "${GENOME_DIR}/exons.txt"

    echo "Running hisat2-build..."
    hisat2-build -p ${THREADS} \
                 --ss "${GENOME_DIR}/splicesites.txt" \
                 --exon "${GENOME_DIR}/exons.txt" \
                 "${FASTA}" \
                 "${INDEX_PREFIX}" \
                 2>&1 | tee "${LOG_DIR}/hisat2_index.log"
    
    echo "[$(date)] HISAT2 index complete"
}

# =============================================================================
# Main Execution
# =============================================================================
main() {
    # Check if index exists (check for .1.ht2 file)
    if [[ ! -f "${INDEX_PREFIX}.1.ht2" ]]; then
        echo "HISAT2 index not found. Please:"
        echo "1. Update GENOME_DIR, GTF, and FASTA paths in this script"
        echo "2. Run: bash $0 --build-index"
        exit 1
    fi
    
    # Check chromosome naming
    check_chromosome_names
    
    mkdir -p "${TRIM_DIR}" "${ALIGN_DIR}"
    
    # Get sample list
    if [[ -f "${META_DIR}/srr_accessions.txt" ]]; then
        SAMPLES=$(cat "${META_DIR}/srr_accessions.txt")
    elif [[ -f "${META_DIR}/correct_srr_accessions.txt" ]]; then
        SAMPLES=$(cat "${META_DIR}/correct_srr_accessions.txt")
    else
        # Fallback to finding files
        SAMPLES=$(ls -1 "${RAW_DIR}"/*.fastq.gz | sed 's/.*\///; s/\..*//' | sort -u)
    fi
    
    # Process each sample
    for SRR in ${SAMPLES}; do
        # Trim - check both naming conventions
        if [[ ! -f "${TRIM_DIR}/${SRR}_1_trimmed.fastq.gz" ]] && [[ ! -f "${TRIM_DIR}/${SRR}_trimmed.fastq.gz" ]]; then
            trim_reads "${SRR}"
        else
            echo "[SKIP] ${SRR} already trimmed"
        fi
        
        # Align
        if [[ ! -f "${ALIGN_DIR}/${SRR}_Aligned.sortedByCoord.out.bam" ]]; then
            align_reads "${SRR}"
        else
            echo "[SKIP] ${SRR} already aligned"
        fi
    done
    
    # Post-alignment QC with MultiQC
    echo "[$(date)] Running post-alignment MultiQC..."
    multiqc \
        "${ALIGN_DIR}" "${TRIM_DIR}" \
        --outdir "${QC_DIR}" \
        --filename "multiqc_alignment" \
        --title "GSE117357 Alignment Quality" \
        --force
    
    echo "[$(date)] Preprocessing complete!"
    echo "BAM files: ${ALIGN_DIR}"
    echo "QC report: ${QC_DIR}/multiqc_alignment.html"
}

# Handle arguments
if [[ "${1:-}" == "--build-index" ]]; then
    build_hisat2_index
elif [[ "${1:-}" == "--sample" ]]; then
    SRR="$2"
    if [[ -z "${SRR}" ]]; then
        echo "Usage: $0 --sample <SRR_ID>"
        exit 1
    fi
    # Setup directories
    mkdir -p "${TRIM_DIR}" "${ALIGN_DIR}" "${LOG_DIR}"
    
    # Process single sample (trim and align)
    # Check both naming conventions for skip logic
    if [[ ! -f "${TRIM_DIR}/${SRR}_1_trimmed.fastq.gz" ]] && [[ ! -f "${TRIM_DIR}/${SRR}_trimmed.fastq.gz" ]]; then
        trim_reads "${SRR}"
    fi
    
    if [[ ! -f "${ALIGN_DIR}/${SRR}_Aligned.sortedByCoord.out.bam" ]]; then
        align_reads "${SRR}"
    fi
else
    main
fi
