#!/bin/bash
# ============================================================
# DEPRECATED — DO NOT USE FOR CANONICAL ANALYSIS
# ============================================================
# This script runs STAR alignment locally using GENCODE vM25 annotation.
# The canonical, validated results were generated on the Yale HPC
# (Bouchet/Gibbs Cluster) using Ensembl 102 GTF.
#
# Local BAMs in data/aligned_star_reprocess/ are NOT compatible with
# the Ensembl 102-based rMATS outputs (coordinate/gene ID mismatch).
#
# For any re-runs, use the HPC pipeline with:
#   GTF: /gpfs/gibbs/project/neurosplice_gse117357/references/Mus_musculus.GRCm38.102.gtf
#   BAMs output to: /gpfs/gibbs/project/neurosplice_gse117357/alignments/
# ============================================================
# Script 43: Re-process GSE117357 (60 Samples) using STAR + rMATS
# Objective: Standardize analysis pipeline with the new dataset (GSE173926).
# Input: Single-End Reads (76bp).
# Aligner: STAR (v2.7+).
# Splicing: rMATS (v4.3+).

PROJECT_DIR="/Volumes/Untitled/NeuroSplice"
RAW_DIR="$PROJECT_DIR/data/raw"
FASTQ_DIR="$PROJECT_DIR/data/trimmed"
ALIGNED_DIR="$PROJECT_DIR/data/aligned_star_reprocess"
RMATS_BASE="$PROJECT_DIR/results/rmats_star_reprocess"
STAR_INDEX="$PROJECT_DIR/reference/mm10_STAR_index"
GTF="$PROJECT_DIR/reference/gencode.vM25.annotation.gtf"
MAP_FILE="$PROJECT_DIR/data/metadata/sample_map.txt"
THREADS=8

# Use x86_64 STAR binary from GitHub (ARM64 STAR has spawn bug on macOS)
STAR_BIN="$PROJECT_DIR/tools/STAR_x86_64"
# Verify STAR binary exists
if [[ ! -f "$STAR_BIN" ]]; then
    echo "ERROR: x86_64 STAR binary not found at $STAR_BIN"
    echo "Download from: https://github.com/alexdobin/STAR/releases/download/2.7.11b/STAR_2.7.11b.zip"
    exit 1
fi

# Activate Conda Env (for samtools, rmats)
source /opt/miniconda3/bin/activate combio-rmats
export PATH=$PATH:/opt/miniconda3/envs/combio-rmats/bin:/opt/miniconda3/envs/combio/bin

mkdir -p "$ALIGNED_DIR"
mkdir -p "$RMATS_BASE"

echo "--- Step 1: Alignment (STAR) ---"
while read -r line; do
    SRR=$(echo "$line" | awk '{print $1}')
    SAMPLE_NAME=$(echo "$line" | awk '{print $2}')
    # Parse Condition and Tissue from "WT_Prefrontal-Cortex_Rep1"
    CONDITION=$(echo "$SAMPLE_NAME" | cut -d_ -f1)
    TISSUE=$(echo "$SAMPLE_NAME" | cut -d_ -f2)
    
    echo "Processing $SRR ($SAMPLE_NAME)..."
    # Input File
    R1="$FASTQ_DIR/${SRR}_trimmed.fastq.gz"
    OUT_PREFIX="$ALIGNED_DIR/${SRR}."
    
    if [[ ! -f "$R1" ]]; then
        echo "  [WARNING] Raw file $R1 not found. Skipping."
        continue
    fi
    
    # Check if BAM exists and is valid (size > 1MB)
    if [[ -f "${OUT_PREFIX}Aligned.sortedByCoord.out.bam" ]]; then
        SIZE=$(stat -f%z "${OUT_PREFIX}Aligned.sortedByCoord.out.bam")
        if (( SIZE > 1000000 )); then
            echo "  Already aligned and valid size ($SIZE bytes). Skipping."
            continue
        else
            echo "  BAM exists but too small ($SIZE bytes). Re-running."
        fi
    fi


    # STAR Alignment (x86_64 via Rosetta, gzcat for decompression)
    echo "  Starting STAR (x86_64 via Rosetta)..."
    
    # Ensure local tmp dir is clean
    rm -rf "/tmp/${SRR}_STARtmp"
    
    # Note: Single-end mode, so alignMatesGapMax is omitted
    arch -x86_64 "$STAR_BIN" --genomeDir "$STAR_INDEX" \
         --readFilesIn "$R1" \
         --readFilesCommand gzcat \
         --outFileNamePrefix "$OUT_PREFIX" \
         --outSAMtype BAM SortedByCoordinate \
         --outTmpDir "/tmp/${SRR}_STARtmp" \
         --runThreadN "$THREADS" \
         --outSAMstrandField intronMotif \
         --outFilterMultimapNmax 10 \
         --alignSJoverhangMin 8 \
         --alignSJDBoverhangMin 1 \
         --outFilterMismatchNmax 999 \
         --outFilterMismatchNoverLmax 0.04 \
         --alignIntronMin 20 \
         --alignIntronMax 1000000
         
    STAR_EXIT=$?
    
    # Cleanup
    rm -rf "/tmp/${SRR}_STARtmp"
    
    if (( STAR_EXIT == 0 )); then
        samtools index "${OUT_PREFIX}Aligned.sortedByCoord.out.bam"
    else
        echo "  STAR failed for $SRR. Exit code $STAR_EXIT."
    fi
done < "$MAP_FILE"

echo "--- Step 2: rMATS Configuration ---"
# DEPRECATED: rMATS now runs via 45_rmats_gse117357.sh
# The rMATS pipeline here has been removed to avoid duplication.

echo "Reprocessing Complete."
