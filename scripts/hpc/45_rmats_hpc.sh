#!/bin/bash
# ============================================================
# rMATS Alternative Splicing (HPC / CANONICAL / COMPLETED)
# ============================================================
# Dataset:    GSE117357 (Adgrl3 KO, 60 samples)
# Library:    SINGLE-END, 76bp, fr-firststrand
# Tool:       rMATS v4.3.0
# Annotation: Ensembl 102 (Mus_musculus.GRCm38.102.gtf)
# Environment: Yale HPC, Bouchet/Gibbs Cluster
# 
# NOTE: This script serves as the authoritative reproducibility
# record for the splicing runs that generated the validated CSVs.
# This step has already completed on the HPC.
#
# CRITICAL PARAMETER: --allow-clipping is required for 76bp 
# single-end reads to rescue short overhang junctions.
# ============================================================

#SBATCH --job-name=rmats_gse117357
#SBATCH --partition=day
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=logs/rmats_%j.out
#SBATCH --error=logs/rmats_%j.err

set -euo pipefail

module load miniconda
conda activate rmats

PROJECT="/gpfs/gibbs/project/neurosplice_gse117357"
BAM_LISTS="${PROJECT}/rmats_input"
GTF="${PROJECT}/references/Mus_musculus.GRCm38.102.gtf"
RMATS_OUT="${PROJECT}/rmats_out"
THREADS=16

mkdir -p "${RMATS_OUT}"

for tissue in hippocampus prefrontal striatum; do
    echo "Running rMATS for ${tissue}"
    OUT_DIR="${RMATS_OUT}/${tissue}"
    mkdir -p "${OUT_DIR}"
    
    rmats.py \
        --b1 "${BAM_LISTS}/wt_${tissue}_bams.txt" \
        --b2 "${BAM_LISTS}/ko_${tissue}_bams.txt" \
        --gtf "${GTF}" \
        --od "${OUT_DIR}" \
        --tmp "${OUT_DIR}/tmp" \
        -t single \
        --libType fr-firststrand \
        --novelSS 1 \
        --allow-clipping \
        --readLength 76 \
        --variable-read-length \
        --nthread "${THREADS}" \
        --cstat 0.0001 \
        --tstat 6
        
    echo "Completed rMATS for ${tissue}"
done
