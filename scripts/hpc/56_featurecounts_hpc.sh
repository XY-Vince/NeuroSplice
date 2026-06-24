#!/bin/bash
# ============================================================
# featureCounts — Gene Count Matrix (HPC / CANONICAL / TARGET)
# ============================================================
# Dataset:    GSE117357 (Adgrl3 KO, 60 samples, all male)
# Library:    SINGLE-END, reverse-stranded (-s 2)
# Aligner:    STAR v2.7.10b
# Annotation: Ensembl 102 (Mus_musculus.GRCm38.102.gtf)
# Environment: Yale HPC, Bouchet/Gibbs Cluster
# ============================================================

#SBATCH --job-name=featurecounts_gse117357
#SBATCH --partition=day
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/featurecounts_%j.out
#SBATCH --error=logs/featurecounts_%j.err

set -euo pipefail

module load subread

PROJECT="/gpfs/gibbs/project/neurosplice_gse117357"
BAM_DIR="${PROJECT}/alignments"
GTF="${PROJECT}/references/Mus_musculus.GRCm38.102.gtf"
OUT_DIR="${PROJECT}/featurecounts"
THREADS=8

mkdir -p "${OUT_DIR}"

shopt -s nullglob
BAMS=("${BAM_DIR}"/*.bam)

if (( ${#BAMS[@]} == 0 )); then
    echo "ERROR: No BAM files found in ${BAM_DIR}"
    exit 1
fi

echo "============================================================"
echo "featureCounts — Gene Count Matrix"
echo "============================================================"
echo "BAM directory:  ${BAM_DIR}"
echo "GTF annotation: ${GTF}"
echo "Output:         ${OUT_DIR}/gene_counts.txt"
echo "Samples:        ${#BAMS[@]}"
echo "Strand:         reverse (-s 2, confirmed via SRA metadata)"
echo "Library:        single-end (no -p flag)"
echo "Threads:        ${THREADS}"
echo "============================================================"

# KEY FLAGS:
#   -s 2        : reverse-stranded
#   -t exon     : count reads overlapping exons
#   -g gene_id  : summarize at gene level
#   -T ${THREADS}: threads
#   NO -p       : data is SINGLE-END (not paired)
#   -a          : Ensembl 102 GTF annotation
#   -o          : output file

featureCounts \
    -a "${GTF}" \
    -o "${OUT_DIR}/gene_counts.txt" \
    -s 2 \
    -t exon \
    -g gene_id \
    -T "${THREADS}" \
    --verbose \
    "${BAMS[@]}"

echo "featureCounts completed successfully"
