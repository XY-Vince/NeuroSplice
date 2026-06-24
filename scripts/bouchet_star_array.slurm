#!/bin/bash
#SBATCH --job-name=star_align
#SBATCH --partition=day
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --array=0-57
#SBATCH --output=logs/star_%A_%a.out
#SBATCH --error=logs/star_%A_%a.err

# STAR Alignment Array Job for NeuroSplice Project
# Submits 58 remaining samples (SRR7540622 - SRR7540679)
# SRR7540620 and SRR7540621 already completed locally

# Load modules
module load STAR/2.7.11b-GCC-13.3.0
module load SAMtools/1.21-GCC-13.3.0

# Paths (using scratch to avoid home quota limit)
PROJECT_DIR="/nfs/roberts/scratch/pi_km585/xw539/NeuroSplice"
FASTQ_DIR="${PROJECT_DIR}/data/trimmed"
ALIGNED_DIR="${PROJECT_DIR}/data/aligned_star_reprocess"
STAR_INDEX="${PROJECT_DIR}/reference/mm10_STAR_index"

# Create output directory
mkdir -p ${ALIGNED_DIR}
mkdir -p logs

# Sample list (remaining 58 samples)
SAMPLES=(
SRR7540622 SRR7540623 SRR7540624 SRR7540625 SRR7540626 SRR7540627 SRR7540628 SRR7540629
SRR7540630 SRR7540631 SRR7540632 SRR7540633 SRR7540634 SRR7540635 SRR7540636 SRR7540637
SRR7540638 SRR7540639 SRR7540640 SRR7540641 SRR7540642 SRR7540643 SRR7540644 SRR7540645
SRR7540646 SRR7540647 SRR7540648 SRR7540649 SRR7540650 SRR7540651 SRR7540652 SRR7540653
SRR7540654 SRR7540655 SRR7540656 SRR7540657 SRR7540658 SRR7540659 SRR7540660 SRR7540661
SRR7540662 SRR7540663 SRR7540664 SRR7540665 SRR7540666 SRR7540667 SRR7540668 SRR7540669
SRR7540670 SRR7540671 SRR7540672 SRR7540673 SRR7540674 SRR7540675 SRR7540676 SRR7540677
SRR7540678 SRR7540679
)

# Get current sample from array index
SRR=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
R1="${FASTQ_DIR}/${SRR}_trimmed.fastq.gz"
OUT_PREFIX="${ALIGNED_DIR}/${SRR}."

echo "Processing sample: ${SRR} (Array ID: ${SLURM_ARRAY_TASK_ID})"
echo "Input: ${R1}"
echo "Output prefix: ${OUT_PREFIX}"

# Check if input exists
if [[ ! -f "${R1}" ]]; then
    echo "ERROR: Input file not found: ${R1}"
    exit 1
fi

# Check if already completed
if [[ -f "${OUT_PREFIX}Aligned.sortedByCoord.out.bam" ]]; then
    SIZE=$(stat -c%s "${OUT_PREFIX}Aligned.sortedByCoord.out.bam")
    if (( SIZE > 100000000 )); then
        echo "Already aligned (${SIZE} bytes). Skipping."
        exit 0
    fi
fi

# Run STAR
STAR --genomeDir ${STAR_INDEX} \
     --readFilesIn ${R1} \
     --readFilesCommand zcat \
     --outFileNamePrefix ${OUT_PREFIX} \
     --outSAMtype BAM SortedByCoordinate \
     --runThreadN ${SLURM_CPUS_PER_TASK} \
     --outSAMstrandField intronMotif \
     --outFilterMultimapNmax 10 \
     --alignSJoverhangMin 8 \
     --alignSJDBoverhangMin 1 \
     --outFilterMismatchNmax 999 \
     --outFilterMismatchNoverLmax 0.04 \
     --alignIntronMin 20 \
     --alignIntronMax 1000000 \
     --alignMatesGapMax 1000000 \
     --limitBAMsortRAM 50000000000

STAR_EXIT=$?

if (( STAR_EXIT == 0 )); then
    echo "Indexing BAM..."
    samtools index ${OUT_PREFIX}Aligned.sortedByCoord.out.bam
    echo "Done: ${SRR}"
else
    echo "STAR failed for ${SRR}. Exit code: ${STAR_EXIT}"
    exit ${STAR_EXIT}
fi
