#!/bin/bash
set -euo pipefail

# Activate STAR env
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate star_env

PROJECT_DIR="/Volumes/Untitled/NeuroSplice"
REF_DIR="${PROJECT_DIR}/reference"
STAR_INDEX_DIR="${REF_DIR}/star_index_sparse"
FASTQ_DIR="${PROJECT_DIR}/data/fastq"
OUT_DIR="${PROJECT_DIR}/aligned/star_pilot"
mkdir -p "${OUT_DIR}"

# Samples: 1 WT (SRR7540620), 1 KO (SRR7540650)
SAMPLES=("SRR7540620" "SRR7540650")

echo "--- Running STAR 2-Pass Pilot ---"

for SAMPLE in "${SAMPLES[@]}"; do
    echo "Aligning ${SAMPLE}..."
    
    # Check if FastQ exists
    R1="${PROJECT_DIR}/data/trimmed/${SAMPLE}_trimmed.fastq.gz"
    
    if [ ! -f "$R1" ]; then
        echo "Error: FastQ $R1 not found."
        continue
    fi
    
    # Decompress to temp file
    TMP_FASTQ="${OUT_DIR}/${SAMPLE}.fastq"
    echo "Decompressing to ${TMP_FASTQ}..."
    gunzip -c "${R1}" > "${TMP_FASTQ}"
    
    # Run STAR (Single End, Uncompressed Input)
    STAR --runThreadN 8 \
         --genomeDir "${STAR_INDEX_DIR}" \
         --readFilesIn "${TMP_FASTQ}" \
         --outFileNamePrefix "${OUT_DIR}/${SAMPLE}_" \
         --outSAMtype BAM SortedByCoordinate \
         --twopassMode Basic \
         --limitBAMsortRAM 4000000000
         
    # Cleanup
    rm "${TMP_FASTQ}"
    
    echo "Finished ${SAMPLE}"
done

echo "Pilot Complete. Check ${OUT_DIR}"
ls -lh "${OUT_DIR}"
