#!/bin/bash
set -euo pipefail

# Activate STAR env
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate star_env

PROJECT_DIR="/Volumes/Untitled/NeuroSplice"
REF_DIR="${PROJECT_DIR}/reference"
STAR_INDEX_DIR="${REF_DIR}/star_index_sparse"
MKDIR_P="mkdir -p"
$MKDIR_P "${STAR_INDEX_DIR}"

GNOME_FA="${REF_DIR}/mm10.fa"
GTF="${REF_DIR}/gencode.vM25.annotation.gtf"

# STAR Indexing (Sparse Mode for 16GB RAM)
# --genomeSAindexNbases 12 (Down from default 14) drastically reduces RAM usage
# --genomeChrBinNbits 14 (Can also help)

echo "--- Building STAR Sparse Index (16GB RAM Limit) ---"
echo "Genome: ${GNOME_FA}"
echo "GTF: ${GTF}"
echo "Output: ${STAR_INDEX_DIR}"

STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir "${STAR_INDEX_DIR}" \
     --genomeFastaFiles "${GNOME_FA}" \
     --sjdbGTFfile "${GTF}" \
     --sjdbOverhang 99 \
     --genomeSAindexNbases 12 \
     --genomeChrBinNbits 14 \
     --limitGenomeGenerateRAM 15000000000 # 15GB Safety Cap

echo "Index Generation Complete."
du -sh "${STAR_INDEX_DIR}"
