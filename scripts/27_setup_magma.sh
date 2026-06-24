#!/bin/bash
set -euo pipefail

# Directories
PROJECT_DIR="/Volumes/Untitled/NeuroSplice"
MAGMA_DIR="${PROJECT_DIR}/results/magma"
mkdir -p "${MAGMA_DIR}"

cd "${MAGMA_DIR}"

echo "--- 1. Downloading MAGMA ---"
if [ ! -f "magma" ]; then
    echo "Downloading MAGMA zip..."
    curl -L -O https://ctg.cncr.nl/software/MAGMA/program/magma_v1.10.zip
    unzip -o magma_v1.10.zip
    chmod +x magma
    rm magma_v1.10.zip
fi

echo "--- 2. Downloading Reference Panel (1000 Genomes EUR) ---"
if [ ! -f "g1000_eur.bim" ]; then
    echo "Downloading Reference zip..."
    curl -L -O https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip
    unzip -o g1000_eur.zip
    rm g1000_eur.zip
fi

echo "--- 3. Downloading ADHD GWAS (Demontis 2019) ---"
if [ ! -f "adhd_eur_jul2017.gz" ]; then
    echo "Downloading GWAS stats..."
    curl -L -O https://ctg.cncr.nl/documents/p1651/adhd_eur_jul2017.gz
fi

if [ ! -f "NCBI37.3.gene.loc" ]; then
    echo "Downloading Gene Locations..."
    curl -L -O https://ctg.cncr.nl/software/MAGMA/aux_files/NCBI37.3.zip
    unzip -o NCBI37.3.zip
    rm NCBI37.3.zip
fi

echo "Setup Complete."
ls -lh "${MAGMA_DIR}"
