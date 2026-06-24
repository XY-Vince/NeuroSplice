#!/bin/bash
# ============================================================
# SLURM wrapper for DESeq2 R script
# ============================================================
#SBATCH --job-name=deseq2_gse117357
#SBATCH --partition=day
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=logs/deseq2_%j.out
#SBATCH --error=logs/deseq2_%j.err

set -euo pipefail

# Load required modules on Yale HPC
module load R/4.2.0-foss-2020b

# Run the R script
echo "Starting DESeq2 Analysis..."
Rscript 58_deseq2_hpc.R
echo "DESeq2 Analysis Complete!"
