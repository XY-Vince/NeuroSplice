#!/bin/bash
# Transfer script for Bouchet HPC setup
# Using scratch directory due to home quota limit (126GB < 155GB needed)

NETID="xw539"
SSH_KEY="$HOME/.ssh/id_rsa_bouchet"
SCRATCH_DIR="/nfs/roberts/scratch/pi_km585/xw539/NeuroSplice"

echo "=== Bouchet Data Transfer Script ==="
echo "Transferring to: ${NETID}@bouchet.ycrc.yale.edu"
echo "Destination: ${SCRATCH_DIR}"
echo ""

# Create directories on Bouchet
echo "Creating directories..."
ssh -i ${SSH_KEY} ${NETID}@bouchet.ycrc.yale.edu "mkdir -p ${SCRATCH_DIR}/data/{trimmed,aligned_star_reprocess} ${SCRATCH_DIR}/reference ${SCRATCH_DIR}/scripts ${SCRATCH_DIR}/logs"

# Transfer FASTQ files (exclude already-completed samples)
echo "Transferring FASTQ files (~130GB)..."
rsync -avP -e "ssh -i ${SSH_KEY}" --exclude='SRR7540620*' --exclude='SRR7540621*' \
    /Volumes/Untitled/NeuroSplice/data/trimmed/*.fastq.gz \
    ${NETID}@bouchet.ycrc.yale.edu:${SCRATCH_DIR}/data/trimmed/

# Transfer STAR index (25GB)
echo "Transferring STAR index..."
rsync -avP -e "ssh -i ${SSH_KEY}" \
    /Volumes/Untitled/NeuroSplice/reference/mm10_STAR_index \
    ${NETID}@bouchet.ycrc.yale.edu:${SCRATCH_DIR}/reference/

# Transfer SLURM script
echo "Transferring job script..."
rsync -avP -e "ssh -i ${SSH_KEY}" \
    /Volumes/Untitled/NeuroSplice/scripts/bouchet_star_array.slurm \
    ${NETID}@bouchet.ycrc.yale.edu:${SCRATCH_DIR}/scripts/

echo ""
echo "=== Transfer Complete! ==="
echo "Next steps:"
echo "1. SSH to Bouchet: ssh -i ${SSH_KEY} ${NETID}@bouchet.ycrc.yale.edu"
echo "2. Submit job: cd ${SCRATCH_DIR}/scripts && sbatch bouchet_star_array.slurm"
