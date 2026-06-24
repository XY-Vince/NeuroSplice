#!/bin/bash
# Test STAR with its NATIVE --readFilesCommand option
# Using full paths to avoid any PATH issues

PROJECT_DIR="/Volumes/Untitled/NeuroSplice"
FASTQ_DIR="$PROJECT_DIR/data/trimmed"
ALIGNED_DIR="$PROJECT_DIR/data/aligned_star_native_test"
STAR_INDEX="$PROJECT_DIR/reference/mm10_STAR_index"
SRR="SRR7540620"
R1="$FASTQ_DIR/${SRR}_trimmed.fastq.gz"

STAR_BIN="/opt/miniconda3/envs/star_env/bin/STAR"

mkdir -p $ALIGNED_DIR
rm -rf /tmp/${SRR}_STARtmp

echo "Testing NATIVE readFilesCommand for $SRR..."
echo "STAR binary: $STAR_BIN"
echo "Input: $R1"

# Using STAR's native readFilesCommand with zcat
$STAR_BIN --genomeDir $STAR_INDEX \
    --readFilesIn "$R1" \
    --readFilesCommand /usr/bin/zcat \
    --outFileNamePrefix "$ALIGNED_DIR/${SRR}." \
    --outSAMtype BAM SortedByCoordinate \
    --outTmpDir /tmp/${SRR}_STARtmp \
    --runThreadN 8 \
    --outFilterMultimapNmax 10

echo "Exit Code: $?"
ls -lh $ALIGNED_DIR
cat "$ALIGNED_DIR/${SRR}.Log.final.out" | head -n 15
