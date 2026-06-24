#!/bin/bash
# Debug script for FIFO logic

PROJECT_DIR="/Volumes/Untitled/NeuroSplice"
FASTQ_DIR="$PROJECT_DIR/data/trimmed"
SRR="SRR7540620"
R1="$FASTQ_DIR/${SRR}_trimmed.fastq.gz"
FIFO="/tmp/${SRR}.debug.fifo"

echo "Checking input file: $R1"
ls -l "$R1"

echo "Creating FIFO $FIFO"
rm -f "$FIFO"
mkfifo "$FIFO"

echo "Starting gzip in background..."
# Redirect stderr to a log file to see errors
/usr/bin/gzip -dc "$R1" > "$FIFO" 2> /tmp/gzip_error.log &
GZIP_PID=$!

echo "Reading from FIFO (first 10 lines)..."
head -n 10 "$FIFO"
READ_EXIT=$?

echo "Head exit code: $READ_EXIT"
echo "Gzip PID: $GZIP_PID"

wait $GZIP_PID
GZIP_EXIT=$?
echo "Gzip exit code: $GZIP_EXIT"

echo "Gzip stderr:"
cat /tmp/gzip_error.log

rm -f "$FIFO"
