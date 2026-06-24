#!/bin/bash

# Create Directory
mkdir -p data/gse173926/raw
cd data/gse173926/raw

# Download using curl
echo "Starting download of 12 samples (24 FASTQ files) from GSE173926..."
echo "URLs loaded from data/metadata/GSE173926_urls.txt"

# Check if file exists
if [ ! -f ../../metadata/GSE173926_urls.txt ]; then
    echo "URL file missing!"
    exit 1
fi

# Use xargs to download in parallel (upto 4 processes)
cat ../../metadata/GSE173926_urls.txt | xargs -n 1 -P 4 curl -O -L -f -s

echo "Download complete."
