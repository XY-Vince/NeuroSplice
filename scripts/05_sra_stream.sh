#!/bin/bash
# 05_sra_stream.sh - Process SRA files as they finish downloading

DATA_DIR="/Volumes/Untitled/NeuroSplice"
RAW_DIR="${DATA_DIR}/data/raw"
SRA_ROOT="${DATA_DIR}" # prefetch puts folders here

mkdir -p "${RAW_DIR}"

echo "[$(date)] Starting SRA streaming processor..."

while true; do
    # Find SRR folders that have a .sra file and ARE NOT being actively downloaded
    # We check if the .sra file size hasn't changed in the last 60 seconds
    for srr_dir in ${SRA_ROOT}/SRR[0-9]*; do
        if [[ ! -d "${srr_dir}" ]]; then continue; fi
        SRR=$(basename "${srr_dir}")
        SRA_FILE="${srr_dir}/${SRR}.sra"
        
        if [[ -f "${SRA_FILE}" ]]; then
            # Check if already processed
            if [[ -f "${RAW_DIR}/${SRR}_1.fastq.gz" ]] || [[ -f "${RAW_DIR}/${SRR}.fastq.gz" ]]; then
                # Already converted to FASTQ, we can remove the SRA folder to save space
                echo "[$(date)] ${SRR} already converted. Cleaning up SRA folder..."
                rm -rf "${srr_dir}"
                continue
            fi
            
            # Check if file size is stable (not being written to)
            size1=$(stat -f %z "${SRA_FILE}")
            sleep 5
            size2=$(stat -f %z "${SRA_FILE}")
            
            if [[ "${size1}" == "${size2}" ]] && [[ "${size1}" -gt 1000000 ]]; then
                echo "[$(date)] ${SRR} download complete (${size1} bytes). Converting to FASTQ..."
                
                source /opt/miniconda3/bin/activate combio
                
                # Convert SRA to FASTQ
                fasterq-dump --outdir "${RAW_DIR}" --split-3 --threads 4 "${SRA_FILE}"
                
                if [[ $? -eq 0 ]]; then
                    echo "[$(date)] ${SRR} conversion successful. Compressing..."
                    # Find what files were created (single or paired)
                    if [[ -f "${RAW_DIR}/${SRR}_1.fastq" ]]; then
                        gzip "${RAW_DIR}/${SRR}_1.fastq"
                        [[ -f "${RAW_DIR}/${SRR}_2.fastq" ]] && gzip "${RAW_DIR}/${SRR}_2.fastq"
                    elif [[ -f "${RAW_DIR}/${SRR}.fastq" ]]; then
                        gzip "${RAW_DIR}/${SRR}.fastq"
                    fi
                    
                    echo "[$(date)] ${SRR} ready for preprocessing. Triggering script..."
                    bash scripts/03_preprocess.sh --sample "${SRR}"
                    
                    # Cleanup SRA
                    rm -rf "${srr_dir}"
                else
                    echo "[ERROR] fasterq-dump failed for ${SRR}"
                fi
            fi
        fi
    done
    
    echo "[$(date)] Sleeping for 2 minutes..."
    sleep 120
done
