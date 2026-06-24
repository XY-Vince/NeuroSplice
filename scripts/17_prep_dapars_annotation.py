import csv
import sys

# Input/Output
GTF_FILE = "reference/gencode.vM25.annotation.gtf"
BED_FILE = "reference/gencode.vM25.bed12"
SYMBOL_FILE = "reference/gencode.vM25.symbol_map.txt"

def parse_gtf():
    print("Parsing GTF...")
    transcripts = {}
    gene_map = {}
    
    with open(GTF_FILE, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            
            chrom = parts[0]
            feature = parts[2]
            start = int(parts[3]) - 1 # BED is 0-indexed
            end = int(parts[4])
            strand = parts[6]
            info = parts[8]
            
            # Extract ID
            tid = None
            gid = None
            gsym = None
            
            fields = info.strip().split(';')
            for field in fields:
                field = field.strip()
                if not field: continue
                if field.startswith('transcript_id'):
                    tid = field.split(' ')[1].strip('"')
                elif field.startswith('gene_id'):
                    gid = field.split(' ')[1].strip('"')
                elif field.startswith('gene_name'):
                    gsym = field.split(' ')[1].strip('"')
            
            if feature == 'transcript':
                if tid:
                    transcripts[tid] = {
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'exons': [],
                        'gid': gid
                    }
                if gid and gsym:
                    gene_map[tid] = gsym # Map Transcript ID to Gene Symbol because DaPars uses Transcript ID in BED
                    
            elif feature == 'exon':
                if tid and tid in transcripts:
                    transcripts[tid]['exons'].append((start, end))

    print(f"Loaded {len(transcripts)} transcripts.")
    
    # Write BED12
    print("Writing BED12...")
    with open(BED_FILE, 'w') as out:
        for tid, data in transcripts.items():
            chrom = data['chrom']
            strand = data['strand']
            exons = sorted(data['exons'])
            
            if not exons: continue
            
            tx_start = exons[0][0]
            tx_end = exons[-1][1]
            
            # BED12 fields
            # 1. chrom, 2. start, 3. end, 4. name, 5. score, 6. strand
            # 7. thickStart, 8. thickEnd, 9. itemRgb, 10. blockCount, 11. blockSizes, 12. blockStarts
            
            block_count = len(exons)
            block_sizes = ",".join([str(e[1] - e[0]) for e in exons])
            block_starts = ",".join([str(e[0] - tx_start) for e in exons])
            
            # thickStart/End usually coding region, but for DaPars standard usage whole transcript is often fine or CDS. 
            # We'll set thick = tx limits (non-coding styling) if we don't have CDS info, or just use tx limits.
            thick_start = tx_start
            thick_end = tx_end
            
            out.write(f"{chrom}\t{tx_start}\t{tx_end}\t{tid}\t0\t{strand}\t{thick_start}\t{thick_end}\t0\t{block_count}\t{block_sizes}\t{block_starts}\n")

    # Write Symbol Map
    print("Writing Symbol Map...")
    with open(SYMBOL_FILE, 'w') as out:
        for tid, gsym in gene_map.items():
            out.write(f"{tid}\t{gsym}\n")

if __name__ == "__main__":
    parse_gtf()
