import os
import glob
import sys

def clean_gtf(filepath):
    """Reads a GTF, removes lines with end < start OR format errors."""
    temp_path = filepath + ".tmp"
    removed_count = 0
    total_count = 0
    
    with open(filepath, 'r') as fin, open(temp_path, 'w') as fout:
        for line in fin:
            total_count += 1
            if line.startswith('#'):
                fout.write(line)
                continue
            
            # Strict TAB separation
            parts = line.strip().split('\t')
            
            # Must have at least 8 columns (attributes can be last)
            if len(parts) < 8:
                removed_count += 1
                continue
                
            # Check Feature (col 2) - specific malformed line check
            feature = parts[2]
            if "transcript6" in feature or "exon6" in feature: # Catch the specific corruption
                removed_count += 1
                continue

            try:
                start = int(parts[3])
                end = int(parts[4])
                
                if end < start:
                    removed_count += 1
                    continue
                    
                # Extra check: Start/End sensible?
                if start < 1 or end < 1:
                    removed_count += 1
                    continue

                fout.write(line)
            except ValueError:
                removed_count += 1
                continue
    
    if removed_count > 0:
        print(f"cleaned {os.path.basename(filepath)}: removed {removed_count}/{total_count} lines")
        os.replace(temp_path, filepath)
    else:
        os.remove(temp_path)

if __name__ == "__main__":
    gtf_dir = "results/stringtie_tss"
    if not os.path.exists(gtf_dir):
        print(f"Directory {gtf_dir} not found.")
        sys.exit(1)
        
    gtfs = glob.glob(os.path.join(gtf_dir, "SRR*.gtf"))
    print(f"Scrubbing {len(gtfs)} GTF files in {gtf_dir} (Strict Mode)...")
    
    for gtf in gtfs:
        clean_gtf(gtf)
    print("Strict cleanup complete.")
