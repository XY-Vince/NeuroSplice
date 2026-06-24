
import requests
import xml.etree.ElementTree as ET
import time
import sys

def get_srr_for_gsm(gsm_id):
    search_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term={gsm_id}&retmax=10"
    try:
        r = requests.get(search_url)
        root = ET.fromstring(r.text)
        ids = [id_node.text for id_node in root.findall(".//Id")]
        if not ids:
            return []
        
        # Requesting rettype=runinfo often returns XML by default in these eutils
        fetch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id={','.join(ids)}&rettype=runinfo"
        r = requests.get(fetch_url)
        srrs = []
        # Parse XML to find <Run> tags
        try:
            fetch_root = ET.fromstring(r.text)
            for run_node in fetch_root.findall(".//Run"):
                srrs.append(run_node.text)
        except ET.ParseError:
            # Fallback for CSV-like output if XML parsing fails
            for line in r.text.split('\n'):
                if line.startswith('SRR'):
                    srrs.append(line.split(',')[0])
        return srrs
    except Exception as e:
        print(f"Error fetching {gsm_id}: {e}", file=sys.stderr)
        return []

if __name__ == "__main__":
    with open("data/metadata/gsm_accessions.txt", "r") as f:
        gsms = [line.strip() for line in f if line.strip()]
    
    with open("data/metadata/srr_accessions.txt", "w") as out:
        for gsm in gsms:
            print(f"Processing {gsm}...", file=sys.stderr)
            srrs = get_srr_for_gsm(gsm)
            for srr in srrs:
                out.write(srr + "\n")
            time.sleep(0.5) # Avoid hitting NCBI rate limits
    print("Done!")
