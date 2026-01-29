import requests
import time
import os

def fetch_pubchem_compounds(term="investigational[Filter]", max_records=500, output_file="data/raw/pubchem_real.sdf"):
    """
    Fetches real SDF data from PubChem using PUG REST and E-Utilities.
    """
    print(f"🔍 Searching PubChem for '{term}' (Max: {max_records})...")
    
    # 1. Search for CIDs using E-Utils
    esearch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "pccompound",
        "term": term,
        "retmax": max_records,
        "retmode": "json"
    }
    
    try:
        r = requests.get(esearch_url, params=params)
        r.raise_for_status()
        data = r.json()
        cids = data["esearchresult"]["idlist"]
        print(f"✅ Found {len(cids)} CIDs. Downloading structures...")
    except Exception as e:
        print(f"❌ Search failed: {e}")
        return

    # 2. Download SDF in batches (PubChem limits URL length)
    batch_size = 50
    with open(output_file, "w", encoding="utf-8") as f_out:
        for i in range(0, len(cids), batch_size):
            batch = cids[i:i+batch_size]
            cid_str = ",".join(batch)
            
            # PUG REST URL
            pug_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid_str}/SDF"
            
            try:
                r_sdf = requests.get(pug_url)
                if r_sdf.status_code == 200:
                    f_out.write(r_sdf.text)
                    print(f"   Downloaded batch {i//batch_size + 1}/{(len(cids)//batch_size)+1}")
                else:
                    print(f"   ⚠️ Failed batch {i}: {r_sdf.status_code}")
                
                # Polite delay
                time.sleep(0.5) 
            except Exception as e:
                print(f"   ❌ Error downloading batch: {e}")

    print(f"🎉 Done! Saved {len(cids)} compounds to {output_file}")

if __name__ == "__main__":
    if not os.path.exists("data/raw"):
        os.makedirs("data/raw")
    # Using a broad filter to get diverse real chemistry
    fetch_pubchem_compounds(term="drug", max_records=500)
