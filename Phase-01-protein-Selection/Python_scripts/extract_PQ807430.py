"""
FMDV Polyprotein Boundary Extractor for Target: PQ807430.1
Author: Md. Shakil Mahmud Supto
Phase: 01 - Protein Selection
Description: Extracts 14 functional viral sub-proteins based on specific structural 
             coordinates optimized for the PQ807430.1 strain dataset.
"""

import os
import urllib.request
import pandas as pd
from Bio import SeqIO

# 1. Path Management and Download Verification
file_id = "PQ807430.1"
fasta_file = f"../raw_neucliotide_wgs/{file_id}.fasta"

print(f"🌐 Fetching target {file_id} from NCBI Database...")
if not os.path.exists("../raw_neucliotide_wgs"):
    os.makedirs("../raw_neucliotide_wgs", exist_ok=True)

try:
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={file_id}&rettype=fasta"
    urllib.request.urlretrieve(url, fasta_file)
    print("✅ Genome file successfully verified!")
except Exception as e:
    print(f"❌ Network fallback issue: {e}. Attempting to proceed with local storage copy.")

if os.path.exists(fasta_file):
    record = SeqIO.read(fasta_file, "fasta")
    nucleotide_seq = record.seq

    # 2. Longest Open Reading Frame (ORF) Extraction Loop
    def find_longest_polyprotein(seq):
        longest_p = ""
        for strand, n_seq in [(1, seq), (-1, seq.reverse_complement())]:
            for frame in range(3):
                trans = n_seq[frame:].translate(to_stop=False)
                proteins = trans.split("*")
                for p in proteins:
                    if len(p) > len(longest_p):
                        longest_p = p
        return longest_p

    full_poly = find_longest_polyprotein(nucleotide_seq)
    p_str = str(full_poly)
    p_len = len(p_str)
    
    print(f"📊 Polyprotein Total Sequence Length: {p_len} AA")
    print("--- 🧬 Executing Strain-Specific Coordinate Mapping (coordinates_2) ---")

    # 3. Optimized Cleavage Structural Coordinates (Zero-Indexed Mapping)
    # Applied strictly based on the newly updated coordinates_2 block
    smart_coordinates = {
        "Leader_pro": (0, 201),
        "VP4": (201, 286),
        "VP2": (286, 504),
        "VP3": (504, 724),
        "VP1": (724, 935),
        "2A": (935, 951),
        "2B": (951, 1105),
        "2C": (1105, 1423),
        "3A": (1423, 1555),
        "3B_1": (1555, 1579),
        "3B_2": (1579, 1603),
        "3B_3": (1603, 1627),
        "3C_pro": (1627, 1840), 
        "3D_pol": (1840, p_len) 
    }

    # 4. Writing Multi-Fasta Profiles & Building Dataframes
    final_data = []
    out_dir = "../extracted_proteins"
    os.makedirs(out_dir, exist_ok=True)

    for p_name, (start, end) in smart_coordinates.items():
        p_seq = p_str[start:end]
        
        out_fasta = f"{out_dir}/{file_id}_{p_name}.fasta"
        with open(out_fasta, "w") as f:
            f.write(f">{file_id}_{p_name}\n{p_seq}\n")
            
        final_data.append({
            "Protein Name": p_name, "Start Pos (AA)": start + 1,
            "End Pos (AA)": end, "Length (AA)": len(p_seq), "Sequence": p_seq
        })

    # 5. Excel Table Export Processing
    df_final = pd.DataFrame(final_data)
    excel_file = f"{out_dir}/{file_id}_proteins_summary.xlsx"
    df_final.to_excel(excel_file, index=False)
    
    print(f"\n✅ Pipeline complete: 14 target functional sub-proteins mapped for {file_id}!")
    print(f"📁 Struct spreadsheet saved to path: {excel_file}")
    print("\n--- 📊 Subunit Feature Matrix Summary ---")
    print(df_final[["Protein Name", "Start Pos (AA)", "End Pos (AA)", "Length (AA)"]].to_string(index=False))
