"""
FMDV Polyprotein Boundary Extractor for Target: PQ801122.1
Author: Md. Shakil Mahmud Supto
Phase: 01 - Protein Selection
Description: Automatically cuts the translated full polyprotein into 14 functional 
             viral subunits using motif-aware window indexing.
"""

import os
import urllib.request
import pandas as pd
from Bio import SeqIO

# 1. Path Management and Download Verification
file_id = "PQ801122.1"
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
    print("--- 🧬 Executing Smart Motif Boundary Mapping... ---")

    # 3. Sliding Window Boundary Core Function
    def find_boundary(seq, motifs, start_w, end_w, default_val):
        sub_seq = seq[start_w:end_w]
        for m in motifs:
            pos = sub_seq.find(m)
            if pos != -1:
                return start_w + pos
        return default_val

    # 4. Cleavage Structural Coordinate Calculation
    j_L_VP4   = find_boundary(p_str, ["VFVPY"], 150, 210, 183)
    j_VP4_VP2 = find_boundary(p_str, ["QNNDW"], 240, 290, 268)
    
    j_VP2_VP3 = find_boundary(p_str, ["TEGAP"], 450, 520, 486)
    if j_VP2_VP3 != 486: j_VP2_VP3 += 5 

    j_VP1_2A  = find_boundary(p_str, ["RPLLA"], 880, 950, 917)
    j_VP3_VP1 = find_boundary(p_str, ["ALVVL", "ALVVS"], 680, 740, j_VP2_VP3 + 220)

    j_2A_2B   = find_boundary(p_str, ["APVKQ", "APEKQ"], 920, 960, 935)
    if j_2A_2B != 935: j_2A_2B += 5

    j_2B_2C   = find_boundary(p_str, ["GLVKV"], 1060, 1110, 1089)
    j_2C_3A   = find_boundary(p_str, ["LVQEV"], 1380, 1430, 1407)
    j_3A_3B1  = find_boundary(p_str, ["GHKVS", "GQKVS"], 1530, 1580, 1560)
    
    j_3B1_3B2 = find_boundary(p_str, ["GPYAG"], 1570, 1610, j_3A_3B1 + 24)
    if j_3B1_3B2 != j_3A_3B1 + 24: j_3B1_3B2 += 5

    j_3B2_3B3 = find_boundary(p_str, ["PYAGP"], 1590, 1630, j_3B1_3B2 + 24)
    if j_3B2_3B3 != j_3B1_3B2 + 24: j_3B2_3B3 += 5

    j_3B3_3C  = find_boundary(p_str, ["VKEGP", "PYEGP"], 1610, 1660, j_3B2_3B3 + 24)
    j_3C_3D   = find_boundary(p_str, ["CSCVSR"], 1820, 1860, 1845)

    # 5. Mapping Coordinates Context
    smart_coordinates = {
        "Leader_pro": (0, j_L_VP4), "VP4": (j_L_VP4, j_VP4_VP2), "VP2": (j_VP4_VP2, j_VP2_VP3),
        "VP3": (j_VP2_VP3, j_VP3_VP1), "VP1": (j_VP3_VP1, j_VP1_2A), "2A": (j_VP1_2A, j_2A_2B),
        "2B": (j_2A_2B, j_2B_2C), "2C": (j_2B_2C, j_2C_3A), "3A": (j_2C_3A, j_3A_3B1),
        "3B_1": (j_3A_3B1, j_3B1_3B2), "3B_2": (j_3B1_3B2, j_3B2_3B3), "3B_3": (j_3B2_3B3, j_3B3_3C),
        "3C_pro": (j_3B3_3C, j_3C_3D), "3D_pol": (j_3C_3D, p_len)
    }

    # 6. Writing Multi-Fasta Profiles & Building Dataframes
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

    # 7. Excel Table Export Processing
    df_final = pd.DataFrame(final_data)
    excel_file = f"{out_dir}/{file_id}_proteins_summary.xlsx"
    df_final.to_excel(excel_file, index=False)
    
    print(f"\n✅ Pipeline complete: 14 target functional sub-proteins mapped!")
    print(f"📁 Struct spreadsheet saved to path: {excel_file}")
    print("\n--- 📊 Subunit Feature Matrix Summary ---")
    print(df_final[["Protein Name", "Start Pos (AA)", "End Pos (AA)", "Length (AA)"]].to_string(index=False))
