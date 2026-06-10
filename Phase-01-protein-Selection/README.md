# Phase 01: Protein Selection & Extraction Workflow 

This folder contains the complete pipeline for fetching, validating, and extracting target functional proteins from Foot-and-Mouth Disease Virus (FMDV) Whole Genome Sequences (WGS).

## 🔬 Methodology & Validation Pipeline

Our computational wet-lab-to-dry-lab validation followed these strict phases:
1. **WGS Acquisition**: Three core FMDV Whole Genome Sequences (WGS) were sourced from NCBI (Accession IDs: `PQ801122.1`, `PQ807429.1`, and `PQ807430.1`).
2. **ORF Identification**: Standard NCBI ORF Finder was deployed using the native FASTA sequences to detect and translate the longest stable Open Reading Frames (ORFs).
3. **GenPept Alignment & Cleavage Mapping**: To precisely locate the cleavage points for polyprotein processing, the retrieved ORFs were verified via BLASTp against high-confidence FMDV GenPept references. This locked the exact amino acid window lengths required for slicing.
4. **Automated Subunit Processing (Google Colab/Python)**: Custom, boundary-aware Python scripts were executed to automatically cut the polyprotein into 14 functional viral proteins (Leader_pro to 3D_pol) directly from the raw nucleotide strings.
5. **Post-Extraction BLAST Confirmation**: Every single isolated protein sequence was subjected to a secondary BLASTp search to confirm 100% sequence identity and structural validity before moving into the epitope screening phase.

---

## 📁 Directory Breakdown

* **`raw_neucliotide_wgs/`**: Contains the 3 original full-length nucleotide sequences in `.fasta` format.
* **`extracted_proteins/`**: Houses the generated sub-unit protein outputs (individual protein `.fasta` files and consolidated validation tables).
* **`Python_scripts/`**: Contains dedicated, isolated Python processing files mapped for each specific WGS Accession ID to ensure clear reproducibility.
