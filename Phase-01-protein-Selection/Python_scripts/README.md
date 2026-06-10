# Python Automation Scripts 

This directory contains standalone Python scripts developed to automate the full-length Open Reading Frame (ORF) 
translation and polyprotein subunit extraction for FMDV genome targets.

## Prerequisites & Installation

Before running these scripts, you must install the required biological and data-processing libraries. 
Run the following command in your terminal or a Jupyter/Google Colab notebook cell:

```bash
pip install biopython pandas openpyxl

(Note: Once these libraries are installed, you don't need to download anything manually. The scripts are
fully automated to fetch the correct WGS genome files directly from the NCBI database and extract
all 14 functional protein subunits into the project folders.)

