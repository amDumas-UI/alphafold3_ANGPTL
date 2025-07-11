# AlphaFold ANGPTL-Lipase Interaction Pipeline

**Author:** Alexander Dumas
**Last Updated:** July 2025

## Overview
This repository contains scripts to analyze AlphaFold3 predictions of Angiopoietin like proteins and their ability to interact with each other and members of the pancreatic lipase family

Key features include:
- Automated AlphaFold submission using Playwright (semi-automated for heterotrimer analysis)
- Extraction of interacting residues using PyMOL and distance-based filtering.
- Bar plot generation of residue-level interaction frequencies across multiple models.
- Highlighting of conserved or catalytic residues across multiple lipases.
- Jaccard similarity calculations to assess prediction consistency.

---

## Repository structure
.
├── combinations/
│   ├── ht_combinations.txt            # Heterotrimer combinations to be modeled
│   └── combinations.txt               # Homotrimer/lipase combinations to be modeled
│
├── scripts/
│   ├── slphafold_submitter.py         # Submission pipeline for single-lipase models
│   |── heterotrimer_submitter.py      # Submission pipeline for heterotrimeric models
│   ├── batch_interface.py             # Extracts interface residues (homotrimers)
│   ├── ht_batch_interface.py          # Extracts interface residues (heterotrimers)
│   ├── full_trimer_interaction.py     # A/B chains interacting with chain C (barplots)
│   ├── lipase_pocket.py               # Lipase (chain D) interacting with A/B/C (barplots)
│   ├── lipase_interactions.py         # Alternative method to extract D-chain residues
│   ├── spatial_proximity.py           # Jaccard similarity between models/cutoffs (to determine proper surface area cutoff in ASA calculations
│   └── interfaceResidues.py           # ASA-based interface residue detection (PyMOL script)
│
├── interfaces/                        # Saved interface .pdb files extracted from models
│   └── [heterotrimer]                 # Organized by all ANGPTL heterotrimer combinations tested
|        └── [lipase folders]          # Organized by combination/cutoff/lipase
|    └── [homotrimer]
|        └── [lipase folders]          # Organized by combination/cutoff/lipase
│
├── plots/                             # PNG barplots, heatmaps, Jaccard curves
│   └── [analysis condition folders]   # Matched to interface folders
│
├── alphafold_results/                 # Saved .zip archives with .cif prediction models from AlphaFold
│   └── [heterotrimer]                 # Organized by all ANGPTL heterotrimer combinations tested
    └── [homotrimer]                   # Organized by all ANGPTL homotrimer combinations tested
│
|
├── sequences/                         # list of .txt files containing fasta sequences for ANGPTL and lipase proteins
|
├── README.md                          # This file

---

##Purpose
This toolkit analyzes residue-level interaction between ANGPTL trimers and lipase proteins including:
- Chain-specific interaction mapping (e.g., chain A/B interacting with chain C)
- Lipase pocket analysis (residues in chain D)
- Assessing prediction stability via Jaccard similarity
- Comparative modeling for different trimeric configurations (homotrimers vs heterotrimer)

---

##Features

- **Barpolots** of intereface residue occurrence across multiple AlphaFold models
- **Pocket residue tracking** based on lipase family homology
- **Automated AlphaFold submissions** via Playwright scripting
- **Jaccard similarity metrics** to compare predictions between cutoffs or shuffled chains
- **Color-coded residues** by chain and importance

---

##Usage

###1. Submit Sequences to AlphaFold (While not required, this was used to automate submissions and be more efficient)

'''bash
####Heterotrimeric complex
python heterotrimer_alphafold_submitter_playwright.py -- input_fasta path/to/sequence/.txt -- output_cif alphafold_results/heterotrimer/.cif

####Homotrimeric complex
python alphafold_submitter_playwright.py -- input_fasta path/to/sequence/.txt -- output_cif alphafold_results/homotrimer/
'''

###2. Extract Interacting Residues based on ASA calculation

'''bash
####Heterotrimer
python ht_batch_interface.py -- input_dir alphafold_results/heterotrimer/combination/  --output_pdb interface/heterotrimer/respective combination/respective lipase/.pdb

####Homotrimer
python batch_interface.py -- input_dir alphafold_results/homotrimer/combination/ -- output_pdb interface/homotrimer/respective combination/respective lipase/.pdb
'''

###2. Run Interface Analysis

####For chains A/B vs C:
'''bash
python full_trimer_interactions.py -- input_dir alphafold_results/____trimer/ -- output_png plots/______trimer/respective combination/.png
'''

####For heterotrimer chains A/B/C vs Lipase chain D
'''bash
python lipase_interactions.py -- input_cif_dir input_pdb_dir Path/to/parent/directory/for/cif_predictions Path/to/parent/directory/for/pdb_interface_predictions -- output_png plots/______trimer/respective combination/.png
'''

####For lipase pocket residues (chain D):
'''bash
python lipase_pocket.py -- input_cif_dir input_pdb_dir Path/to/parent/directory/for/cif_predictions Path/to/parent/directory/for/pdb_interface_predictions -- output_png plots/______trimer/respective combination/.png
'''

##Dependencies
- PyMOL
- Biopython
- matplotlib
- pandas
- seaborn
- tqdm (for progress bars)
- Python >= 3.8

- install via conda
    - conda install -c conda-forge biopython pymol matplotlib pandas seaborn tqdm
