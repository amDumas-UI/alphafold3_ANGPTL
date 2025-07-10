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

---

## Purpose
This toolkit analyzes interaction interfaces predicted by AlphaFold, focusing on:

- Chain-specific interaction mapping (e.g., chain A/B interacting with chain C)
- Lipase pocket analysis (residues in chain D)
- Prediction stability via Jaccard similarity
- Comparative modeling for different trimeric configurations

---

## Features

- **Barpolots** of intereface residue occurrence across multiple AlphaFold models
- **Pocket residue tracking** based on lipase family homology
- **Automated AlphaFold submissions** via Playwright scripting
- **Jaccard similarity metrics** to compare predictions between cutoffs or shuffled chains
- **Color-coded residues** by chain and importance

---

## Usage

### 1. Submit Sequences to AlphaFold

'''bash
python heterotrimer_alphafold_submitter_playwright.py -- input_fasta your_input.fasta --output_dir alphafold_results/
'''

### 2. Run Interface Analysis
#### For heterotrimer chains A/B vs C:
'''bash
python full_trimer_interactions.py
'''

####For heterotrimer chains A/B/C vs Lipase chain D
'''bash
python lipase_interactions.py
'''
#### For lipase pocket residues (chain D):
'''bash
python lipase_pocket.py
'''
