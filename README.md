# AlphaFold ANGPTL-Target Interaction Pipeline

**Author:** Alexander Dumas
**Last Updated:** August 2025

## Project Abstract
The angiopoietin-like protein (ANGPTL) family consists of 8 members which are structurally similar to angiopoietins. ANGPTL3, ANGPTL4, and ANGPTL8 play important roles in regulating lipoprotein metabolism by inhibiting lipoprotein lipase (LPL) and endothelial lipase (EL). While the atomic structures of ANGPTL proteins have yet to be solved, ANGPTL3 and ANGPTL4 have been shown to form homotrimers, and ANGPTL3 and ANGPTL8 have been shown to form a 2:1 heterotrimer. We have found that formation of these oligomeric structures is essential for inhibiting lipases. Is has been suggested that ANGPTL4 might inhibit other lipases, including hepatic lipase, pancreatic lipase, and lysosomal acid lipase, but these interactions have not been strongly validated. Moreover, other, less studied, ANGPTL family members such as ANGPTL5 and ANGPTL6 have been implied to be associated with triglyceride metabolism in humans. However, their interactions with lipases and with other ANGPTL proteins have not been characterized. Using AlphaFold3 as a prediction tool to explore these interactions could highlight probable interacting partners within the ANGPTL family and predict an interaction with target lipases. To investigate this idea, I developed a Python and PyMOL-based pipeline using AlphaFold3 to predict ANGPTL complex formation as well as identify interactions between predicted ANGPTL trimers and target lipases. Using a method that calculates the interacting residues from multiple predicted models, I show that AlphaFold3 can successfully predict the interactions between experimentally validated ANGPTL oligomers and the lipases they interact with, including the interaction of ANGPTL4 homotrimers and ANGPTL3/ANGPTL8 heterotrimers with LPL. I also show that AlphaFold3 consistently predicts that ANGPTL4 and ANGPTL5 can form a heterotrimer that interacts within the catalytic pocket of EL akin to the interaction ANGPTL3 homotrimers are known to have with EL. 

## Overview
This repository contains scripts to analyze AlphaFold3 predictions of ANGPTL proteins and their ability to interact with each other and target proteins such as members of the pancreatic lipase family.

##Key features include:
- Automated AlphaFold submission using Playwright (semi-automated for heterotrimer analysis)
      - Note: This utilizes Google DeepMind and Isomorphic labs' AlphaFold Server which is powered by AlphaFold3.
- Jaccard similarity calculations to assess prediction consistency (for SASA cutoff determination).
- Extraction of interacting residues using PyMOL and distance-based filtering using solvent-accessible surface area (∆SASA) calculations.
- Bar plot generation of residue-level interaction frequencies across multiple models.
- Highlighting of conserved or catalytic residues across multiple target lipases.
- Highlighting consensus residues onto the ANGPTL and target models in a PYMOL environment.


---

## Repository structure
 
    ``` 
    . 
    ├── combinations/ 
    │     ├── ht_combinations.txt              # Heterotrimer/target combinations to be modeled 
    │     └── combinations.txt                 # Homotrimer/target combinations to be modeled 
    |
    ├── scripts/ 
    │     ├── ANGPTL_alphafold_submitter_playwright.py        # AlphaFold submission  
    │     ├── batch_interface.py               # Extracts interface residues 
    │     ├── full_trimer_interaction.py       # Barplots: A/B vs C chains 
    |     ├── map_consensus_residues.py        # Maps interaction residues back to ANGPTL and target structures
    │     ├── lipase_pocket.py                 # Barplots: Target lipase interactions 
    │     ├── lipase_interactions.py           # Barplots: ANGPTL interactions with target lipase 
    │     ├── spatial_proximity.py             # Jaccard similarity comparisons 
    │         └── interfaceResidues.py         # PyMOL interface script (∆SASA calculation) 
    │ 
    ├── alphafold_results/                     # .zip results from AlphaFold3 predictions 
    │     ├── [ANGPTL combinations]            # Sorted by ANGPTL complex structure and target interaction
    │
    ├── interfaces/                            # PDBs containing interface residues from ∆SASA calculations  
    │     ├── [ANGPTL combinations]            # Sorted by ANGPTL complex structure and target interaction         
    │ 
    ├── plots/                                 # Plots: barplots, Jaccard curves, PYMOL .pse environments 
    │     └── [matching to interface folders]
    │ 
    ├── sequences/                             # FASTA files (.txt) for ANGPTL and targets 
    │ 
    └── README.md                              
    ``` 

---

##Usage

###1. Submit Sequences to AlphaFold (Automation not required, this was used to automate submissions and be more efficient)

'''bash
####ANGPTL + Target Interaction Combinations
python ANGPTL_alphafold_submitter_playwright.py -- input_seq_dir combinations/path/to/sequence/directory/ -- output_cif alphafold_results/path/to/.cif


###2. Extract Interacting Residues based from ∆SASA Calculations

'''bash
####ANGPTL + Target: Interaction Prediction Models
python batch_interface.py -- input_cif_dir alphafold_results/path/to/cif/directory/ -- output_pdb interface/path/to/.pdb
'''

###3. ANGPTL Complex Consistency Analysis

####For chains A/B vs C:
'''bash
python full_trimer_interactions.py -- input_cif_dir alphafold_results/path/to/cif/directory/ -- output_png plots/path/to/.png
'''

###4. ANGPTL-Target Interaction Consistency Analysis

####For ANGPTL chains A/B/C vs target chain D
'''bash
python lipase_interactions.py -- input_cif_dir input_pdb_dir alphafold_results/path/to/cif/directory/ interface/path/to/pdb/directory -- output_png plots/path/to/.png
'''

####For interactions on target lipase:
'''bash
python lipase_pocket.py -- input_cif_dir input_pdb_dir alphafold_results/path/to/cif/directory/ interface/path/to/pdb/directory -- output_png plots/path/to/.png
'''

'''bash
python map_consensus_residues.py -- input_target_cif input_ANGPTL_cif input_pdb_dir alphafold_results/path/to/target/.zip alphafold_results/path/to/ANGPTL/.zip interface/path/to/interaction/.pdb -- output_pse plots/path/to/.pse
'''


##Dependencies
- PyMOL (with API access)
- Biopython
- matplotlib
- pandas
- seaborn
- tqdm (for progress bars)
- Python >= 3.8

- install via conda
    - conda install -c conda-forge biopython pymol matplotlib pandas seaborn tqdm
