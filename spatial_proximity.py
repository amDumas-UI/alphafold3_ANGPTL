"""
spatial_proximity.py

Author: Alexander Dumas  
Last Edited: 2025-07-09  
Condition: Stable  

Purpose:
    This script evaluates interface consistency across AlphaFold models predicted at various distance cutoffs.
    For each cutoff folder (e.g., 1.0_cutoff, 2.0_cutoff), it computes the average Jaccard similarity between
    all model pairs to quantify prediction consistency. It outputs a plot to help select the most stable cutoff.

Inputs:
    - BASE_DIR: Path to parent directory containing subfolders named by cutoff (e.g., "1.0_cutoff", "2.0_cutoff")
    - Each folder must contain .pdb files of interface residues.

Outputs:
    - A plot of average Jaccard similarity vs. cutoff value.

Dependencies:
    - Biopython
    - matplotlib
    - seaborn

Usage:
    $ python spatial_proximity.py

Notes:
	This 
"""

# ---- CONFIGURATION ----

import os
import glob
import matplotlib.pyplot as plt
from itertools import combinations
from Bio.PDB import PDBParser
from collections import defaultdict, Counter
import seaborn as sns
import pandas as pd

BASE_DIR = "./interfaces/heterotrimer/A3_A8_LPL/"  # FOLDER CONTAINING SUBFOLDERS FOR EACH CUTOFF. ASSUMES FILES ARE NAMED X_CUTOFF 
PLOT_DIR = os.path.join(os.getcwd(), "plots/heterotrimer/A3_A8_LPL/") # FOLDER TO SAVE PLOTS

def extract_residues(pdb_path):
	"""
	Extracts interface residue identifiers from a PDB file using normalized residue label (resnum + resname).

	Args:
		pdb_path (str): Path to the .pdb file.

	Returns:
		set: Residue labels in the form '123_ALA'
	"""
	parser = PDBParser(QUIET=True)
	structure = parser.get_structure("model", pdb_path)

	residues = set()
	for model in structure:
		for chain in model:
			for res in chain:
				if res.id[0] != " ":
					continue
				res_id = res.id[1]
				resname = res.get_resname()
				label = f"{res_id}_{resname}"
				residues.add(label)
	return residues

# ---- JACCARD SIMILARITY ----

def jaccard(set1, set2):
	"""
	Computes the Jaccard similarity between two sets.

	Args:
		set1 (set): First set.
		set2 (set): Second set.

	Returns:
		float: Jaccard index ranging from 0 (no overlap) to 1 (identical).
	"""
	return len(set1 & set2) / len(set1 | set2) if set1 | set2 else 1.0

def avg_jaccard(interface_sets):
	pairs = combinations(interface_sets, 2)
	scores = [jaccard(a,b) for a, b in pairs]
	return sum(scores) / len(scores) if scores else 1.0


# ---- MAIN ----
def main():
	print("Starting Jaccard similarity analysis...")
	os.makedirs(PLOT_DIR, exist_ok=True)
	cutoff_dirs = sorted([d for d in os.listdir(BASE_DIR) if os.path.isdir(os.path.join(BASE_DIR, d))])
	
	cutoff_data = []
	occurrence_data = {}
	
	for cutoff_dir in cutoff_dirs:
		cutoff_path = os.path.join(BASE_DIR, cutoff_dir)
		print(f"\nProcessing folder: {cutoff_path}")
		pdb_files = sorted(glob.glob(os.path.join(cutoff_path, "*.pdb")))
		print(f"   Found {len(pdb_files)} PDB files")
		
		interfaces = [extract_residues(pdb) for pdb in pdb_files]
		score = avg_jaccard(interfaces)

		try:
			cutoff_value = float(cutoff_dir.split("_")[0]) # ASSUMES FILES ARE NAMED X_CUTOFF
		except ValueError:
			cutoff_value = cutoff_dir
	
		print(f"   Cutoff: {cutoff_value}, Avg Jaccard Similarity: {score:.3f}")

		cutoff_data.append((cutoff_value, score))

		# COUNT RESIDUE OCCURENCES
		all_res = []
		for interface in interfaces:
			all_res.extend(interface)
		counter = Counter(all_res)
		occurrence_data[cutoff_value] = counter
		
	# SORT BY NUMERIC CUTOFF VALUE IF POSSIBLE
	try:
		cutoff_data.sort(key=lambda x: float(x[0]))
	except:
		cutoff_data.sort(key=lambda x: str(x[0]))

	cutoff_values, avg_scores = zip(*cutoff_data)	
	
	# PLOT JACCARD SIMILARITY
	print("\nGenerating plot...")	
	plt.figure(figsize=(8,5))
	plt.plot(cutoff_values, avg_scores, marker='o', label='Average Jaccard')

	# HIGHLIGHT MOST STABLE CUTOFF
	max_idx = avg_scores.index(max(avg_scores))
	best_cutoff = cutoff_values[max_idx]
	best_score = avg_scores[max_idx]
	plt.scatter(best_cutoff, best_score, color='red', s=100, zorder=5, label='Most Stable Cutoff')
	plt.annotate(f"{best_cutoff} Å^2\nScore: {best_score:.3f}",
		xy = (best_cutoff, best_score),
		xytext =(best_cutoff, best_score + 0.05),
		ha = 'center', color='red')
	
	plt.xlabel("Cutoff Distance (Å^2)")
	plt.ylabel("Average Jaccard Similarity")
	plt.grid(True)
	plt.legend()
	plt.tight_layout()
	plot_path = os.path.join(PLOT_DIR, "jaccard_vs_cutoff.png")
	plt.savefig(plot_path)
	plt.close()
	print(f"Plot saved to: {plot_path}")

	
if  __name__ == "__main__":
	main()
	
