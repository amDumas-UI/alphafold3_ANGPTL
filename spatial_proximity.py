####
# AUTHOR: ALEXANDER DUMAS
# DATE OF LAST EDIT: 6-17-2025
# CONDITION: IN PROGRESS
# NOTES: This script will identify similarity between interaction .pdb files to determine a good cutoff for InterfaceResidues()
# USING JACCARD SIMILARITY MATRIC: MEASURE HOW SIMILAR TWO SETS ARE. RANGE FROM HIGH SIMILARITY (1) TO LOW SIMILARITY (0)
# IF MODELS AT EACH CUTOFF WILL BE TESTED FOR SEEING HOW SIMILAR THEY ARE. HIGH SIMILARITY AT ONE CUTOFF OR ANOTHER WILL HELP DETERMINE WHICH CUTOFF IS BEST
####

# ---- CONFIGURE ----

import os
import glob
import matplotlib.pyplot as plt
from itertools import combinations
from Bio.PDB import PDBParser
from collections import defaultdict, Counter
import seaborn as sns
import pandas as pd

# ---- PARAMETERS ----

BASE_DIR = "./interfaces/heterotrimer/A3_A8_LPL/"  # FOLDER CONTAINING SUBFOLDERS FOR EACH CUTOFF. ASSUMES FILES ARE NAMED X_CUTOFF 
PLOT_DIR = os.path.join(os.getcwd(), "plots/homotrimer/A4") # FOLDER TO SAVE PLOTS

DIR_A = "./interfaces/homotrimer/A4/h_A4_LPL/"
DIR_B = "./interfaces/homotrimer/A4_LPL_Altered/shuffled/"


# ---- LOAD INTERFACE RESIDUE DATA ----

# NORMALIZING TO RESIDUE NUMBER TO BE ABLE TO JUST COMPARE BASED ON REISDUE NUMBER PRESENCE AND RESIDUE NAME
def extract_residues(pdb_path):
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
	return len(set1 & set2) / len(set1 | set2) if set1 | set2 else 1.0

def avg_jaccard(interface_sets):
	pairs = combinations(interface_sets, 2)
	scores = [jaccard(a,b) for a, b in pairs]
	return sum(scores) / len(scores) if scores else 1.0

altered_files = sorted([f for f in os.listdir(DIR_B) if f.endswith(".pdb")])
unaltered_files = sorted([f for f in os.listdir(DIR_A) if f.endswith(".pdb")])

if len(altered_files) != len(unaltered_files):
	print("[WARNING] Directory lengths differ - truncating to shortest list.")
n = min(len(altered_files), len(unaltered_files))

print(f"{'Altered':<30} | {'Unaltered':<30} | Jaccard | Shared | Total")
print("-" * 80)

for i in range(n):
	f1 = altered_files[i]
	f2 = unaltered_files[i]
	path1 = os.path.join(DIR_B, f1)
	path2 = os.path.join(DIR_A, f2)
		
	r1 = extract_residues(path1)
	r2 = extract_residues(path2)
	score = jaccard(r1, r2)
	print(f"{f1:<30} | {f2:<30} | {score:>7.3} | {len(r1 & r2):>6} | {len(r1 | r2):>5}")
	

#pairwise_scores = []
#for i, (set_A, set_B) in enumerate(zip(interface_A, interface_B)):
#	score = jaccard(set_A, set_B)
#	pairwise_scores.append(score)
#	print(f"Model {i+1}: Jaccard = {score:.3f}")
#
#plt.figure(figsize=(8,5))
#plt.plot(range(1, len(pairwise_scores)+1), pairwise_scores, marker='o', label='Jaccard Similarity')
#plt.xlabel("Model Index")
##plt.ylabel("Jaccard Similarity")
#plt.title("Pairwise Jaccard Similarity: Original vs. Shuffled Chain Order")
#plt.grid(True)
#plt.tight_layout()
#plt_path = os.path.join(PLOT_DIR, "jaccard_per_model.png")
#plt.savefig(plt_path)
#plt.close()
#print(f"Saved plot to {plt_path}")

#avg_score = sum(pairwise_scores) / len(pairwise_scores) if pairwise_scores else 1.0
#print(f"\nAverage Jaccard similarity across models: {avg_score:.3f}")

# ---- MAIN ----
#def main():
#	print("Starting Jaccard similarity analysis...")
#	os.makedirs(PLOT_DIR, exist_ok=True)
#	cutoff_dirs = sorted([d for d in os.listdir(BASE_DIR) if os.path.isdir(os.path.join(BASE_DIR, d))])
#	
#	cutoff_data = []
#	occurrence_data = {}
#	
#	for cutoff_dir in cutoff_dirs:
#		cutoff_path = os.path.join(BASE_DIR, cutoff_dir)
#		print(f"\nProcessing folder: {cutoff_path}")
#		pdb_files = sorted(glob.glob(os.path.join(cutoff_path, "*.pdb")))
#		print(f"   Found {len(pdb_files)} PDB files")
#		
#		interfaces = [extract_residues(pdb) for pdb in pdb_files]
#		score = avg_jaccard(interfaces)
#
#		try:
#			cutoff_value = float(cutoff_dir.split("_")[0]) # ASSUMES FILES ARE NAMED X_CUTOFF
#		except ValueError:
#			cutoff_value = cutoff_dir
#	
#		print(f"   Cutoff: {cutoff_value}, Avg Jaccard Similarity: {score:.3f}")
#
#		cutoff_data.append((cutoff_value, score))
#
#		# COUNT RESIDUE OCCURENCES
#		all_res = []
#		for interface in interfaces:
#			all_res.extend(interface)
#		counter = Counter(all_res)
#		occurrence_data[cutoff_value] = counter
#		
#	# SORT BY NUMERIC CUTOFF VALUE IF POSSIBLE
#	try:
#		cutoff_data.sort(key=lambda x: float(x[0]))
#	except:
#		cutoff_data.sort(key=lambda x: str(x[0]))
#
#	cutoff_values, avg_scores = zip(*cutoff_data)	
#	
#	# PLOT JACCARD SIMILARITY
#	print("\nGenerating plot...")	
#	plt.figure(figsize=(8,5))
#	plt.plot(cutoff_values, avg_scores, marker='o', label='Average Jaccard')
#
#	# HIGHLIGHT MOST STABLE CUTOFF
#	max_idx = avg_scores.index(max(avg_scores))
#	best_cutoff = cutoff_values[max_idx]
#	best_score = avg_scores[max_idx]
#	plt.scatter(best_cutoff, best_score, color='red', s=100, zorder=5, label='Most Stable Cutoff
#	plt.annotate(f"{best_cutoff} Å^2\nScore: {best_score:.3f}",
#		xy = (best_cutoff, best_score),
#		xytext =(best_cutoff, best_score + 0.05),
#		ha = 'center', color='red')
#	
#	plt.xlabel("Cutoff Distance (Å^2)")
#	plt.ylabel("Average Jaccard Similarity")
#	plt.grid(True)
#	plt.legend()
#	plt.tight_layout()
#	plot_path = os.path.join(PLOT_DIR, "jaccard_vs_cutoff.png")
#	plt.savefig(plot_path)
#	plt.close()
#	print(f"Plot saved to: {plot_path}")
#
#	# PLOT RESIDUE OCCURENCE HEATMAP
#	print("Generating per-residue occurrence heatmap...")
#	all_residues = sorted({res for counts in occurrence_data.values() for res in counts})
#	heatmap_data = pd.DataFrame(index=all_residues, columns=cutoff_values)
#
#	for cutoff in cutoff_values:
#		for res in all_residues:
#			heatmap_data.loc[res, cutoff] = occurrence_data[cutoff].get(res, 0)
#
#	heatmap_data = heatmap_data.astype(float)
#	plt.figure(figsize=(12, max(6, len(all_residues) * 0.2)))
#	sns.heatmap(heatmap_data, cmap="viridis", cbar_kws={'label': 'Occurrence Count'})
#	plt.xlabel("ASA Cutoff (Å²)")
#	plt.ylabel("Residue")
#	heatmap_path = os.path.join(PLOT_DIR, "residue_occurrence_heatmap.png")
#	plt.tight_layout()
#	plt.savefig(heatmap_path)
#	plt.close()
#	print(f"Saved heatmap to: {heatmap_path}")
	
#if  __name__ == "__main__":
#	main()
	
