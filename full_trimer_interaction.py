"""
full_trimer_interaction.py

Author: Alexander Dumas
Date: 2025-07-09

Purpose:
    This script processes AlphaFold-predicted models of heterotrimeric ANGPTL proteins 
    and identifies interface residues from chains A and B that interact with chain C. 
    The goal is to quantify and visualize how consistently the trimer interface is formed 
    across multiple prediction models.

Inputs:
    - Zip archives from AlphaFold structure predictions (.cif files inside)
    - Each subdirectory in BASE_DIR should contain AlphaFold .zip files per lipase

Outputs:
    - Per-lipase bar plots of interface residue frequency (A/B interacting with C)
    - Combined bar plot across lipases highlighting conserved residues

Dependencies:
    - Biopython
    - matplotlib, pandas, seaborn
    - tqdm

Run:
    python full_trimer_interaction.py

Notes:
    Residues are labeled as "Chain:RESNAME_RESNUM" (e.g., A:ARG_24) and sorted by chain and index.
    Shared residues are highlighted in red if they appear in ≥1/3 of models and in all lipases.
"""

# ---- CONFIGURATION ----

import os, zipfile, tempfile, re
import glob
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser
from collections import Counter
from Bio.PDB import MMCIFParser
import seaborn as sns
import pandas as pd
import math
from tqdm import tqdm


BASE_DIR = "./alphafold_results/homotrimer/A4/"  # PARENT FOLDER CONTAINING LIPASE FOLDERS WITH .ZIPs INSIDE 
PLOT_DIR = os.path.join(os.getcwd(), "plots/homotrimer/A4/") # FOLDER TO SAVE PLOTS
HETEROTRIMER_CHAINS = {"A","B","C"}
LIPASE_CHAIN = "D"
DEBUG_EXTRACTION = True
EXPECTED_CHAINS = ["A", "B", "C", "D"]
palette = sns.color_palette("Set2", n_colors=len(EXPECTED_CHAINS)) # COLOR PALETTE FOR BARPLOT
chain_color_map = {ch: palette[i] for i, ch in enumerate(EXPECTED_CHAINS)}
SHARED_COLOR = "#d62728" #RED FROM MATPLOTLIB'S DEFAULT PALETTE
TRIMER_THRESHOLD = 3 # RATIO OF HOW MANY PREDICTIONS NEED TO CONTAIN THE SAME RESIDUE. EX. 3:  IF >= 1/3 OF TOTAL PREDICTION MODELS
CUTOFF_DISTANCE = 5 # Å DISTANCE CUTOFF FOR INTERACTIONS 
chain_AB = ["A","B"]
chain_C = "C"
		

def extract_residues(zip_path):
	"""
	Extracts residues from chains A and B that are within a 5 Å distance of chain C across all .cif models within a provided AlphaFold zip archive. 

	Args:
		zip_path (str): Path to the .zip file containing AlphaFold .cif model predictions.

	Returns:
		tuple:
			- dict (str, set[str]): A dictionary mapping each model name to a set of interacting residue labels in the format 'Chain:RESNAME_RESNUM'.
			- list (str): A sorted list of all residues from chains A and B (reference residues) used to ensure consistent x-axis plotting across models.
	"""
	parser = MMCIFParser(QUIET=True)
	residues_per_model = {}
	ref_labels = []
	residues = set()
	chain_AB = ["A","B"]
	chain_C = "C"

#	zip_files = [f for f in os.listdir(lipase_folder) if f.endswith(".zip")]
#	for zip_file in zip_files:
#		zip_path = os.path.join(lipase_folder, zip_file)
#		if DEBUG_EXTRACTION:
#			print(f"\n[INFO] Processing zip: {zip_file}")
#
#		# CLEAN TEMP FOLDER BEFORE EXTRACTION
#		if os.path.exists(tmp_dir):
#			shutil.rmtree(tmp_dir)
#		os.makedirs(tmp_dir)

	print(f"\n[INFO] Extracting residues from ZIP: {os.path.basename(zip_path)}")

	with tempfile.TemporaryDirectory() as tmpdir:		
		with zipfile.ZipFile(zip_path, 'r') as zip_ref:
			zip_ref.extractall(tmpdir)

		print(f"[INFO] Extracted files: {len(os.listdir(tmpdir))}")
		#for fname in os.listdir(tmpdir):
		cif_files = [f for f in os.listdir(tmpdir) if f.endswith(".cif")]
		print(f"[INFO] Extracted {len(cif_files)} .cif files. Starting processing...")

		for fname in cif_files:
			if not fname.endswith(".cif"):
				print(f"[SKIP] Not a .cif file: {fname}")
				continue
			model_path = os.path.join(tmpdir, fname)
			model_name = os.path.splitext(fname)[0]
			structure = parser.get_structure(model_name, model_path)
			print(f"[PROCESS] Reading model: {model_name}")

			try:
				structure = parser.get_structure(model_name, model_path)
				chain = {chain.id: chain for model in structure for chain in model}
				print(f" - Chains found: {list(chain.keys())}")
			except Exception as e:
				print(f"[ERROR] Could not parse structure {fname}: {e}")
				continue

			if chain_C not in chain:
				print(f"[ERROR] Chain C not found in {zip_path}")
				continue

			c_atoms = [atom for res in chain[chain_C] for atom in res.get_atoms()]	
			model_residues = set()

			for chain_id in chain_AB:
				if chain_id not in chain:
					print(f"[WARNING] Chain {chain_id} not found in {fname} - skipping this chain.")
				for res in chain[chain_id]:
					resname = res.get_resname()
					resid = res.get_id()[1]
					label = f"{chain_id}:{resname}_{resid}"
					residues.add(label)

					for atom in res.get_atoms():
						if any((atom - c_atom) < CUTOFF_DISTANCE for c_atom in c_atoms):  # DISTANCE CUTOFF IN Å
							#if DEBUG_EXTRACTION:
							#	print(f"DEBUG: Interaction found: {res_label_full} (atom {atom.get_id()}) <-> Lipase (atom {lip_atom.get_id()}) Distance: {distance: .2f} Å")
							model_residues.add(label)
							break # FOUND INTERACTING ATOM -> MOVE TO NEXT
			
			print(f" - Total interacting residues in {model_name}: {len(model_residues)}")
			residues_per_model[model_name] = model_residues
	
	sorted_residue_list = sorted(residues, key=extract_sort_key)
	print(f"[DONE] Finished extracting from {os.path.basename(zip_path)}. Models parsed: {len(residues_per_model)}")			
	return residues_per_model, sorted_residue_list			

def extract_sort_key(res):

	"""
	Generates a sorting key from a residue label of the form 'Chain:RESNAME_RESNUM'. 

        Args:
                res (str): The residue label string (e.g., "A:ALA_34").

        Returns:
		tuple: (chain_id, residue_number) used to sort residues for plotting.
	"""

	# EXPECTS FORMAT "A:ALA_34"
	match = re.match(r"([A-Z]):([A-Z]{3})_(\d+)", res)
	if match:
		chain, resname, num = match.groups()
		return (chain, int(num))
	else:
		return ("2", float('inf'))

def process_lipase_data(lipase_path, lipase_name):
	"""
	Processes AlphaFold .zip model files for a given lipase to identify interface residues between chains A/B and chain C of the heterotrimer.
		
	Args:
		Lipase_path (str): Path to the folder containing PDB files for this lipase
		lipase_name (str): The name of the lipase

	Returns:
		tuple:
			- lipase_name (str): The name of the lipase for tracking and plotting.
			- plot_df (pd.DataFrame): DataFrame containing 'Residue', 'Count', 'Chain', 'Percent' columns.
			- total_models (int): Total number of Alphafold models analyzed from .zip files.
	"""
	print(f"\n--- Processsing data for lipase: {lipase_name} ---")
	zip_files = [f for f in os.listdir(lipase_path) if f.endswith(".zip")]
	if not zip_files:
		print(f"WARNING: No zip files found in {lipase_path}")
		return None
	
	# COUNT RESIDUE OCCURRENCE
	residue_matrix = {}	
	total_models = 0

	print(f"[INFO] Found {len(zip_files)} zip files in {lipase_path}")
	for zip_file in tqdm(zip_files, desc=f"Processing zips for {lipase_name}", unit="zip"):
		zip_path = os.path.join(lipase_path, zip_file)
		model_res_map, ref_residues = extract_residues(zip_path)
		for model, residues in model_res_map.items():
			residue_matrix[model] = residues
			total_models += 1

	if not residue_matrix:
		print(f"No residues found for {lipase_name}.")
		return None
		
	# DEBUGGING PRINT: BEFORE SORTING ALL_RESIDUES
	print("--- All unique Residues (Before Sort) (Debug) ---")
	for r_unsorted in list(ref_residues)[:20]:
		print(f"Unsorted Residue: {r_unsorted}, Raw Key: {extract_sort_key(r_unsorted)}")
	if len(ref_residues) > 20: print("...(showing firt 20 out of", len(ref_residues), " total)")
	print("---------------------------")

	sorted_residues = sorted(ref_residues, key=extract_sort_key)

	print("--- All unique Residues (After Sort) (Debug) ---")
	for r_sorted in list(sorted_residues)[:20]:
		print(f"Sorted Residue: {r_sorted}, Raw Key: {extract_sort_key(r_sorted)}")
	if len(sorted_residues) > 20: print("...(showing firt 20 out of", len(sorted_residues), " total)")
	print("---------------------------")
	
	df = pd.DataFrame(index=ref_residues)

	for model, residues in residue_matrix.items():
		df[model] = [1 if res in residues else 0 for res in ref_residues]
	
	if DEBUG_EXTRACTION:
		print(f"\nDEBUG: DataFrame 'df' head for {lipase_name} after population:")
		print(df.head().to_string())
		print(f"\nDEBUG: DataFrame 'df' column sums for {lipase_name} after population:")
		print(df.sum(axis=0).to_string())

	#BUILD PLOT_DF
	plot_df = df.sum(axis=1).reset_index()
	plot_df.columns = ["Residue","Count"]
	plot_df["Chain"] = plot_df["Residue"].apply(lambda r: r.split(":")[0])
	plot_df["Percent"] = plot_df["Count"] / total_models * 100
	
	 # Debugging print to show the final ordered residues for the bar plot
	print("--- Final Ordered Residues for Barplot (Debug) ---")
	print(plot_df.head(20))
	if len(plot_df) > 20:
		print("...(showing first 20 out of", len(plot_df), "total)")
	print("plot_df shape:", plot_df.shape)
	print(plot_df["Residue"].tolist()[-10:])
	print("-------------------------------------------")
		
	print(f"Generating residue occurence bar plot for {lipase_name}...")

	return lipase_name, plot_df, total_models
	

def main():
	"""
	Orchestrates batch processing of multiple lipase folders containing AlphaFold predictions. 
	Generates bar plots showing residues from chains A and B that interact with chain C across models.
	Highlights residues that are consistently shared across lipases and above frequency thresholds.

	Outputs:
		- Individual bar plots for each lipase folder.
		- Acombinded plot highlighting shared interface residues

	"""

	print(f"Starting batch processing for lipases in '{BASE_DIR}'...")
	os.makedirs(PLOT_DIR, exist_ok=True)

	lipase_folders = sorted([os.path.join(BASE_DIR, d) for d in os.listdir(BASE_DIR) if os.path.isdir(os.path.join(BASE_DIR, d))])
	
	all_data = []
	for folder in lipase_folders:
		lipase_name = os.path.basename(folder)
		result = process_lipase_data(folder, lipase_name)
		if result:
			all_data.append(result)
	
		if not all_data:
			print("No data collected.")
			return

	# COLLECT RESIDUES FOR TRACKING ACROSS LIPASES: MUST APPEAR IN ALL LIPASES AND IN OVER HALF THE PREDICTIONS
	global_residue_counts = Counter()
	global_model_res_map = {}
	lipase_presence = Counter()
	total_models = 0
	
	for lipase_name, plot_df, n_models in all_data:
		total_models += n_models
		residues_seen = set()
		for _, row in plot_df.iterrows():
			if row["Count"] > 0:
				residues_seen.add(row["Residue"])
				global_residue_counts[row["Residue"]] += row["Count"]
		for res in residues_seen:
			lipase_presence[res] += 1
		global_model_res_map[lipase_name] = residues_seen
	
	# FIND RESIDUE PRESENT IN 5/5 LIPASES
	#shared_across_lipase = set.intersection(*global_model_res_map.values())
	required_lipases = 5

	# APPLY 1/3 MODEL CUTOFF
	highlight_residues = {
		residue for residue in global_residue_counts
		if lipase_presence[residue] >= required_lipases
		and global_residue_counts[residue] >= total_models / TRIMER_THRESHOLD
	}

	
	#SETUP PLOT GRID
	cols = 1
	rows = len(all_data)
	fig_height = rows * 4
	fig_width = max(15, max(entry[1].shape[0] for entry in all_data) * 0.05)

	fig, axes = plt.subplots(rows, cols, figsize=(fig_width, fig_height), squeeze = False)
	axes = axes.flatten()

	for i, (lipase_name, plot_df, _) in enumerate(all_data):
		ax = axes[i]
		x = list(range(len(plot_df)))
		
		# EXTRACT COMMON RESIDUES
		res_keys = plot_df["Residue"].apply(lambda r: r.split(":")[1])
		plot_df["ResKey"] = res_keys
		res_key_counts = res_keys.value_counts()
		shared_keys = set(res_key_counts[res_key_counts > 1].index) # ONLY THOSE APPEARING IN MULTIPLE CHAINS

		highlight_flags = plot_df["ResKey"].isin({
			res.split(":")[1] for res in highlight_residues
		}) 
		
		bar_colors = [
			SHARED_COLOR if is_highlight else chain_color_map[row["Chain"]]
			for is_highlight, (_, row) in zip(highlight_flags, plot_df.iterrows()) 
		]

		for i in range(len(plot_df)):
			ax.bar(
				x[i],
				plot_df["Percent"].iloc[i],
				color=bar_colors[i],
				width=0.6
			)
		label_step = max(1, len(plot_df) // 50)
		plot_df["ResLabel"] = plot_df["Residue"].apply(lambda r: r.replace("_", ""))
		ax.set_xticks(x)
		ax.set_xticklabels(plot_df["ResLabel"], rotation=90, ha='right', fontsize=3)
		for tick, is_shared in zip(ax.get_xticklabels(), highlight_flags):
			tick.set_color(SHARED_COLOR if is_shared else "black")
		ax.set_title(f"{lipase_name}: A+B vs C Interface Residues", fontsize=30)	
		ax.set_ylabel("Percent of Models Interacting")
		ax.set_xlabel("Residue")
	
	handles = [plt.Rectangle((0,0),1,1, color=chain_color_map[ch]) for ch in EXPECTED_CHAINS]
	fig.legend(handles, [f"Chain {ch}" for ch in EXPECTED_CHAINS], title="Chain", loc='upper right')
	plt.tight_layout(rect=[0, 0, 0.95, 0.98])
	highlight_patch = plt.Rectangle((0,0),1,1, facecolor='white', edgecolor='black', linewidth=1.2, label='Shared across chains')
	handles.append(highlight_patch)
	
	filename = f"full_trimer_combined_residue_occurrence.png"
	out_path = os.path.join(PLOT_DIR, filename)
	plt.savefig(out_path, dpi=300)
	plt.close()
	print(f"Saved combined plot to: {out_path}")
	print("Batch processing complete")
	
if  __name__ == "__main__":
	main()
	
