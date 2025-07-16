"""
full_trimer_interaction.py

Author: Alexander Dumas  
Date: 2025-07-01

Purpose:
    This script analyzes AlphaFold-predicted protein complexes to identify and visualize residues from a heterotrimer
    (chains A, B, and C) that interact with a lipase (chain D). It extracts .cif files from AlphaFold output zip archives,
    scans PDBs for proximity-based interactions using a 5.0 Å distance threshold, and generates residue-level bar plots.

Inputs:
    - A parent directory containing one folder per lipase.
    - Each lipase folder must contain:
        - AlphaFold .zip files with .cif models for residue reference.
        - .pdb files representing interface-specific models to scan for chain D interactions.

Outputs:
    - Per-lipase bar plots of interacting residues from chain D.
    - A combined plot comparing residue occurrence across all lipases.

Dependencies:
    - Biopython (Bio.PDB)
    - matplotlib, seaborn, pandas
    - tqdm (for progress visualization)

Run:
    $ python full_trimer_interaction.py

Notes:
    - This script is tailored to systems with fixed chain conventions:
        - Chains A–C: Heterotrimer
        - Chain D: Lipase
    - For reversed interaction direction (e.g., chain D interacting with chains A/B/C), see the alternative script.
	- lipase_pocket.py

"""

import os, zipfile, tempfile
import glob
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, MMCIFParser
from collections import Counter
import seaborn as sns
import pandas as pd
import re
import math

# ---- CONFIGURATION ----

PDB_DIR = "./interfaces/heterotrimer/A4_A6/A4x1_A6x2/"  # PARENT FOLDER CONTAINING LIPASE FOLDERS WITH .PDBs INSIDE 
CIF_DIR = "./alphafold_results/heterotrimer/A4_A6/A4x1_A6x2/" # PARENT FOLDER CONTAINING LIPASE FOLDERS WITH .ZIP containing .CIF MODELS INSIDE (FOR REFERENCE RESIDUES)
PLOT_DIR = os.path.join(os.getcwd(), "plots/heterotrimer/A4_A6/A4x1_A6x2/") # FOLDER TO SAVE PLOTS
HETEROTRIMER_CHAINS = {"A","B","C"}
LIPASE_CHAIN = "D"
DEBUG_EXTRACTION = True
EXPECTED_CHAINS = ["A", "B", "C", "D"]
palette = sns.color_palette("Set2", n_colors=len(EXPECTED_CHAINS))
chain_color_map = {ch: palette[i] for i, ch in enumerate(EXPECTED_CHAINS)}
SHARED_COLOR = "#d62728" #RED FROM MATPLOTLIB'S DEFAULT PALETTE
CUTOFF_DISTANCE = 5.0 # Å DISTANCE FOR INTERACTION CUTOFF

def get_reference_residues(zip_path, chain_ids=["A","B","C"]):
	"""
	Extracts all residues from specified chains in the full model .cif file.
	
	Args:
		zip_path (str): Path to .cif models.
		chain_ids (list of str): Which chains to collect reference residues for. 

	Returns:
		residue_labels (list of str): list of reside labels like 'A:ALA_23'.
	"""
	
	parser = MMCIFParser(QUIET = True)
	
	# LOOP THROUGH EACH ZIP ARCHIVE IN THE LIPASE FOLDER
	with zipfile.ZipFile(zip_path, 'r') as zip_ref:
		cif_files = [f for f in zip_ref.namelist() if f.endswith(".cif")]
		if not cif_files:
			print(f"[WARNING] No .cif found in: {zip_path}")
		# EXTRACT EACH .CIF FILE FOR REFERENCE RESIDUE COLLECTION
		with zip_ref.open(cif_files[0]) as cif_handle:
			with tempfile.NamedTemporaryFile(delete=False, suffix=".cif") as tmp_cif:
				tmp_cif.write(cif_handle.read())
				tmp_cif.flush()
				structure = parser.get_structure("ref", tmp_cif.name)
	
	residue_labels = []

	for model in structure:
		for chain in model:
			if chain.id not in chain_ids:
				continue
			for res in chain:
				if res.id[0] != " ":
					continue
				resname = res.get_resname()
				resid = res.id[1]
				label = f"{chain.id}:{resname}_{resid}"
				residue_labels.append(label)
		break # ONLY USE FIRST MODEL

	return residue_labels



def extract_residues(pdb_path):
	"""
	Extracts residues from heterotrimer chains that are within 5.0 Å of the lipase chain

	Args:
		pdb_path (str): Path to the PDB files

	Returns:
		set: A set of strings, where each string represents an interacting residue in the format "Chain:ResName_ResID"

	"""
	parser = PDBParser(QUIET=True)
	structure = parser.get_structure("model", pdb_path)
	residues = set()

	chain = {chain.id: chain for model in structure for chain in model}
	
	if LIPASE_CHAIN not in chain:
		print(f"[WARNING]: Lipase chain {LIPASE_CHAIN} not found in {pdb_path}")
		return residues

	if DEBUG_EXTRACTION:
		print(f"[DEBUG]: Processing PDB: {os.path.basename(pdb_path)}")
		print(f"[DEBUG]: Available chains in PDB: {list(chain.keys())}")
	
	lipase_chain = chain[LIPASE_CHAIN]
	lipase_atoms = [atom for res in lipase_chain for atom in res.get_atoms()]

	if not lipase_atoms and DEBUG_EXTRACTION:
		print(f"[DEBUG]: Lipase chain '{lipase_chain_id}' in {os.path.basename(pdb_path)}")
		return residues
	
	found_interactions_in_file = False
	for chain_id in HETEROTRIMER_CHAINS:
		if chain_id not in chain:
			if DEBUG_EXTRACTION:
				print(f"[DEBUG]: Heterotrimer chain '{chain_id}' not found in {os.path.basename(pdb_path)}. Skipping.")
			continue
		
		heterotrimer_chain_obj = chain[chain_id]
		if not heterotrimer_chain_obj.get_list() and DEBUG_EXTRACTION:
			print(f"[DEBUG]: Heterotrimer chain '{chain_id}' has no residues. Skipping.")
			continue
	
		for res in heterotrimer_chain_obj:
			res_label_full = f"{chain_id}:{res.get_resname()}_{res.get_id()[1]}"
			interaction_found_for_residue = False
			for atom in res.get_atoms():
				for lip_atom in lipase_atoms:
					try:
						distance = atom - lip_atom
						if distance < CUTOFF_DISTANCE:  # DISTANCE CUTOFF IN Å
							# COMMENTED OUT AS THIS SHOWS EACH RESIDUE INTERACTION FOUND
							#if DEBUG_EXTRACTION:
							#	print(f"[DEBUG]: Interaction found: {res_label_full} (atom {atom.get_id()}) <-> Lipase (atom {lip_atom.get_id()}) Distance: {distance: .2f} Å")
							resname = res.get_resname()
							resid = res.get_id()[1]
							label = f"{chain_id}:{resname}_{resid}"
							residues.add(label)
							found_interactions_in_file = True
							interaction_found_for_residue = True
							break # FOUND INTERACTING ATOM -> MOVE TO NEXT
					
					except Exception as e:
						if DEBUG_EXTRACTION:
							print(f"[DEBUG]: Error calculating distance for {res_label_full} in {os.path.basename(pdb_path)}: {e}")
				if interaction_found_for_residue:
					break
		
	if DEBUG_EXTRACTION:
		if found_interactions_in_file:
			print(f"[DEBUG]: Total interacting residues found in {os.path.basename(pdb_path)}: {len(residues)}")
		else:
			print(f"[DEBUG]: No interactions found in {os.path.basename(pdb_path)} within 5.0 Å cutoff.")
	
	if DEBUG_EXTRACTION and residues:
		involved_chains = sorted({r.split(":")[0] for r in residues})
		print(f"[DEBUG]: Residues found from chains: {involved_chains}")
	return residues

def extract_sort_key(res):

	"""
	Extracts sorting keys from a residue label
       

        Args:
                res (str): The residue label string

        Returns:
		tuple: A tuple containing the chain ID and residue number for sorting
	"""
	# EXPECTS FORMAT "A:ALA_34"
	match = re.match(r"([A-Z]):([A-Z]{3})_(\d+)", res)
	if match:
		chain, resname, num = match.groups()
		return (chain, int(num))
	else:
		return ("2", float('inf'))

def process_lipase_data(lipase_path, zip_dir, lipase_name):
	"""
	Processes PDB files for a single lipase, generates a bar plot of interacting residues.

	Args:
		Lipase_path (str): Path to the directory containing .pdb files representing models for the lipase.
		zip_dir (str): Path to the directory containing AlphaFold result .zip archives with .cif models for reference residue collection
		lipase_name (str): Name identifier for the current lipase (usd in output labeling).

	Returns:
		tuple:
			lipase_name (str): The name of the liapse processed.
			plot_df (pd.DataFrame): DataFrame containing residue labels, interaction counts, and chain identifiers.
			total_models (int): Total number of .pdb models analyzed for the lipase
	"""
	print(f"\n--- Processsing data for lipase: {lipase_name} ---")
	pdb_files = sorted(glob.glob(os.path.join(lipase_path, "*.pdb")))
	zip_files = [f for f in os.listdir(zip_dir) if f.endswith(".zip")]

	total_models = len(pdb_files)
	print(f"Found {len(pdb_files)} PDB files in {PDB_DIR}")

	if not pdb_files:
		print(f"[WARNING] No .pdb files found for {lipase_name}. Skipping plot generation for this lipase.")
		return
	if not zip_files:
		print(f"[WARNING] No .zip files found for {lipase_name}. Skipping plot generation for this lipase.")
		return

	# COLLECT REFERENCE RESIDUE SET FROM CHAINS A, B, AND C FROM THE FIRST ZIP FOR EACH LIPASE
	zip_path = os.path.join(zip_dir, zip_files[0])
	reference_residues = get_reference_residues(zip_path)

	if not reference_residues:
		print(f"[ERROR] No trimer chain reference residues found for {lipase_name}")
		return None


	# COUNT RESIDUE OCCURANCE ACROSS MODELS
	all_residues = set() 
	residue_matrix = {}
	
	for pdb_file in pdb_files:
		model_name = os.path.splitext(os.path.basename(pdb_file))[0]
		residues = extract_residues(pdb_file)
		all_residues.update(residues)
		residue_matrix[model_name] = residues
	
		if DEBUG_EXTRACTION:
			print(f"[DEBUG]: Stored for model '{model_name}' in matrix (sample repr): {[repr(r) for r in list(residue_matrix[model_name])[:5]]}")
	
	# DEBUGGING PRINT: BEFORE SORTING ALL_RESIDUES
	print("--- All unique Residues (Before Sort) (Debug) ---")
	for r_unsorted in list(all_residues)[:20]:
		print(f"Unsorted Residue: {r_unsorted}, Raw Key: {extract_sort_key(r_unsorted)}")
	if len(all_residues) > 20: print("...(showing firt 20 out of", len(all_residues), " total)")
	print("---------------------------")

	sorted_residues = sorted(all_residues, key=extract_sort_key)

	print("--- All unique Residues (After Sort) (Debug) ---")
	for r_sorted in list(sorted_residues)[:20]:
		print(f"Sorted Residue: {r_sorted}, Raw Key: {extract_sort_key(r_sorted)}")
	if len(sorted_residues) > 20: print("...(showing firt 20 out of", len(sorted_residues), " total)")
	print("---------------------------")
	
	df = pd.DataFrame(index=reference_residues)

	for model, residues in residue_matrix.items():
		df[model] = [1 if res in residues else 0 for res in reference_residues]
	
	if DEBUG_EXTRACTION:
		print(f"\n[DEBUG]: DataFrame 'df' head for {lipase_name} after population:")
		print(df.head().to_string())
		print(f"\n[DEBUG]: DataFrame 'df' column sums for {lipase_name} after population:")
		print(df.sum(axis=0).to_string())
		print(f"\n[DEBUG]: Using cutoff: {CUTOFF_DISTANCE} Å")
		print(f"\n[DEBUG]: -> {len(residues)} interacting residues")

	#BUILD PLOT_DF
	plot_df = df.sum(axis=1).reset_index()
	plot_df.columns = ["Residue","Count"]
	plot_df["Chain"] = plot_df["Residue"].apply(lambda r: r.split(":")[0])
	plot_df["Percent"] = plot_df["Count"] / total_models * 100

	 # DEBUGGING PRINT TO SHOW THE FINAL ORDER FOR THE PLOT
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
	Main function to iterate through different lipase folders and process data for each.
	"""

	print(f"Starting batch processing for lipases in '{PDB_DIR}'...")
	os.makedirs(PLOT_DIR, exist_ok=True)

	lipase_folders = sorted([os.path.join(PDB_DIR, d) for d in os.listdir(PDB_DIR) if os.path.isdir(os.path.join(PDB_DIR, d))])
	
	all_data = []
	for folder in lipase_folders:
		lipase_name = os.path.basename(folder)
		zip_folder = os.path.join(CIF_DIR, lipase_name)
		lipase_name, plot_df, total_models = process_lipase_data(folder, zip_folder, lipase_name)
		if plot_df is not None:
			all_data.append((lipase_name, plot_df, total_models))
	
		if not all_data:
			print("No data found.")
			return

	#SETUP PLOT GRID
	cols = 1
	rows = len(all_data)
	fig_height = rows * 4
	fig_width = max(15, max(entry[1].shape[0] for entry in all_data) * 0.05)

	fig, axes = plt.subplots(rows, cols, figsize=(fig_width, fig_height), squeeze = False)
	axes = axes.flatten()

	for i, (lipase_name, plot_df, total_models) in enumerate(all_data):
		ax = axes[i]
		x = list(range(len(plot_df)))
		
		# EXTRACT COMMON RESIDUES
		res_keys = plot_df["Residue"].apply(lambda r: r.split(":")[1])
		plot_df["ResKey"] = res_keys
		res_key_counts = res_keys.value_counts()
		shared_keys = set(res_key_counts[res_key_counts > 1].index) # ONLY THOSE APPEARING IN MULTIPLE CHAINS
		
		model_threshold = int(total_models * 0.33) # HAVE TO APPEAR IN > 0.33 OF THE TOTAL MODELS
		valid_keys = plot_df[plot_df["ResKey"].isin(shared_keys)]
		valid_keys = valid_keys[valid_keys["Count"] >= model_threshold]
		final_highlight_keys = set(valid_keys["ResKey"])

		highlight_flags = plot_df["ResKey"].isin(final_highlight_keys)
		
		bar_colors = [
			SHARED_COLOR if is_shared else chain_color_map[plot_df["Chain"].iloc[i]]
			for i, is_shared in enumerate(highlight_flags)
		]

		highlight_edge_colors = ["black" if flag else "none" for flag in highlight_flags]
	
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
		ax.set_title(f"{lipase_name} Interface Residues", fontsize=25)	
		ax.set_ylabel("Percent of Models Interacting", fontsize=20)
		ax.set_ylim(0,100)
		ax.tick_params(axis='y', labelsize=15)
		ax.set_xlabel("Residue")
	
	handles = [plt.Rectangle((0,0),1,1, color=chain_color_map[ch]) for ch in EXPECTED_CHAINS]
	fig.legend(handles, [f"Chain {ch}" for ch in EXPECTED_CHAINS], title="Chain", loc='upper right')
	plt.tight_layout(rect=[0, 0, 0.95, 0.98])
	highlight_patch = plt.Rectangle((0,0),1,1, facecolor='white', edgecolor='black', linewidth=1.2, label='Shared across chains')
	handles.append(highlight_patch)
	
	filename = f"combined_residue_occurrence.png"
	out_path = os.path.join(PLOT_DIR, filename)
	plt.savefig(out_path, dpi=300)
	plt.close()
	print(f"Saved combined plot to: {out_path}")
	print("Batch processing complete")
	
if  __name__ == "__main__":
	main()
	
