"""
lipase_pocket.py

Author: Alexander Dumas
Date: 2025-09-01

Purpose:
    This script identifies and visualizes lipase (chain D) interface residues across multiple AlphaFold predictions 
    involving heterotrimeric ANGPTL complexes (chains A, B, and C). It extracts residues from chain D that are 
    spatially close to any residue in chains A–C using a distance-based cutoff.

    The script collects all interacting residues from .pdb models, uses a single .cif model per lipase to retrieve 
    reference residue labels from chain D, and generates bar plots of residue-level interaction frequencies.

Inputs:
    - Directory containing folders per lipase, each with multiple .pdb prediction models.
    - Corresponding AlphaFold .zip files for each lipase, used to extract chain D residue references from .cif models.

Outputs:
    - Bar plots per lipase showing interacting residues from chain D across multiple predictions.
    - Consistent residue indexing based on chain D references ensures comparability across models.

Dependencies:
    - Biopython (Bio.PDB)
    - matplotlib, seaborn, pandas
    - tqdm (optional, for progress bars)

Usage:
    $ python lipase_pocket.py

Note:
    - The script assumes consistent chain labeling (chains A–C for heterotrimer, D for lipase).
    - Only residues from chain D are plotted, even though proximity is assessed relative to chains A–C.
"""


import os
import glob
import zipfile, tempfile, os
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, MMCIFParser
from collections import Counter
import seaborn as sns
import pandas as pd
import re
import math

# ---- CONFIGURATION ----

BASE_DIR = "./interfaces/heterotrimer/A8_A5/A8x1_A5x2/"  # PARENT FOLDER CONTAINING LIPASE FOLDERS WITH .PDBs INSIDE
CIF_DIR = "./alphafold_results/heterotrimer/A8_A5/A8x1_A5x2/" # PARENT FOLDER CONTAINING A ZIP FILE FOR .CIF EXTRACTION
PLOT_DIR = os.path.join(os.getcwd(), "plots/heterotrimer/A8_A5/A8x1_A5x2/") # FOLDER TO SAVE PLOTS
HETEROTRIMER_CHAINS = {"A","B","C"}
LIPASE_CHAIN = "D"
DEBUG_EXTRACTION = True
EXPECTED_CHAINS = ["A", "B", "C", "D"]
palette = sns.color_palette("Set2", n_colors=len(EXPECTED_CHAINS))
chain_color_map = {ch: palette[i] for i, ch in enumerate(EXPECTED_CHAINS)}
SHARED_COLOR = "#1f77b4" #Blue FROM MATPLOTLIB'S DEFAULT PALETTE
CUTOFF_DISTANCE = 5.0 # Å DISTANCE FOR INTERACTION CUTOFF
POCKET_RANGES = {
	"LPL": set(range(50, 71)) | set(range(215, 246)), # TAKEN FROM A4 HOMOTRIMER INTERACTION
	"EL": set(range(67, 88)) | set(range(231, 263)), # TAKEN FROM BLAST COMPARED WITH LPL
	"PL": set(range(69, 90)) | set(range(237, 268)),
	"HL": set(range(63, 74)) | set(range(231,262))  
}	

FREQUENCY_THRESHOLD = 0.33

def get_reference_residues(zip_path, chain_ids="D"):
	"""
	Extracts all residues from specified chains in the full model .cif file

	Args:
		zip_path (str): Path to .cif models.
		chain_ids (str): Which chain to collect reference residues from.
		
	Returns:
		residue_labels (list of str): list of residue labels like 'A:ALA_23'.
	"""
	
	parser = MMCIFParser(QUIET=True)
	
	with zipfile.ZipFile(zip_path, 'r') as zip_ref:
		cif_files = [f for f in zip_ref.namelist() if f.endswith(".cif")]
		if not cif_files:
			print(f"[WARNING] No .cif found in: {zip_path}")

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

def extract_residues(pdb_path, target_chain="D", query_chains={"A","B","C"}):
	"""
	Extracts residues from the lipase chain (chain D) that are within 5.0 Å of the any trimer chain (A, B, or C)

	Args:
		pdb_path (str): Path to the PDB files

	Returns:
		set: A set of strings, where each string represents an interacting residue in the format "Chain:ResName_ResID"

	"""
	parser = PDBParser(QUIET=True)
	structure = parser.get_structure("model", pdb_path)

	chain = {chain.id: chain for model in structure for chain in model}
	
	if target_chain not in chain:
		print(f"WARNING: Lipase chain {target_chain} not found in {pdb_path}")
		return set()

	if DEBUG_EXTRACTION:
		print(f"[DEBUG]: Processing PDB: {os.path.basename(pdb_path)}")
		print(f"[DEBUG]: Available chains in PDB: {list(chain.keys())}")
	
	target_atoms = [atom for res in chain[target_chain] for atom in res.get_atoms()]
	interacting = set()
	found_interactions_in_file = False
	
	lip_chain = chain[target_chain]
	lip_residues = lip_chain.get_residues()
	
	query_atoms = []
	for q in query_chains:
		if q in chain:
			for res in chain[q]:
				query_atoms.extend(res.get_atoms())
				
	for res in lip_residues:
		resname = res.get_resname()
		resid = res.id[1]
		label = f"{target_chain}:{resname}_{resid}"
		for atom in res.get_atoms():
			if any(atom - q_atom < CUTOFF_DISTANCE # RETURN TRUE IF WITHIN CUTOFF_DISTANCE AND ANY OF THE BOTTOM LOOPS ARE TRUE
				for q_atom in query_atoms):
				interacting.add(label)	
				found_interactions_in_file = True
				break # FOUND INTERACTING ATOM -> MOVE TO NEXT
	if DEBUG_EXTRACTION:
		if found_interactions_in_file:
			print(f"[DEBUG]: Total interacting residues found in {os.path.basename(pdb_path)}: {len(interacting)}")
		else:
			print(f"[DEBUG]: No interactions found in {os.path.basename(pdb_path)} within 5.0 Å cutoff.")
	
	if DEBUG_EXTRACTION and interacting:
		involved_chains = sorted({r.split(":")[0] for r in interacting})
		print(f"[DEBUG]: Residues found from chains: {involved_chains}")
	return interacting


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

def process_lipase_data(pdb_dir, zip_dir, lipase_name):
	"""
	Processes structural models of a single lipase, identifies residues in chain D (lipase) 
	that interact with chains A, B, or C of a heterotrimer, and constructs a dataframe 
	for plotting interaction frequencies.
	
	Args:
		pdb_dir (str): Path to the directory containing .pdb models for this lipase
		zip_dir (str): Path to the directory containing .zip archives with full AlphaFold server models (.cif)
		lipase_name (str): The lipase identifier (used for labeling plots and outputs)
	
	Returns:
		tuple: (lipase_name, plot_df, total_models)
			- plot_df is a DataFrame with residue names, interaction count, and chain info.
			- total_models is the number of .pdb models analyzed.
	"""
	print(f"\n--- Processsing data for lipase: {lipase_name} ---")
	pdb_files = sorted(glob.glob(os.path.join(pdb_dir, "*.pdb")))
	zip_files = [f for f in os.listdir(zip_dir) if f.endswith(".zip")]

	total_models = len(pdb_files)
	print(f"Found {len(pdb_files)} PDB files in {BASE_DIR}")

	if not pdb_files:
		print(f"[WARNING] No PDB files found for {lipase_name}. Skipping plot generation for this lipase.")
		return
	if not zip_files:
		print(f"[WARNING] No zip files found for {lipase_name}. Skipping plot generation for this lipase.")
		return

	# COLLECT REFERENCE RESIDUE SET FROM CHAIN D IN THE FIRST ZIP 
	zip_path = os.path.join(zip_dir, zip_files[0])
	reference_residues = get_reference_residues(zip_path)

	if not reference_residues:
		print(f"[ERROR] No chain D reference residues found for {lipase_name}")
		return None

	parser = PDBParser(QUIET=True)
	
	# COUNT RESIDUE OCCURANCE ACROSS MODELS
	residue_matrix = {}
	total_models = 0
	
	for pdb_file in pdb_files:
		model_name = os.path.splitext(os.path.basename(pdb_file))[0]
		interacting_residues = extract_residues(pdb_file)
		residue_matrix[model_name] = interacting_residues
		total_models += 1



	# THIS IS BASICALLY IN LINE CODE OF THE EXTRACT RESIDUES FUNCTION BUT IT HAS NOT BEEN CHANGED TO COLLECT CHAIN D RESIDUES		
#		chain_map = {chain.id: chain for model in structure for chain in model}
#
#		interacting_residues = set()	
#
#		if "D" not in chain_map:
#			print(f"[SKIP] Chain D not found in {pdb_file}")
#			continue
#		
#		chain_D_atoms = [atom for res in chain_map["D"] for atom in res.get_atoms()]
#		if not chain_D_atoms:
#			continue
#	
#		for het_chain in ["A","B","C"]:
#			if het_chain not in chain_map:
#				continue
#			for res in chain_map[het_chain]:
#				for atom in res.get_atoms():
#					for lip_atom in chain_D_atoms:
#						try:
#							if atom - li_atom < CUTOFF_DISTANCE:
#								resname = lip_atom.get_parent().get_resname()
#								resid = lip_atom.get_parent().id[1]
#								label = f"D:{resname}_{resid}"
#								interacting_residues.add(label)
#								break # FOUND INTERACTION
#						except:
#							continue
#	
	
		if DEBUG_EXTRACTION:
			print(f"[DEBUG]: Stored for model '{model_name}' in matrix (sample repr): {[repr(r) for r in list(residue_matrix[model_name])[:5]]}")
	
	# DEBUGGING PRINT: BEFORE SORTING ALL_RESIDUES
	print("--- All unique Residues (Before Sort) (Debug) ---")
	for r_unsorted in list(interacting_residues)[:20]:
		print(f"Unsorted Residue: {r_unsorted}, Raw Key: {extract_sort_key(r_unsorted)}")
	if len(interacting_residues) > 20: print("...(showing firt 20 out of", len(interacting_residues), " total)")
	print("---------------------------")

	sorted_residues = sorted(interacting_residues, key=extract_sort_key)

	print("--- All unique Residues (After Sort) (Debug) ---")
	for r_sorted in list(sorted_residues)[:20]:
		print(f"Sorted Residue: {r_sorted}, Raw Key: {extract_sort_key(r_sorted)}")
	if len(sorted_residues) > 20: print("...(showing first 20 out of", len(sorted_residues), " total)")
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
		print(f"\n[DEBUG]: -> {len(interacting_residues)} interacting residues")

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
		
	print(f"Generating lipase interaction bar plot for {lipase_name}...")

	return lipase_name, plot_df, total_models
	

def main():
	"""
	Main function to iterate through different lipase folders and process data for each.
	"""

	print(f"Starting batch processing for lipases in '{BASE_DIR}'...")
	os.makedirs(PLOT_DIR, exist_ok=True)

	lipase_folders = sorted([
		d for d in os.listdir(BASE_DIR) 
		if os.path.isdir(os.path.join(BASE_DIR, d))
	])
	
	all_data = []
	for lipase_name in lipase_folders:
		pdb_folder = os.path.join(BASE_DIR, lipase_name)
		zip_folder = os.path.join(CIF_DIR, lipase_name)

		result = process_lipase_data(pdb_folder, zip_folder, lipase_name)
		#lipase_name, plot_df, total_models = process_lipase_data(folder, lipase_name)
		if result is not None:
			all_data.append(result)
	
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
		
		model_threshold = int(total_models * FREQUENCY_THRESHOLD) # HAVE TO APPEAR IN > 0.33 OF THE TOTAL MODELS
		valid_keys = plot_df[plot_df["ResKey"].isin(shared_keys)]
		valid_keys = valid_keys[valid_keys["Count"] >= model_threshold]
		final_highlight_keys = set(valid_keys["ResKey"])
		
		# HIGHLIGHT POCKET RESIDUES
		residue_numbers = plot_df["Residue"].apply(lambda r: int(r.split("_")[1]))
		pocket_ranges = POCKET_RANGES.get(lipase_name, set())

		def is_pocket_residue(label, pocket_set):
			"""
			Determine if a given residue label (e.g., 'D:GLU_242') is in the provided pocket set.

			Args:
				label (str): Residue label in format 'Chain: RESNAME_RESNUM'
				pocket_set (set): Set of integer residue numbers known to be pocket residues

			Returns:
				bool: True if residuee number is in pocket_set
			"""
			try:
				resnum = int(label.split(":")[1].split("_")[1])
				return resnum in pocket_set
			except Exception as e:
				print(f"[ERROR] Failed parsing residue label {label}: {e}")
				return False

		plot_df["IsPocket"] = plot_df["Residue"].apply(lambda label: is_pocket_residue(label, pocket_ranges))
			
		#highlight_flags = plot_df["ResKey"].isin(final_highlight_keys)
		highlight_flags = residue_numbers.isin(pocket_ranges)
	
		print(f"[DEBUG] {lipase_name} \n{plot_df[highlight_flags]}\n*************************")
		
		bar_colors = [
			SHARED_COLOR if is_pocket else chain_color_map[row["Chain"]]
			for is_pocket, (_, row) in zip(plot_df["IsPocket"], plot_df.iterrows())
		]

		# highlight_edge_colors = ["black" if row["IsPocket"]  else "none" for _ in plot_df.iterrows()]
	
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
		for tick, is_pocket in zip(ax.get_xticklabels(), plot_df["IsPocket"]):
			tick.set_color(SHARED_COLOR if is_pocket else "black")
		ax.set_title(f"{lipase_name} Interface Residues", fontsize=20)	
		ax.set_ylabel("Percent of Models Interacting", fontsize=14)
		ax.set_ylim(0,100)
		ax.tick_params(axis='y', labelsize=15)
		ax.set_xlabel("Residue")
	
	handles = [plt.Rectangle((0,0),1,1, color=chain_color_map[ch]) for ch in EXPECTED_CHAINS]
	fig.legend(handles, [f"Chain {ch}" for ch in EXPECTED_CHAINS], title="Chain", loc='upper right')
	plt.tight_layout(rect=[0, 0, 0.95, 0.98])
	highlight_patch = plt.Rectangle((0,0),1,1, facecolor='white', edgecolor='black', linewidth=1.2, label='Shared across chains')
	handles.append(highlight_patch)
	
	filename = f"lipase_interaction_occurrence.png"
	out_path = os.path.join(PLOT_DIR, filename)
	plt.savefig(out_path, dpi=300)
	plt.close()
	print(f"Saved combined plot to: {out_path}")
	print("Batch processing complete")
	
if  __name__ == "__main__":
	main()
	
