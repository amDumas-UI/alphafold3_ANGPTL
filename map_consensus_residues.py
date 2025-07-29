"""
map_consensus_residues.py

Author: Alexander Dumas  
Last Updated: July 2025

Description:
	This script identifies and maps consensus interface residues from predicted lipase–heterotrimer interactions
	onto representative AlphaFold structures. It highlights residues that appear in more than one-third of 
	interface models and visualizes them using PyMOL.

	Two PyMOL sessions are generated:
		1. A session showing only the lipase (chain D) with highlighted interface residues.
		2. A session showing only the heterotrimer (chains A, B, C) with highlighted residues from chains A–C 
         	that interact with the lipase.

Inputs:
	- Interface residue .pdb files (from individual predictions):
		interfaces/heterotrimer/[lipase_folder]/*.pdb

	- Reference AlphaFold .cif files:
		Lipase structure:        alphafold_results/lipase/[lipase_folder]/ranked_0.cif  
		Heterotrimer structure:  alphafold_results/heterotrimer/[lipase_folder]/ranked_0.cif

Outputs:
	- Colored .pdb file with highlighted residues:
        	interfaces/heterotrimer/[lipase_folder]/highlighted_interface.pdb

	- PyMOL session files:
        	plots/heterotrimer/[lipase_folder]/LPL_lipase_only.pse  
        	plots/heterotrimer/[lipase_folder]/LPL_heterotrimer_only.pse

Dependencies:
	- Python 3.8+
	- PyMOL (with Python API access)
	- Standard Python libraries: os, sys, zipfile, collections

Usage:
	$ python map_consensus_residues.py

    (Ensure working directory contains the correct folder structure for input files.)
"""

#---- CONFIGURATION ----

import pymol
import os
import zipfile
import tempfile
from collections import Counter, defaultdict
from pymol import cmd

PDB_DIR = "./interfaces/homotrimer/A3/"
TRIMER_ZIP_DIR = "./alphafold_results/Trimer_alone/A3/"
LIPASE_ZIP_DIR = "./lipase_predictions/"
OUTPUT_DIR = "./plots/homotrimer/A3/"
THRESHOLDS = [
	(0.9, "red"),	# Highest confidence
	(0.6, "orange"),
	(0.3, "yellow"),	# 1/3 of total models. Similar to threshold I have used in the other scripts
	# < 0.3 = white
]

os.makedirs(OUTPUT_DIR, exist_ok=True)


def extract_interface_residues(pdb_dir):
	"""
	Extracts interface residues that appear in more than 'threshold' fraction of models.
	
	Args:
		pdb_dir (str): Path to directory containing interface .pdb files.
	Returns:
		set[str, float]: Residue identifiers and their frequency (0-1) (e.g., D/123) that meet the thresholds.
	"""
	model_residue_map = defaultdict(set)
	for pdb_file in os.listdir(pdb_dir):
		if not pdb_file.endswith(".pdb"):
			continue
		with open(os.path.join(pdb_dir, pdb_file)) as f:
			for line in f:
				if line.startswith("ATOM"):
					resn = line[17:20].strip()
					chain = line[21].strip()
					resi = line[22:26].strip()
					res_id = f"{chain}:{resn}_{resi}"
					model_residue_map[pdb_file].add(res_id)

	all_residues = Counter()
	for residues in model_residue_map.values():
		all_residues.update(residues)

	num_models = len(model_residue_map)
	residue_frequencies = {
		res: count / num_models for res, count in all_residues.items()}
	
	return residue_frequencies


def extract_first_cif(zip_folder):
	"""
	Unzips and extracts the first .cif file found inside a .zip archive.

	Args:
		zip_folder (str): Path to parent directory containing lipase folders with interface .pdb files.

	Returns:
		
	"""
	zip_files = [f for f in os.listdir(zip_folder) if f.endswith(".zip")]
	if not zip_files:
		print(f"[WARNING] No zip files in {zip_folder}")
		return None

	zip_path = os.path.join(zip_folder, zip_files[0])
	with zipfile.ZipFile(zip_path, 'r') as zip_ref:
		cif_files = [f for f in zip_ref.namelist() if f.endswith(".cif")]
		if not cif_files:
			print(f"[WARNING] No .cif found in {zip_folder}")
			return None
		
		temp_dir = os.path.join("temp_extracted")
		os.makedirs(temp_dir, exist_ok=True)
		extracted_file = zip_ref.extract(cif_files[0], path=temp_dir)
		
		return extracted_file	

def color_and_export(cif_path, highlighted_residues, output_path, chain_filter):
	"""
	Loads structures in PYMOL and color residues matching the highlighted list.
	Saves two .pse files: one full trimer, one lipase only.

	Args:
		trimer_cif (str): Path to trimer .cif reference structure
		lipase_cif (str): Path to lipase .cif reference structure
		highlighted_residues (set[str]): Residues to highlight (e.g., D/123).
		lipase_name (str): Name of the lipase to be plotted
	"""
	cmd.reinitialize()
	cmd.load(cif_path, "af_model")
	print(f"[DEBUG] Loaded  model from: {cif_path}")
	
	chains = cmd.get_chains("af_model")
	print(f"[DEBUG] Chains in model: {chains}")

	# Select chains of interest
	chain_selection = " + ".join([f"chain {ch}" for ch in chain_filter])
	cmd.select("visible_chains", f"(af_model and ({chain_selection}))")
	cmd.create("my_model", "visible_chains")
	if cmd.count_atoms("my_model") == 0:
		print(f"[ERROR] no atoms selected for chains {chain_filter} in {cif_path}")
		return
	cmd.delete("af_model")

	# Highlight residues
	def get_color_for_frequency(freq):
		for threshold, color in THRESHOLDS:
			if freq >= threshold:
				return color
		return "white" # default fallback

	print(f"[DEBUG] Total residues to highlight on lipase: {len(highlighted_residues)}")
	
	print(sorted(list(highlighted_residues.items())[:5]))
	for res_id, freq in highlighted_residues.items():
		if res_id.startswith("D:"):
			color = get_color_for_frequency(freq)
			chain, rest = res_id.split(":")
			resn, resi = rest.split("_")
			translated_chain = "A" # Since using a lipase monomer the chain is now chain A in the lipase only pse
			selection = f"my_model and chain {translated_chain} and resi {resi} and resn {resn}"
			cmd.select("residue_sel", selection)
			cmd.show("cartoon", "residue_sel")
			cmd.color(color, "residue_sel")
			cmd.delete("residue_sel")
	
	for res_id, freq in highlighted_residues.items():
		if res_id.startswith(("A:","B:","C:")):
			color = get_color_for_frequency(freq)
			chain, rest = res_id.split(":")
			resn, resi = rest.split("_")
			selection = f"my_model and chain {chain} and resi {resi} and resn {resn}"
			cmd.select("residue_sel", selection)
			cmd.show("cartoon", "residue_sel")
			cmd.color(color, "residue_sel")
			cmd.delete("residue_sel")

	cmd.save(output_path)
	print(f"[INFO] Saved PyMOL session: {output_path}")


def main():
	"""
	Batch-process all lipase folders: map shared interface residues onto .cif model and export colored .pse.
	"""
	os.makedirs(OUTPUT_DIR, exist_ok=True)
	
	# Load trimer .cif once
	trimer_cif = extract_first_cif(TRIMER_ZIP_DIR)
	if not trimer_cif:
		print(f"[ERROR] Could not extract trimer .cif from {TRIMER_ZIP_DIR}")
		return
	
	# Loop over each lipase folder
	for lipase_folder in os.listdir(PDB_DIR):
		pdb_dir = os.path.join(PDB_DIR, lipase_folder)
		if not os.path.isdir(pdb_dir):
			print(f"[WARNING] No pdb files found for {lipase_folder}. Skipping.")
			continue
		
		print(f"\n[PROCESSING] {lipase_folder}")

		consensus_residues = extract_interface_residues(pdb_dir)
		print(f" -> {len(consensus_residues)} residues found within the thresholds")

		lipase_zip_path = os.path.join(LIPASE_ZIP_DIR, lipase_folder)
		lipase_cif = extract_first_cif(lipase_zip_path)

		if not lipase_cif or not trimer_cif:
			print(f"[SKIP] Missing .cif structure for {lipase_folder}")
			continue

		# Output paths
		lipase_out = os.path.join(OUTPUT_DIR, f"{lipase_folder}_lipase_only.pse")
		trimer_out = os.path.join(OUTPUT_DIR, f"{lipase_folder}_trimer_only.pse")

		# Make lipase-only visualization (Chain D)
		lipase_residues = {res: freq for res, freq in consensus_residues.items() if res.startswith("D:")}
		print(f"[DEBUG] Residues on chain D (lipase): {len(lipase_residues)}")
		color_and_export(lipase_cif, lipase_residues, lipase_out, chain_filter=["A"])

		# Make trimer-only visualization (chains A/B/C)
		trimer_residues = {res: freq for res, freq in consensus_residues.items() if res.startswith(("A:","B:","C:"))}
		print(f"[DEBUG] Residues on chains A, B, C (trimer chains): {len(trimer_residues)}")
		color_and_export(trimer_cif, trimer_residues, trimer_out, chain_filter=["A", "B", "C"])
	
	print("\n[FINISH] Batch processesing complete.")

if __name__ == "__main__":
	main()

