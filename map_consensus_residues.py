"""
map_consensus_residues.py

Author: Alexander Dumas  
Last Updated: August 2025

Description:
	This script identifies and maps consensus interface residues from predicted ANGPTL-target interactions
	onto representative AlphaFold structures within PyMOL. It highlights residues based on a threshold of 
	the percent of interaction models that contain that interaction. 

	Two PyMOL sessions are generated:
		1. A session showing only the target with highlighted interface residues.
		2. A session showing only the ANGPTL complex with highlighted residues that interact with the target. 

Inputs:
	- Interface residue .pdb files (from individual predictions):
		interfaces/path/to/*.pdb

	- Reference AlphaFold .cif files:
		Target structure:  alphafold_results/path/to/.zip  
		ANGPTL structure:  alphafold_results/path/to/.zip

Outputs:
	- PyMOL session files:

Dependencies:
	- Python 3.8+
	- PyMOL (with Python API access)
	- Standard Python libraries: os, sys, zipfile, collections

Usage:
	$ python map_consensus_residues.py
	
Notes:
	Make sure to alter the thresholds values if needed.
	Make sure to change which chains are present within the ANGPTL and target structures based on the complex organization.

    (Ensure working directory contains the correct folder structure for input files.)
"""

#---- CONFIGURATION ----

import pymol
import os
import zipfile
import tempfile
import shutil
from collections import Counter, defaultdict
from pymol import cmd

PDB_DIR = "./interfaces/integrin/A3/cleaved_monomer/" # Folder containing .pdb interface models (not zipped)
ANGPTL_ZIP_DIR = "./alphafold_results/monomer/A3/cleaved_monomer/" # Folder containing a single a .zip with a .cif
TARGET_ZIP_DIR = "./alphafold_results/integrin_predictions/integrin_a5B3/"
OUTPUT_DIR = "./plots/integrin/A3/cleaved_monomer/"
THRESHOLDS = [
	(0.9, "red"),	# Highest confidence
	(0.6, "orange"),
	(0.3, "yellow"),	# 1/3 of total models. Similar to threshold I have used in the other scripts
	# < 0.3 = white
]
ANGPTL_CHAINS = ("A",)		# Alter these based on complex organization
TARGET_CHAINS = ("B","C")

os.makedirs(OUTPUT_DIR, exist_ok=True)


def extract_interface_residues(pdb_dir):
	"""
	Extracts interface residues that appear in more than 'threshold' fraction of models.
	
	Args:
		pdb_dir (str): Path to directory containing interface .pdb files.
	Returns:
		dict[str, float]: Residue identifiers and their frequency (0-1) (e.g., D/123) that meet the thresholds.
	"""
	model_residue_map = defaultdict(set)
	
	pdb_files = [f for f in os.listdir(pdb_dir) if f.endswith(".pdb")]
	print(f"[DEBUG] Found {len(pdb_files)} PDB files in {pdb_dir}")
	
	for pdb_file in pdb_files:
		file_path = os.path.join(pdb_dir, pdb_file)
		print(f"[DEBUG] Processing {pdb_file}")
	
		with open(file_path) as f:
			count_residues = 0
			for line in f:
				if line.startswith("ATOM"):
					resn = line[17:20].strip()
					chain = line[21].strip()
					resi = line[22:26].strip()
					res_id = f"{chain}:{resn}_{resi}"
					model_residue_map[pdb_file].add(res_id)
					count_residues += 1
			print(f"[DEBUG] Extracted {count_residues} residues from {pdb_file}")

	all_residues = Counter()
	for residues in model_residue_map.values():
		all_residues.update(residues)

	num_models = len(model_residue_map)
	print(f"[DEBUG] Total models processed: {num_models}")
	print(f"[DEBUG] Total unique residues across all models: {len(all_residues)}")

	residue_frequencies = {res: count / num_models for res, count in all_residues.items()}
	
	return residue_frequencies


def extract_first_cif(zip_folder):
	"""
	Unzips and extracts the first .cif file found inside a .zip archive.

	Args:
		zip_folder (str): Path to .zip archive containing .cif models.

	Returns:
		extracted_file (str): Path to the extracted .cif file in the temporary directory.
	"""
	print(f"[DEBUG] Starting extraction from folder: {zip_folder}")

	zip_files = [f for f in os.listdir(zip_folder) if f.endswith(".zip")]
	if not zip_files:
		raise FileNotFoundError(f"No zip files in {zip_folder}")
	
	print(f"[DEBUG] ZIP files found in {zip_folder}: {zip_files}")
	zip_path = os.path.join(zip_folder, zip_files[0])
	print(f"[DEBUG] using ZIP files: {zip_path}")
	
	# Create a temporary directory that cleans up automatically
	temp_dir_obj = tempfile.TemporaryDirectory()
	temp_dir = temp_dir_obj.name
	print(f"[DEBUG] Created temp directory: {temp_dir}")
	
	with zipfile.ZipFile(zip_path, 'r') as zip_ref:
		cif_files = [f for f in zip_ref.namelist() if f.endswith(".cif") and "_model_" in f]
		print(f"[DEBUG] CIF files in zip: {cif_files}")
		if not cif_files:
			raise FileNotFoundError(f"No .cif files found in {zip_folder}")
		
		# Extract the first model .cif file
		extracted_file = zip_ref.extract(cif_files[0], path=temp_dir)
		print(f"[DEBUG] Extracted CIF file: {extracted_file}")

	# Check that extracted_file actually exists on disk
	if not os.path.isfile(extracted_file):
		print(f"[ERROR] Extracted file does NOT exist: {extracted_file}")
	else:
		print(f"[DEBUG] Verified extracted file exists.")

	return extracted_file, temp_dir_obj 	

def color_and_export(cif_path, highlighted_residues, output_path, chain_filter):
	"""
	Loads structures in PYMOL and color residues matching the highlighted list.
	Saves two .pse files: one ANGPTL only, one target only.

	Args:
		cif_path (str): Path to .cif reference structure.
		highlighted_residues (dict[str, int]): Residues to highlight (e.g., D/123).
		output_path (str): Path to desired output for .pse models
		chain_filter (str): Which chains to collect interaction residues from.
	"""
	
	def build_chain_map_dynamic(chain_filter, cif_chains):
		"""
		Map original chain IDs to CIF chain IDs dynamically based on order.

		Args:
			chain_filter (tuple[str]): Chain IDs from the original residue data.
			cif_chains (list[str]): Chain IDs present in the loaded CIF reference model.
	
		Returns:
			dict[str, int]: Mapping from original chain IDs to CIF reference chain IDs.
		"""
		if len(chain_filter) != len(cif_chains):
			print(f"[WARNING] chain_filter ({len(chain_filter)}) and cif_chains ({len(cif_chains)}) differ in length.")
		return {orig: new for orig, new in zip(chain_filter, cif_chains)}

	def remap_chain_ids_multi(residues, chain_map):
		"""
		Replace chain IDs in residues according to chain_map.

		Args:
			residues (dict[str, int]): Mapping of residue IDs (e.g., "D:123") to frequencies.
			chain_map (dict[str, int]): Mapping from original chain IDs to CIF chain IDs.

		Returns:
			dict[str, int]: Residue mapping with chain IDs replaced by those from the CIF.
		"""
		remapped = {}
		for res, freq in residues.items():
			try:
				chain_id, rest = res.split(":", 1) # split at firt colon
			except ValueError:
				print(f"[WARNING] Unexpected residue format: {res}")
				continue

			new_chain = chain_map.get(chain_id, chain_id) # default to original if not found
			new_res = f"{new_chain}:{rest}"
			remapped[new_res] = freq

		return remapped
	
	print(f"[DEBUG] Starting color_and_export with cif_path: {cif_path}")
	print(f"[DEBUG] Number of residues to highlight: {len(highlighted_residues)}")
	
	cmd.reinitialize()
	cmd.delete("af_model")
	print(f"[DEBUG] Objects in PyMOL before load (should be none): {cmd.get_object_list()}")
	
	print(f"[DEBUG] Absolute cif_path: {os.path.abspath(cif_path)}")
	print(f"[DEBUG] File exists: {os.path.exists(cif_path)}, Size: {os.path.getsize(cif_path) if os.path.exists(cif_path) else 'N/A'} bytes")
	cmd.load(cif_path, "af_model")

	print(f"[DEBUG] Objects in PyMOL after load: {cmd.get_object_list()}")
	chains = cmd.get_chains("af_model")
	print(f"[DEBUG] Chains in loaded model: {chains}")
	print(f"[DEBUG] Chains specified in chain_filter: {chain_filter}")

	# Translate chain ids from .pdb to match .cif chain ids
	if list(chains) != list(chain_filter):
		print(f"[DEBUG] Chain filter mismatch reason: {list(chains)} != {list(chain_filter)}")
		chain_map = build_chain_map_dynamic(chain_filter, chains)
		print(f"[DEBUG] Built chain map for mismatch case: {chain_map}")
	else:
		print("Already matching... no change to chain filter")
		chain_map = {ch: ch for ch in chain_filter}
	
	# Remap residues
	highlighted_residues = remap_chain_ids_multi(highlighted_residues, chain_map)
	print(f"[DEBUG] Remapped residues (first 10):")
	for i, (res, freq) in enumerate(highlighted_residues.items()):
		if i >= 10:
			break
		print(f"{res}: {freq}")	

	# Select chains of interest
	mapped_chains = list(chain_map.values())
	chain_selection = " + ".join([f"chain {ch}" for ch in mapped_chains])
	
	cmd.select("visible_chains", f"(af_model and ({chain_selection}))")
	print(f"[DEBUG] PyMOL selection string: (af_model and ({chain_selection}))")
	print(f"[DEBUG] Atom count in visible_chains selection: {cmd.count_atoms('visible_chains')}")

	cmd.create("my_model", "visible_chains")

	print(f"[DEBUG] Objects after creating my_model: {cmd.get_object_list()}")
	print(f"[DEBUG] Chains in my_model; {cmd.get_chains('my_model')}")

	if cmd.count_atoms("my_model") == 0:
		print(f"[ERROR] no atoms selected for chains {chain_filter} in {cif_path}")
		return
	cmd.delete("af_model")

	# Highlight residues
	def get_color_for_frequency(freq):
		"""
		Determine PyMOL color name for a residue based on interaction frequency.
		
		Args:
			freq (int): Number of times the residue appears across models.

		Returns:
			str: PyMOL color name.
		"""

		for threshold, color in THRESHOLDS:
			if freq >= threshold:
				return color
		return "white" # default fallback

	print(f"[DEBUG] Total residues to highlight: {len(highlighted_residues)}")
	
	print(sorted(list(highlighted_residues.items())[:5]))

	for res_id, freq in highlighted_residues.items():
		color = get_color_for_frequency(freq)
		chain, rest = res_id.split(":")
		resn, resi = rest.split("_")
		chain = chain
		selection = f"my_model and chain {chain} and resi {resi} and resn {resn}"
		cmd.select("residue_sel", selection)
		cmd.show("cartoon", "residue_sel")
		cmd.color(color, "residue_sel")
		cmd.delete("residue_sel")

	cmd.save(output_path)
	print(f"[INFO] Saved PyMOL session: {output_path}")


def main():
	"""
	Single directory processing of ANGPTL interaction with target
	"""
	os.makedirs(OUTPUT_DIR, exist_ok=True)
	
	# Load ANGPTL .cif once
	angptl_cif, angptl_temp_dir = extract_first_cif(ANGPTL_ZIP_DIR)
	target_cif, target_temp_dir = extract_first_cif(TARGET_ZIP_DIR)
	print("Extracted TARGET cif:", target_cif)
	with open(target_cif) as f:
		for i, line in enumerate(f):
			if i < 20:
				print(line.strip())
	
	if not angptl_cif:
		print(f"[ERROR] Could not extract ANGPTL .cif from {ANGPTL_ZIP_DIR}")
		return
	if not target_cif:
		print(f"[ERROR] Could not extract target .cif from {TARGET_ZIP_DIR}")
		
	print(f"\n[PROCESSING] ANGPTL CIF path: {angptl_cif}")
	print(f"\n[PROCESSING] TARGET CIF path: {target_cif}")

	# Pring file sizes
	print(f"[DEBUG] ANGPTL CIF size: {os.path.getsize(angptl_cif)} bytes")
	print(f"[DEBUG] TARGET CIF size: {os.path.getsize(target_cif)} bytes")

	# Printing first 20 lines of each cif file for more debugging
	with open(angptl_cif) as f:
		print(f"[DEBUG] First 20 lines of ANGPTL CIf:")
		for i, line in enumerate(f):
			print(line.strip())
			if i >= 19:
				break
	
	with open(target_cif) as f:
		print(f"[DEBUG] First 20 lines of TARGET CIF:")
		for i, line in enumerate(f):
			print(line.strip())
			if i >= 19:
				break

	consensus_residues = extract_interface_residues(PDB_DIR)
	print(f"[INFO] -> {len(consensus_residues)} residues found within the thresholds")
	
	# Split residues by complex
	angptl_residues = {res: freq for res, freq in consensus_residues.items()
		if res.startswith((ANGPTL_CHAINS))}
	target_residues = {res: freq for res, freq in consensus_residues.items()
		if res.startswith((TARGET_CHAINS))}

	# Output paths
	target_out = os.path.join(OUTPUT_DIR, "target_only.pse")
	angptl_out = os.path.join(OUTPUT_DIR, "angptl_only.pse")

	# Make target-only visualization 
	print(f"[DEBUG] Residues on target chains: {len(target_residues)}")
	color_and_export(target_cif, target_residues, target_out, chain_filter=TARGET_CHAINS)

	# Make ANGPTL-only visualization
	print(f"[DEBUG] Residues on ANGPTL chains: {len(angptl_residues)}")
	color_and_export(angptl_cif, angptl_residues, angptl_out, chain_filter=ANGPTL_CHAINS)
	
	angptl_temp_dir.cleanup()
	target_temp_dir.cleanup()
	
	print("\n[FINISH] Processesing complete.")

if __name__ == "__main__":
	main()

