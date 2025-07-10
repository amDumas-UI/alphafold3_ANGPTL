"""
batch_interface.py

Author: Alexander Dumas
Date: 2025-06-13
Purpose: This script processes AlphaFold-predicted structures of homotrimers and lipases.
It extracts .cif files from zip archives, detects interacting residues between homotrimer chains and a lipase,
and saves interface residues as new PDB files using PyMOL and a custom interfaceResidues function.

Dependencies:
- PyMOL with Python API
- InterfaceResidues.py in Python path

Inputs:
- Directory containing AlphaFold output zip files
- Sequences assumed to have consistent chain IDs per model

Outputs:
- Interface PDB files saved in `output_dir`

Usage:
Run inside a PyMOL-enabled Python environment.

Note:
This script is designed for homotrimer-lipase systems. For ANGPTL heterotrimers, see the alternative script (ht_batch_interface.py).
"""

import os 
import shutil
import zipfile
from pymol import cmd
from InterfaceResidues import interfaceResidues
from collections import defaultdict


# ---- CONFIGURATION ----
input_dir = "/Volumes/rdss_bsdavies/Lab_Notes_and_Data/A_Dumas_Notes_and_Data/alphafold3_automation/alphafold_results/Altered_chains/A4_LPL/"
output_dir = "interfaces/homotrimer/A4_LPL_Altered/"
temp_dir = "temp_cif"
cutoff = 6 # CUTOFF FOR HOW MUCH SURFACE AREA IS LOST TO INTERACTION. SET IN SQUARE ANGSTROMS
os.makedirs(output_dir, exist_ok = True)

def get_chain_sequence(model, chain):
	"""
	Extract 1-letter amino acid sequence of a chain from a PyMOL object.

	Args:
		model (str): Name of the loaded PyMOL object.
		chain (str): Chain identifier.

	Returns:
		str: 1-letter amino acid sequence
	"""
	seq = []
	cmd.iterate(f"{model} and chain {chain} and name CA", "seq.append(resn)", space={'seq': seq})

	# Convert 3-letter to 1-letter codes
	three_to_one = {
		'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C',
		'GLU':'E','GLN':'Q','GLY':'G','HIS':'H','ILE':'I',
		'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P',
		'SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V',
	}
	return ''.join([three_to_one.get(res.upper(), 'X') for res in seq])


def group_chains_by_sequence(model_name, threshold=0.99):
	"""
	Group chains in a PyMOL object by sequence identity.
	
	Args:
		model_name (str): Name of the loaded PyMOL object.
		threshold: Minimum fraction of identity for chains to be grouped.
	
	Returns: 
		A list of chain groups, where each group is a list of chain IDs.
	"""

	chains = cmd.get_chains(model_name)
	seqs = {chain: get_chain_sequence(model_name, chain) for chain in chains}

	def seq_identity(seq1, seq2):
		if len(seq1) != len(seq2) or len(seq1) == 0:
			return 0
		matches = sum(a == b for a, b in zip(seq1, seq2))
		return matches / len(seq1)

	groups = []
	for chain, seq in seqs.items():
		placed = False
		for group in groups:
			# COMPARE WITH FIRST CHAIN SEQUENCE IN THE GROUP
			group_seq = seqs[group[0]]
			if seq_identity(seq,group_seq) >= threshold:
				group.append(chain)
				placed = True
				break
		if not placed:
			groups.append([chain])
	return groups



# LOOP OVER ALL ZIP FILES
def main():
	"""
	Main function to loop through lipase directories containing .zips and extract .cif models for structure analysis. Call InterfaceResidues() to calculate the residues within a specific interaction distance based on ASA calculations
	"""

	print(f"Starting PyMol batch interface script with a cutoff of {cutoff}")
	
	# LOOP OVER ALL ALPHAFOLD RESULT ZIP FILES
	for zipname in os.listdir(input_dir): 
		if not zipname.endswith(".zip"):
			continue
	
		print(f"Processing archive: {zipname}")
	
		# CLEAN TEMP FOLDER BEFORE EXTRACTION
		if os.path.exists(temp_dir):
			shutil.rmtree(temp_dir)
		os.makedirs(temp_dir)

		zip_path = os.path.join(input_dir, zipname)
		with zipfile.ZipFile(zip_path, 'r') as zip_ref:
			zip_ref.extractall(temp_dir)
		print("[DEBUG]: Extracted files: ", os.listdir(temp_dir))
	
		# PROCESS ALL .CIF FILES IN TEMP FOLDER
		for fname in os.listdir(temp_dir):
			if not fname.endswith(".cif"):
				continue
		
			cif_path = os.path.join(temp_dir, fname)
			model_name = os.path.splitext(fname)[0]
			print(f"[DEBUG]:  Loading CIF: {fname}")
			cmd.load(cif_path, model_name)
		
			try: 
				groups = group_chains_by_sequence(model_name)
				groups_sorted = sorted(groups, key=len, reverse=True)
				
				if len(groups_sorted) < 2:
					print(f"[WARNING]: Only one chain group found. Skipping {model_name}.")
					cmd.delete("all")
					continue

				homotrimer_chains = groups_sorted[0]
				lipase_chains = [ch for group in groups_sorted[1:] for ch in group]

				print(f"[DEBUG]: Detected homotrimer chains: {homotrimer_chains}")
				print(f"[DEBUG]: Detected lipase chains: {lipase_chains}")
			
				sel1 = f"chain {'+'.join(homotrimer_chains)}"
				sel2 = f"chain {'+'.join(lipase_chains)}"
				combined_int_name = f"{model_name}_int_{''.join(homotrimer_chains)}_{''.join(lipase_chains)}"
		
				cmd.create("MyComplex", f"{model_name} and chain {'+'.join(homotrimer_chains + lipase_chains)}")
				
				print(f"[DEBUG] Selection 1 atom count: {cmd.count_atoms(sel1)}")
				print(f"[DEBUG] Selection 2 atom count: {cmd.count_atoms(sel2)}")
		
				# CALL HELPER FUNCTION THAT CALCULATES THE RESIDUES INTERACTING BASED ON ASA CALCULATIONS	
				interfaceResidues("MyComplex", sel1, sel2, cutoff, combined_int_name)
	
				count = cmd.count_atoms(combined_int_name)
				print(f" -> Interface atom count: {count}")
				if count == 0:
					print("[WARNING] No atoms found in interface selection! Saving entire chain groups instead.")
					fallback_sel = f"{model_name}_fallback_{''.join(homotrimer_chains)}_{''.join(lipase_chains)}"
					cmd.select(fallback_sel, f"MyComplex and ({sel1} or {sel2})")
					out_path = os.path.join(output_dir, f"{fallback_sel}.pdb")
					cmd.save(out_path, fallback_sel)
					print(f"[FALLBACK] Saved fallback chain group  to {out_path}")
				else:
					out_path = os.path.join(output_dir, f"{combined_int_name}.pdb")
					cmd.save(out_path, combined_int_name)
					print(f"[SAVE]: Saved interface selection to {out_path}")

			except Exception as e:
				print(f"[ERROR]: Failed combined interfaceResidues (A+B+C vs D): {e}")
			
			cmd.delete("all") # RESET FOR NEXT STRUCTURE	

if __name__ =="__main__":
	main()
