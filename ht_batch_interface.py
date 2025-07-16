"""
ht_batch_interface.py

Author: Alexander Dumas
Date: 2025-06-30

Purpose:
    This script processes AlphaFold structure predictions for heterotrimeric ANGPTL complexes interacting with a lipase.
    It extracts CIF files from AlphaFold output zip archives, loads them into PyMOL, and identifies interacting residues
    between the ANGPTL heterotrimer (chains A, B, C) and the lipase (chain D). Output includes interface PDBs and plots.

Inputs:
    - Directory of folders, each containing one or more AlphaFold result zip files.
    - CIF files must follow the naming and chain-labeling conventions described in the methods.

Outputs:
    - CSV and bar plots of lipase interface residues interacting with chains A, B, or C.
    - Results organized per lipase folder in the specified output directory.

Dependencies:
    - PyMOL with Python API (cmd)
    - InterfaceResidues.py
    - matplotlib, pandas

Run:
    $ python ht_batch_interface.py 

Note:
This script is designed for heterotrimer-lipase systems due to hardcoded chain identifiers within the script. 
This script expects chains A-C to be the ANGPTL heterotrimer and chain D to be the lipase. Alteration of this is accepted with the script.
For ANGPTL homotrimers, see the alternative script (batch_interface.py) which identifies trimer chains based on sequence similarity.
"""

import os 
import shutil
import zipfile
from pymol import cmd
from InterfaceResidues import interfaceResidues
from collections import defaultdict

# ---- CONFIGURATION ----
input_dir = "/Volumes/rdss_bsdavies/Lab_Notes_and_Data/A_Dumas_Notes_and_Data/alphafold3_automation/alphafold_results/heterotrimer/A8_A5/A8x1_A5x2/"  # DIRECTORY CONTAINING LIPASE FOLDERS WITH .ZIP FILES OF STRUCTURE PREDICTIONS
output_root = "interfaces/heterotrimer/A8_A5/A8x1_A5x2/"
temp_dir = "temp_cif"
cutoff = 6 # ANGSTROM SQUARED DIFFERENCE IN ASA. CUTOFF FOR HOW MUCH SURFACE AREA IS LOST TO INTERACTION. HIGHER THE CUTOFF THE MORE STRICT 
os.makedirs(output_root, exist_ok = True)
heterotrimer_chains = ['A','B','C']
lipase_chain = ['D']

def main():
	"""
	Main function that will loop through each .zip archive for each lipase prediction and extract the .cif models. This function will call the helper function InterfaceResidues() to use ASA calculations to collect interacting residues for .pdb model creation
	"""
	print(f"Starting PyMol batch interface script with a cutoff of {cutoff}")

	# LOOP THROUGH EACH LIPASE FOLDER
	for subfolder in sorted(os.listdir(input_dir)):
		subfolder_path = os.path.join(input_dir, subfolder)
		if not os.path.isdir(subfolder_path):
			continue
	
		print(f"\n--- Processing lipase folder: {subfolder} ---")
		output_dir = os.path.join(output_root, subfolder)
		os.makedirs(output_dir, exist_ok=True)

		zip_files = [f for f in os.listdir(subfolder_path) if f.endswith("zip")]
		if not zip_files:
			print(f"WARNING: No zip files found in: {subfolder_path}")
			continue
		
		# LOOP THROUGH EACH ZIP ARCHIVE INSIDE THE LIPASE FOLDERS
		for zipname in zip_files: 
			print(f"Processing archive: {zipname}")
			zip_path = os.path.join(subfolder_path, zipname)

			# CLEAN TEMP FOLDER BEFORE EXTRACTION
			if os.path.exists(temp_dir):
				shutil.rmtree(temp_dir)
			os.makedirs(temp_dir)

			with zipfile.ZipFile(zip_path, 'r') as zip_ref:
				zip_ref.extractall(temp_dir)
			print("Extracted files: ", os.listdir(temp_dir))

			# PROCESS ALL .CIF FILES IN TEMP FOLDER
			for fname in os.listdir(temp_dir):
				if not fname.endswith(".cif"):
					continue
		
				cif_path = os.path.join(temp_dir, fname)
				model_name = os.path.splitext(fname)[0]
		
				print(f"[DEBUG]: Loading CIF: {fname}")
				cmd.load(cif_path, model_name)
		
				try: 
					print(f"[DEBUG]: Detected heterotrimer chains: {heterotrimer_chains}")
					print(f"[DEBUG]: Detected lipase chains: {lipase_chain}")
			
					sel1 = f"chain {'+'.join(heterotrimer_chains)}"
					sel2 = f"chain {'+'.join(lipase_chain)}"
					combined_int_name = f"{model_name}_int_{''.join(heterotrimer_chains)}_{''.join(lipase_chain)}"
		
					cmd.create("MyComplex", f"{model_name} and chain {'+'.join(heterotrimer_chains + lipase_chain)}")
					
					print(f"[DEBUG] Selection 1 atom count: {cmd.count_atoms(sel1)}")
					print(f"[DEBUG] Selection 2 atom count: {cmd.count_atoms(sel2)}")
			
					# CALL INTERFACE RESIDUE HELPER FUNCTION
					interfaceResidues("MyComplex", sel1, sel2, cutoff, combined_int_name)
	
					count = cmd.count_atoms(combined_int_name)
					print(f" -> Interface atom count: {count}")
			
					if count == 0:
						print("[WARNING] No atoms found in interface selection! Saving entire chain groups instead.")
						fallback_sel = f"{model_name}_fallback_{''.join(heterotrimer_chains)}_{''.join(lipase_chain)}"
						cmd.select(fallback_sel, f"MyComplex and ({sel1} or {sel2})")
						out_path = os.path.join(output_dir, f"{fallback_se1}.pdb")
						cmd.save(out_path, fallback_sel)
						print(f"[FALLBACK] Saved fallback chain group  to {out_path}")
					else:
						out_path = os.path.join(output_dir, f"{combined_int_name}.pdb")
						cmd.save(out_path, combined_int_name)
						print(f"[SAVE]: Saved interface selection to {out_path}")

				except Exception as e:
					print(f"[ERROR]: Failed combined interfaceResidues for {model_name}: {e}")
		
				cmd.delete("all") # RESET FOR NEXT STRUCTURE	

	print("\nAll lipase folders processed.")

if __name__ == "__main__":
	main()
