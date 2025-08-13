"""
ANGPTL_alphafold_submitter_playwright.py

Author: Alexander Dumas  
Last updated: August 2025  

Purpose:
    Automates the submission of heterotrimeric and homotrimeric ANGPTL protein combinations to AlphaFoldServer.com using Playwright.
    The script reads a list of protein combinations (with copy numbers), automatically fills in the form,
    waits for the model to finish (within a specified wait time), and downloads the resulting prediction files (.zip). Output files are
    renamed for clarity and unzipped into organized folders.

Inputs:
    - A combination file listing the protein FASTA file names and copy numbers (e.g., h_A3_mature.txt:2,...)
    - Individual FASTA files in a specified sequence directory

Outputs:
    - .zip files downloaded from AlphaFold server
    - Unzipped result folders stored in a structured output directory

Dependencies:
    - Playwright (Python API)
    - Google login session (handled by first-time setup)
    - FASTA files containing only valid amino acid sequences

Usage:
    $ python ANGPTL_alphafold_submitter_playwright.py

Notes:
    - This script is configured for macOS by default. Adjust paths and settings if running on Windows.
    - Requires a one-time login to Google via browser to initiate a session.
    - Supports multiple model submissions and tracks progress per combination.
    - This script requires manual entry of the third combination entry due to playwright's inability to distinguish the third fasta 
      sequence text box. The sequence will be printed in the terminal for ease.
    - For more efficient use of following scripts, inputing the combination in terms of ANGPTL:2, ANGPTL:1, Lipase, will yield the proper 
      chain order that is hardcoded in the downstream analysis.
"""


# ---- CONFIGURATION ----
import os
import time
import shutil
import asyncio
import sys
import select
import zipfile
import platform
from zipfile import ZipFile
from datetime import datetime
from glob import glob
from playwright.sync_api import TimeoutError
from playwright.sync_api import sync_playwright
from pathlib import Path

SEQUENCE_DIR = "sequences"
DOWNLOAD_DIR = Path.home() / "Downloads" # CHANGE WHERE THE RESULTS ARE INITIALLY DOWNLOADED HERE
RESULTS_DIR = os.path.join(os.getcwd(), "alphafold_results/heterotrimer/")	# CHANGE WHERE THE RESULTS WILL BE MOVED TO HERE
COMBINATIONS_FILE = "./combinations/ht_combinations.txt"	# CHANGE THE .TXT FILE THAT CONTAINS THE COMBINATIONS IN QUESTION
ALPHAFOLD_URL = "https://alphafoldserver.com/"

# THIS IS MAC SPECIFIC. FOR WINDOWS INPUT (USER_DATA_DIR = os.path.expandvars(r"%APPDATA%\Playwright\alphafold_profile"))
USER_DATA_DIR = os.path.expanduser("~/Library/Applications Support/Playwright/alphafold_profile")

#-- MAKE SURE TO CHANGE BACK TO MAC SETTING BEFORE EXITING, THANK YOU --

SESSION_FILE = "google_state.json"
os.makedirs(DOWNLOAD_DIR, exist_ok=True)

def read_sequence(filename):
	"""
	Reads a FASTA file and returns the cleaned amino acid sequence.

	Args:
		filename (str): Path to the FASTA file containing the protein sequence.

	Returns:
		str: Cleaned amino acid sequence (uppercase, with invalid characters removed).
	"""

	path = os.path.join(SEQUENCE_DIR, filename)
	print(f"Looking for file at: {path}")
	if not os.path.exists(path):
		raise FileNotFoundError(f"File not found: {path}")

	with open(path, "r") as f:
		lines = f.readlines()
	
	sequence = "".join(line.strip().upper() for line in lines if not line.startswith(">"))
	
	valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
	cleaned = "".join([aa for aa in sequence if aa in valid_aa])
	return cleaned

	
def expand_combo(combo_line):  # USE THE COMBINATION NAME TO SAY HOW MANY OLIGOMERS EXIST OF THE PROTEIN
	"""
	Parses a line from the combination file (ht_combinations.txt) and expands it into a list of (filename, count) tuples.

	Args:
		combo_line (str): A string from the combination file (e.g., "h_A3_mature.txt:2,h_A5_mature.txtx1").

	Returns:
		expanded (list of tuples): List of (filename, int) pairs indicating the sequence file and its copy number.
	"""
	items = combo_line.strip().split(",")
	expanded = []
	for item in items:
		if ":" in item:
			fname,count = item.split(":")
			expanded.append((fname.strip(), int(count)))
		else: 
			expanded.append((item.strip(), 1))
	return expanded

# USE THE NAME OF THE COMBINATION THAT ENDS WITH ":X" TO DESIGNATE HOW MANY OLIGOMERS. THIS FUNCTION SEPARATES h_A3_mature.txt:3 into h_A3_mature.txt, 3 TO SPECIFY HOW MANY ENTITIES SHOULD BE ENTERED INTO ALPHAFOLD


def submit_sequences(page, combo,result_index):
	"""
	Submits a protein combination to AlphaFoldServer.com using Playwright.

	Args:
		page (playwright.sync_api.Page): The Playwright page object.
		combo (list): List of (filename, count) pairs.
		result_index (int): Index of the current result used for naming and logging.

	Returns:
		None
	"""
	page.goto(ALPHAFOLD_URL)

	# ADD ENTITIES AS NEEDED 
	for _ in range(len(combo) - 1):
		page.locator("button:has-text('Add entity')").click()
		time.sleep(0.5)
			
	# WAIT UNTIL ALL INPUT TEXTAREAS ARE PRESENT
	expected_inputs = len(combo)
	timeout = 10
	for _ in range(timeout * 2):
		textareas = page.locator("textarea")
		try:
			if textareas.count() >= expected_inputs:
				break
		except Exception:
			pass
		page.wait_for_timeout(500)
	else:
		raise TimeoutError(f"Timed out waiting for {expected_inputs} input boxes.")
	
	for idx, (filename, _) in enumerate(combo):
		sequence = read_sequence(filename)
		box = textareas.nth(idx)

		# TRY FILLING AUTOMATICALLY
		success = False
		for attempt in range(1):
			try:
				box.scroll_into_view_if_needed()
				box.click(force = True)
				page.wait_for_timeout(300)
				box.fill("")
				box.fill(sequence)
				success = True
				print(f"Auto-filled box #{idx + 1}")
				break
			except Exception as e:
				print(f"Attempt {attempt+1} to fill box #{idx + 1} failed: {e}")
				time.sleep(1)

		if not success:
			print(f"\n[Manual Entry Required] Could not auto-fill box #{idx + 1}.")
			print(f"-> Paste this sequence into the input box #{idx + 1} on the website:\n\n{sequence}")
			input("\nThen press ENTER here in the terminal to continue...")
		
	# SET COPY NUMBERS
	number_inputs = page.query_selector_all("input[type='number']")
	for i, (_, count) in enumerate(combo):
		if i < len(number_inputs):
			page.evaluate(
				"""(args) => {
					const element = args.element;
					const value = args.value;
					element.value = value;
					element.dispatchEvent(new Event('input', { bubbles: true }));
					element.dispatchEvent(new Event('change',{ bubbles: true }));
				}""",
				{
					"element": number_inputs[i],
					"value": count
				}
			)

	# SUBMIT THE JOB
	page.locator("button:has-text('Continue and Preview job')").click()
	page.locator("button:has-text('Confirm and submit job')").click()
	timestamp = datetime.now().strftime("%Y-%m-%d_%H:%M")

	if wait_for_completion(page, combo, timestamp):
		print("Download completed.")
	else:
		print("Download skipped due to timeout.")

# -- CAN SET A SHORTER OR LONGER WAIT TIME HERE --
def wait_for_completion(page, combo, timestamp, wait_before_click_minutes=40):
	"""
	Waits for a job to complete on AlphaFoldServer.com and downloads the resulting model.

	Args:
		page (playwright.sync_api.Page): The Playwright page object.
		combo (list): The list of (filename, count) pairs submitted.
		timestamp (str): Timestamp used for identifying the correct result card.
		wait_before_click_minutes (int): Max wait time before forcing result check.

	Returns:
		bool: True if download succeeded, False otherwise.
	"""
	
	wait_seconds = wait_before_click_minutes * 60
	check_interval = 30 
	elapsed = 0
	
	print(f"(Waiting up to {wait_before_click_minutes} minutes before opening job result...)")
	time.sleep(2)
	print("Press ENTER to open result immediately when it's ready.")

	enter_pressed = False

	while elapsed < wait_seconds:
		if sys.stdin in select.select([sys.stdin], [], [], 0)[0]:
			sys.stdin.readline()
			print("Enter pressed. Proceeding to open job result early.")
			enter_pressed = True
			break
		
		print(f"... waiting ({elapsed // 60} min elapsed)")
		time.sleep(check_interval)
		elapsed += check_interval
		
	if not enter_pressed:
		print(f"{wait_before_click_minutes} minutes elapsed. Proceeding to open job result.")
# FIND THE RESULT CARD FOR THE LATEST SUBMISSION BASED ON SUBMISSION TIME
	try:
		result_card = page.locator(f"td:has-text('{timestamp}')")
		if result_card.is_visible() and result_card.is_enabled():
			result_card.click(force = True)
			print(f"Opened job result card for {timestamp}.")
			time.sleep(5)			
			
			# FIND AND CLICK THE DOWNLOAD BUTTON
			download_link = page.locator("a[download]:has-text('Download')")
			download_link.wait_for(state="visible", timeout=10000)
			
			#DOWNLOAD THE FILE FOR THE STRUCTURE AND RENAME THEM BASED ON COMBO NAME AND MOVE THEM FROM DOWNLOAD_DIR TO RESULTS_DIR (THESE CAN BE CHANGED AT THE TOP OF THE SCRIPT)
			with page.expect_download() as download_info:
				download_link.click()
			download = download_info.value
			
			suggested_name = download.suggested_filename
			if not suggested_name.endswith(".zip"):
				print(f"WARNGING: unexpected file type: {suggested_name}")			

			combo_name = "_".join([f.replace(".txt","") + f"x{count}" for f, count in combo])
			base_name = f"{combo_name}"
			new_name = f"{base_name}.zip"
			destination = os.path.join(RESULTS_DIR, new_name)
			
			index = 1
			while os.path.exists(destination):
				new_name = f"{base_name}_{index}.zip"
				destination = os.path.join(RESULTS_DIR, new_name)
				index += 1

			os.makedirs(RESULTS_DIR, exist_ok=True)
			download.save_as(destination)

			print(f"Downloaded and moved to: {destination}")		
		
			 # UNZIP THE MOVED FILE INTO A FOLDER WITH THE SAME NAME
			try:
				unzip_dir = os.path.join(RESULTS_DIR, combo_name)
				os.makedirs(unzip_dir, exist_ok=True)
				with ZipFile(destination, 'r') as zip_ref:
                      			zip_ref.extractall(unzip_dir)
				print(f"Unzipped contents to: {unzip_dir}")
			except Exception as e:
				print(f"Failed to unzip file: {e}")
	
			try:
				back_links = page.locator("a:has-text('Back')")
				count = back_links.count()
				for i in range(count):
					candidate = back_links.nth(i)
					if candidate.is_visible() and candidate.is_enabled():
						candidate.click(force = True)
						print("Back link clicked.")
						time.sleep(3)

			except Exception as e:
				print(f"ERROR while trying to click back link: {e}")
				return True

			try:
				clear_button = page.locator("button:has-text('Clear')")
				if clear_button.is_visible() and clear_button.is_enabled():
					clear_button.click(force = True)
					print("Previous submission cleared.")
					page.wait_for_selector("textarea")
					return True
				else:
					print("WARNING: No visible and enabled 'Clear' button found.")
			except Exception as e:
				print(f"ERROR while trying to click clear button: {e}")
				return True
		else:	
			print("WARNING: Result card not visible or not enabled.")
	
	except TimeoutError:
		print("ERROR: Timeout while trying to interact with result card.")
		return False

	except Exception as e:
		print(f"ERROR during wait_for_completion: {e}")
		return False


#---- MAIN LOOP ----
if __name__ == "__main__":
	# GO THROUGH EACH LINE OF THE COMBINATION FILE AND RUN THE SCRIPT ONE LINE AT A TIME EXCEPT COMMENTED OUT LINES
	with open(COMBINATIONS_FILE) as f:
		combo_lines = [line.strip() for line in f if line.strip() and not line.strip().startswith("#")]

		# LAUNCH CHROME THROUGH PLAYWRIGHT AND ALLOW THE USER TO SEE WHAT IS HAPPENING
	with sync_playwright() as p:
		browser = p.chromium.launch(headless = False) # LEAVE AS FALSE. CAN CHANGE TO TRUE ONLY AFTER MANUAL GOOGLE LOG IN
		
		#CREATE A TEMP CONTEXT WITHOUT A SAVED SESSION
		if not os.path.exists(SESSION_FILE):
			print("[INFO] No saved session found. Starting login process...")
			temp_context = browser.new_context(accept_downloads=True)
			temp_page = temp_context.new_page()
			temp_page.goto(ALPHAFOLD_URL)
			time.sleep(2)

			#CHECK IF LOGIN IS REQUIRED
			needs_login = temp_page.locator("a:has-text('Continue with Google')").is_visible()
			if needs_login:
				print("Google login required. Please log in using the browser window.")
				input("-> After logging in, press ENTER here to save your session...")
			
				for _ in range(30):
					if "alphafoldserver.com" in temp_page.url:
						break
					print(f"Waiting for redirect back to Alphafold... (Current URL: {temp_page.url})")
					time.sleep(1)

				else:
					print(f"WARNING: Still not on Alphafold after login. Saving session anyway.")
			# SAVE THE WORKING SESSION (WHETER IT IS NEW OR ALREADY LOGGED IN)
			temp_context.storage_state(path=SESSION_FILE)
			temp_context.close()
			print("[INFO] Session saved.")

		# USE THIS SESSION FOR REAL WORK
		context = browser.new_context(
			accept_downloads=True,
			storage_state=SESSION_FILE
		)
		page = context.new_page()
		page.goto(ALPHAFOLD_URL)

		for idx, line in enumerate(combo_lines):
			combo = expand_combo(line)
			print(f"\n-> Submitting combo {idx + 1}: {combo}")	
			try:
				submit_sequences(page,combo,idx)
				time.sleep(5)
			except Exception as e:
				print(f"X Submission {idx + 1} failed: {e}")
				continue

		browser.close()
