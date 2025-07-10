#### AUTHOR: ALEXANDER DUMAS
#### LAST UPDATED: 05/30/25
#### MOST RECENT STATUS: FULLY FUNCTIONAL

#### NOTES:
# RIGHT NOW THIS IS SET FOR A MAC, SEE CONFIGURATION SECTION TO ALTER FOR WINDOWS. IF CHANGED TO WINDOWS CONFIGURATION, CHANGE BACK TO MAC WHEN DONE. 
# I HAVE YET TO TRY IF WE CAN RUN THIS ON MULTIPLE MACHINES AT THE SAME TIME TO SPEED UP PROGRESS. AGAIN I DONT KNOW IF THIS WILL WORK BUT THE CHANGES BELOW SHOULD HELP IT WORK.
# IF RUNNING ON MULTIPLE MACHINES MAKE SURE LOGGING IN WITH DIFFERENT GOOGLE ACCOUNTS (HAVE DIFFERENT SESSION_FILES PER USER) / NOT USE THE SAME RESULTS_DIR / ONE USER ON WINDOWS AND THE OTHER ON MAC SO THE USER_DATA_DIR IS DIFFERENT PER USER / ALSO USE INDEPENDENT COMBINATION_FILE PER USER
# WOULD BE BEST IF THE SECOND USER COPIES ALL NEEDED FILES TO THEIR SECTION OF THE RDSS DRIVE AS TO BOTH NOT BE PULLING FROM A_DUMAS


# --------- REMINDER THIS IS SET FOR HOMOTRIMER TESTING. SEE OTHER SCRIPT FOR HETEROTRIMER SUBMISSIONS ---------

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

#---- CONFIGURE ----
SEQUENCE_DIR = "sequences"
DOWNLOAD_DIR = Path.home() / "Downloads" # CHANGE WHERE THE RESULTS ARE INITIALLY DOWNLOADED HERE
RESULTS_DIR = os.path.join(os.getcwd(), "alphafold_results")	# CHANGE WHERE THE RESULTS WILL BE MOVED TO HERE
COMBINATIONS_FILE = "combinations.txt"	# CHANGE THE .TXT FILE THAT CONTAINS THE COMBINATIONS IN QUESTION
ALPHAFOLD_URL = "https://alphafoldserver.com/"

# THIS IS MAC SPECIFIC. FOR WINDOWS INPUT (USER_DATA_DIR = os.path.expandvars(r"%APPDATA%\Playwright\alphafold_profile"))
USER_DATA_DIR = os.path.expanduser("~/Library/Applications Support/Playwright/alphafold_profile")
#== MAKE SURE TO CHANGE BACK TO MAC SETTING BEFORE EXITING, THANK YOU ==

SESSION_FILE = "google_state.json"
os.makedirs(DOWNLOAD_DIR, exist_ok=True)

#---- READ SEQUENCES ----
def read_sequence(filename):
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

# CLEAN UP THE FASTA SEQUENCES TO YIELD STRICT AMINO ACID SEQUENCE
	

#---- RENAME ALPHAFOLD DOWNLOADS ---- NOT CURRENTLY USING THIS FUNCTION ---- RUNNING THIS WITHIN WAIT_FOR_COMPLETION()----
#def rename_latest_download( combo):
#	time.sleep(5)
#	
#	zip_files = glob(os.path.join(DOWNLOAD_DIR, "fold_*.zip"))
#	if not zip_files:
#		print("WARNING: No downloaded ZIP files found.")
#		return
#
#	zip_files.sort(key=os.path.getmtime, reverse=True)
#	latest_file = zip_files[0]
#
#	combo_name = "_".join([f.replace(".txt", "") + f"x{count}" for f, count in combo])
#	new_name = f"{combo_name}.zip"
#
#	os.makedirs(RESULTS_DIR, exist_ok=True)
#	destination = os.path.join(RESULTS_DIR, new_name)
#
#	try:
#		shutil.move(latest_file, destination)
#		print(f"Moved and renamed to: {destination}")
#	except Exception as e:
#		print(f"Failed to move/rename file: {e}")
#	
#	# UNZIP THE MOVED FILE INTO A FOLDER WITH THE SAME NAME
#	try:
#		unzip_dir = os.path.join(RESULTS_DIR, combo_name)
#		os.makedirs(unzip_dir, exist_ok=True)
#		with ZipFile(destination, 'r') as zip_ref:
#			zip_ref.extractall(unzip_dir)
#		print(f"Unzipped contents to: {unzip_dir}")
#	except Exception as e:
#		print(f"Failed to unzip file: {e}")



	# CHANGE THE NAME NAMES OF THE ALPHAFOLD DOWNLOADS TO MIMIC THAT OF THE COMBO USED AND ALSO PLACE IN RESPECTIVE RESULTS DIRECTORY
	
# ---- EXPAND THE COMBINATIONS ----
def expand_combo(combo_line):  # USE THE COMBINATION NAME TO SAY HOW MANY OLIGOMERS EXIST OF THE PROTEIN
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

#---- SUBMIT SEQUENCES ----
def submit_sequences(page, combo,result_index):
	page.goto(ALPHAFOLD_URL)

	# ADD ENTITIES AS NEEDED 
	for _ in range(len(combo) - 1):
		page.locator("button:has-text('Add entity')").click()
		time.sleep(0.5)
	
	# FILL FASTA SEQUENCES
	textareas = page.locator("textarea[placeholder*='Input']")
	for i, (filename, _) in enumerate(combo):
		sequence = read_sequence(filename)
		textareas.nth(i).fill("")
		textareas.nth(i).fill(sequence)

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

#---- WAIT FOR COMPLETION OF THE JOB ----
# CAN SET A SHORTER OR LONGER WAIT TIME HERE
# ALLOW THE USER TO HIT ENTER IF THE JOB IS DONE BEFORE THE ALLOTED TIME
def wait_for_completion(page, combo, timestamp, wait_before_click_minutes=20):
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

				print("WARNING: No visible and enabled 'Back' link found.")
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

#def wait_for_download_file(download_dir, timeout=30):
#	"""Wait until a new .zip file appears in the download directory."""
#	start_time = time.time()
#	while time.time() - start_time < timeout:
#		zip_files = sorted(download_dir.glob("*.zip"), key=lambda f: f.stat().st_mtime, reverse=True)
#		if zip_files:
#		# CHECK IF FILE IS STILL DOWNLOADING
#			latest_file = zip_files[0]
#			size1 = latest_file.stat().st_size
#			time.sleep(1)
#			size2 = latest_file.stat().st_size
#			if size1 == size2 and size1 > 0:
#				return True
#		time.sleep(1)
#	return False


#def run_combinations():
#	with open(COMBINATIONS_FILE) as f:
#		combo_lines = [line for line in f if line.strip()][:5]
#
#	with sync_playwright() as p:
#		browser = p.chromium.launch(headless = False)
#		
#		#HANDLE SESSION: LOAD IF EXISTS, ELSE PROMPT LOGIN
#		if os.path.exists(SESSION_FILE):
#			print("Using saved Google session.")
#			context = browser.new_context(
#				accept_downloads=True,
#				storage_state=SESSION_FILE
#			)
#		
#		else:
#			print("No saved session found. Logging in manually...")
#			context = browser.new_context(accept_downloads=True)
#
#		page = context.new_page()
#		page.goto(ALPHAFOLD_URL)
#
#		if not os.path.exists(SESSION_FILE):
#			print("Please log into Google in the browser window.")
#			input("Press Enter here when you're done logging in...")
#			context.storage_state(path=SESSION_FILE)
#			print("Session saved.")
#
#		for idx, line in enumerate(combo_lines):
#			combo = expand_combo(line)
#			print(f"\n-> Submitting combo {idx + 1}: {combo}")
#			try:
#				submit_sequences(page,combo,idx)
#				time.sleep(5)
#			except Exception as e:
#				print(f"X Submission {idx + 1} failed: {e}")
#		
#		browser.close()


#---- MAIN LOOP ----
if __name__ == "__main__":
	# GO THROUGH EACH LINE OF THE COMBINATION FILE AND RUN THE SCRIPT ONE LINE AT A TIME EXCEPT COMMENTED OUT LINES
	with open(COMBINATIONS_FILE) as f:
		combo_lines = [line.strip() for line in f if line.strip() and not line.strip().startswith("#")]
	# GO LINE BY LINE
	for idx, line in enumerate(combo_lines):
		combo = expand_combo(line)

		# LAUNCH CHROME THROUGH PLAYWRIGHT AND ALLOW THE USER TO SEE WHAT IS HAPPENING
		with sync_playwright() as p:
			browser = p.chromium.launch(headless = False) # LEAVE AS FALSE. CAN CHANGE TO TRUE ONLY AFTER MANUAL GOOGLE LOG IN

              		#HANDLE SESSION: LOAD IF EXISTS, ELSE PROMPT LOGIN
			if os.path.exists(SESSION_FILE):
				print("Using saved Google session.")
				context = browser.new_context(
					accept_downloads=True,
					storage_state=SESSION_FILE
				)

			else:
				print("No saved session found. Logging in manually...")
				context = browser.new_context(accept_downloads=True)

			page = context.new_page()
			page.goto(ALPHAFOLD_URL)
		
			if not os.path.exists(SESSION_FILE):
				print("Please log into Google in the browser window.")
				input("Press Enter here when you're done logging in...")
				context.storage_state(path=SESSION_FILE)
				print("Session saved.")
			
			print(f"\n-> Submitting combo {idx + 1}: {combo}")	
			try:
				submit_sequences(page,combo,idx)
				time.sleep(5)
			except Exception as e:
				print(f"X Submission {idx + 1} failed: {e}")
				continue

	browser.close()
