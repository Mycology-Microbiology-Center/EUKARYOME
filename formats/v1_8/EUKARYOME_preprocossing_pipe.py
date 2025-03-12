"""
EUKARYOME pre-processing Pipeline

This script processes FASTA files through multiple stages to ensure data quality and consistency.
The pipeline includes the following stages:

1. Conversion: Converts .txt files to .fasta format using 'latin-1' encoding to handle unusual characters.
2. Preprocessing: Fixes common formatting issues in FASTA files, such as removing quotation marks from headers.
3. Duplicate Removal: Removes duplicate sequences based on taxonomic identifiers and sequence strings, and the longest one will remain.
4. Hyphen Removal and Header Check: Removes hyphens from sequences, cleans headers, and checks for duplicate headers.

Outputs:
- Converted FASTA files in the 'converted_fasta' directory.
- FASTA files with duplicates were removed from the 'duplicate_removed' directory.
- Final processed FASTA files in the 'filtered_fasta' directory.
- Reports on repetitive headers and sequences without headers in the 'filtered_fasta' directory.

The script uses 'latin-1' encoding throughout to handle non-ASCII characters properly and logs all operations to 'pipeline.log'.
"""
import os
import json
import logging
import subprocess
from collections import OrderedDict, defaultdict
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from multiprocessing import Pool
import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='latin-1', errors='replace')
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='latin-1', errors='replace')

# Setup logging to file and console.
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("pipeline.log", mode='a', encoding='latin-1'),
        logging.StreamHandler()
    ]
)

STATE_FILE = "workflow_state.json"

def load_state():
    """Load the workflow state from a JSON file."""
    if os.path.exists(STATE_FILE):
        with open(STATE_FILE, "r", encoding="latin-1") as f:
            return json.load(f)
    return {}

def save_state(state):
    """Save the workflow state to a JSON file."""
    with open(STATE_FILE, "w", encoding="latin-1") as f:
        json.dump(state, f, indent=4)

def update_state(state, filename, stage, success, error_msg=""):
    """Update the state dictionary for a given file and stage."""
    if filename not in state:
        state[filename] = {"converted": False, "duplicates": False, "filtered": False, "errors": {}}
    state[filename][stage] = success
    if not success:
        state[filename]["errors"][stage] = error_msg
    else:
        state[filename]["errors"].pop(stage, None)
    save_state(state)

def convert_txt_to_fasta(input_dir, output_dir):
    """
    Convert .txt files in the input directory to .fasta files in the output directory.
    Uses 'latin-1' encoding to ensure all unusual characters are rendered.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    state = load_state()
    for filename in os.listdir(input_dir):
        if filename.endswith('.txt'):
            base = os.path.splitext(filename)[0] + ".fasta"
            output_filepath = os.path.join(output_dir, base)
            # Skip if already converted successfully.
            if state.get(base, {}).get("converted", False):
                logging.info(f"Skipping conversion for {filename} (already converted).")
                continue
            input_filepath = os.path.join(input_dir, filename)
            try:
                with open(input_filepath, 'r', encoding='latin-1', errors='replace') as infile, \
                     open(output_filepath, 'w', encoding='latin-1', errors='replace') as outfile:
                    outfile.write(infile.read())
                logging.info(f"Converted {filename} to {base}")
                update_state(state, base, "converted", True)
            except Exception as e:
                err = f"Conversion error for {filename}: {e}"
                logging.error(err)
                update_state(state, base, "converted", False, err)

def preprocess_fasta_files(input_dir):
    """
    Preprocess FASTA files to fix common formatting issues like quoted headers
    """
    for filename in os.listdir(input_dir):
        if filename.endswith(".fasta"):
            filepath = os.path.join(input_dir, filename)
            tmp_filepath = os.path.join(input_dir, f"tmp_{filename}")
            
            try:
                with open(filepath, 'r', encoding='latin-1', errors='replace') as infile, \
                     open(tmp_filepath, 'w', encoding='latin-1', errors='replace') as outfile:
                    for line in infile:
                        line = line.strip()
                        # Fix quoted FASTA headers
                        if line.startswith('">'):
                            line = line.strip('"')
                        elif line.startswith('>'):
                            line = line.replace('"', '')
                        outfile.write(line + '\n')
                
                # Replace original with corrected file
                os.replace(tmp_filepath, filepath)
                logging.info(f"Preprocessed {filename} to fix formatting issues")
            except Exception as e:
                logging.error(f"Error preprocessing {filename}: {e}")
                if os.path.exists(tmp_filepath):
                    os.remove(tmp_filepath)

def remove_duplicates(input_fasta, output_fasta):
    """
    Remove duplicate sequences from a FASTA file based on the taxonomic identifier (first field of the header).
    When duplicates are found, retains only the longest sequence.
    Cleans header lines by removing double quotes and hyphens.
    Uses 'latin-1' encoding throughout and handles non-ASCII characters properly.
    """
    unique_sequences = {}  # Will hold {tax_id: (header, sequence)} with longest sequence for each tax_id
    sequence_order = []    # Maintains the order of unique taxonomic IDs for output
    
    try:
        # Read the file manually to ensure full control over encoding
        current_header = None
        current_seq = ""
        with open(input_fasta, "r", encoding="latin-1", errors='replace') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                    
                if line.startswith('>') or line.startswith('">'):
                    # Process previous sequence if any
                    if current_header is not None and current_seq:
                        # Clean sequence (remove hyphens)
                        clean_seq = current_seq.replace('-', '')
                        # Extract taxonomic ID from the cleaned header
                        cleaned_header = current_header.strip('"')
                        tax_id = cleaned_header.split(';')[0].lstrip('>')
                        
                        # Check if we already have this tax_id
                        if tax_id in unique_sequences:
                            # Compare sequence lengths and keep the longer one
                            existing_seq = unique_sequences[tax_id][1]
                            if len(clean_seq) > len(existing_seq):
                                # This sequence is longer, replace the existing one
                                unique_sequences[tax_id] = (cleaned_header, clean_seq)
                                logging.info(f"Replaced shorter sequence for {tax_id} with a longer one ({len(existing_seq)} -> {len(clean_seq)} bases)")
                        else:
                            # First time seeing this tax_id
                            unique_sequences[tax_id] = (cleaned_header, clean_seq)
                            sequence_order.append(tax_id)
                    
                    # Start new sequence, clean header by removing quotes
                    if line.startswith('">'):
                        line = line.strip('"')
                    current_header = line.replace('"', '')
                    current_seq = ""
                else:
                    # Append to current sequence
                    current_seq += line
            
            # Don't forget the last sequence
            if current_header is not None and current_seq:
                clean_seq = current_seq.replace('-', '')
                tax_id = current_header.split(';')[0].lstrip('>')
                
                if tax_id in unique_sequences:
                    existing_seq = unique_sequences[tax_id][1]
                    if len(clean_seq) > len(existing_seq):
                        unique_sequences[tax_id] = (current_header, clean_seq)
                        logging.info(f"Replaced shorter sequence for {tax_id} with a longer one ({len(existing_seq)} -> {len(clean_seq)} bases)")
                else:
                    unique_sequences[tax_id] = (current_header, clean_seq)
                    sequence_order.append(tax_id)
        
    except Exception as e:
        logging.error(f"Error reading {input_fasta}: {str(e)}")
        raise Exception(f"Error processing {input_fasta}: {str(e)}")

    # Check if any sequences were found
    if not unique_sequences:
        logging.warning(f"No valid sequences found in {input_fasta}!")
        # Create an empty file but mark as successful to continue pipeline
        with open(output_fasta, 'w', encoding="latin-1", errors='replace') as out_f:
            pass
        return False
    
    # Write sequences directly to file without using BioPython
    try:
        with open(output_fasta, 'w', encoding="latin-1", errors='replace') as out_f:
            for tax_id in sequence_order:
                header, sequence = unique_sequences[tax_id]
                out_f.write(f"{header}\n{sequence}\n")
    except Exception as e:
        logging.error(f"Error writing to {output_fasta}: {str(e)}")
        raise Exception(f"Error writing to {output_fasta}: {str(e)}")

    # Use seqkit with careful subprocess handling
    try:
        temp_filepath = output_fasta + ".temp"
        # Use binary mode for subprocess to avoid encoding issues
        subprocess.run(["seqkit", "seq", "-u", output_fasta, "-o", temp_filepath], 
                      check=True)
        subprocess.run(["seqkit", "seq", "-w", "0", temp_filepath, "-o", output_fasta], 
                      check=True)
        if os.path.exists(temp_filepath):
            os.remove(temp_filepath)
    except subprocess.SubprocessError as e:
        logging.warning(f"Seqkit processing failed for {output_fasta}: {e}")
    
    return True  # Return True if processing was successful

def process_file_for_duplicates(file_info):
    """
    Wrapper function for multiprocessing that processes duplicate removal for a single FASTA file.
    """
    input_fasta, output_dir = file_info
    base = os.path.basename(input_fasta)
    output_fasta = os.path.join(output_dir, base)
    state = load_state()
    
    # Skip if already processed for duplicates.
    if state.get(base, {}).get("duplicates", False):
        logging.info(f"Skipping duplicate removal for {base} (already processed).")
        return
    
    # Check if input file is empty
    if os.path.getsize(input_fasta) == 0:
        logging.warning(f"Input file {input_fasta} is empty! Skipping duplicate removal.")
        update_state(state, base, "duplicates", False, f"Input file {input_fasta} is empty")
        return
        
    try:
        success = remove_duplicates(input_fasta, output_fasta)
        if success:
            logging.info(f"Processed duplicates for {input_fasta} -> {output_fasta}")
            update_state(state, base, "duplicates", True)
        else:
            msg = f"No sequences found in {input_fasta}"
            logging.warning(msg)
            update_state(state, base, "duplicates", False, msg)
    except Exception as e:
        err = f"Duplicate removal error for {input_fasta}: {e}"
        logging.error(err)
        update_state(state, base, "duplicates", False, err)

def process_duplicates_in_directory(input_dir, output_dir, num_processes=4):
    """
    Process duplicate removal for all FASTA files in the input directory using multiprocessing.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    fasta_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith('.fasta')]
    file_info_list = [(file, output_dir) for file in fasta_files]
    with Pool(processes=num_processes) as pool:
        pool.map(process_file_for_duplicates, file_info_list)

def remove_hyphens_and_check_headers(input_dir, output_dir, report_dir):
    """
    Remove hyphens from sequences and check for duplicate headers in FASTA files.
    Also cleans headers by removing double quotes.
    Generates a report of repetitive headers and flags sequences without headers.
    Uses 'latin-1' encoding to handle unusual characters.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if not os.path.exists(report_dir):
        os.makedirs(report_dir)
    headers = defaultdict(lambda: {'count': 0, 'lines': [], 'sequences': []})
    unique_headers = True
    sequence_without_header = False
    state = load_state()
    empty_files = []

    for filename in os.listdir(input_dir):
        if filename.endswith(".fasta"):
            # Skip if already processed
            if state.get(filename, {}).get("filtered", False):
                logging.info(f"Skipping header processing for {filename} (already processed).")
                continue

            input_filepath = os.path.join(input_dir, filename)
            temp_filepath = os.path.join(input_dir, "temp_" + filename)
            output_filepath = os.path.join(output_dir, filename)
            
            # Check if input file is empty
            if os.path.getsize(input_filepath) == 0:
                logging.warning(f"Input file {input_filepath} is empty! Creating empty output file.")
                with open(output_filepath, 'w') as _:
                    pass  # Create empty file
                empty_files.append(filename)
                update_state(state, filename, "filtered", False, f"Input file {input_filepath} is empty")
                continue
            
            try:
                subprocess.run(["seqkit", "seq", "-u", input_filepath, "-o", temp_filepath])
                subprocess.run(["seqkit", "seq", "-w", "0", temp_filepath, "-o", input_filepath])
            except Exception as e:
                err = f"Seqkit formatting error for {filename}: {e}"
                logging.error(err)
                update_state(state, filename, "filtered", False, err)
                continue

            sequences_found = False
            try:
                with open(input_filepath, 'r', encoding="latin-1", errors='replace') as infile, \
                     open(output_filepath, 'w', encoding="latin-1", errors='replace') as outfile:
                    current_header = None
                    line_number = 0
                    for line in infile:
                        line_number += 1
                        if line.startswith('>') or line.startswith('">'):
                            # Properly handle headers that might start with a quote
                            line = line.strip()
                            if line.startswith('">'):
                                header = line.strip('"')
                            else:
                                header = line.replace('"', '')
                                
                            current_header = header
                            headers[(filename, header)]['count'] += 1
                            headers[(filename, header)]['lines'].append(line_number)
                            if header in headers[(filename, header)]['sequences']:
                                continue
                            headers[(filename, header)]['sequences'].append('')
                        else:
                            if current_header is None:
                                sequence_without_header = True
                                with open(os.path.join(report_dir, 'sequences_without_headers.fasta'), 'a', encoding="latin-1", errors='replace') as seq_no_header_file:
                                    seq_no_header_file.write(f"{filename}: {line}")
                            else:
                                sequence = line.replace('-', '').strip()
                                if sequence in headers[(filename, current_header)]['sequences']:
                                    continue
                                headers[(filename, current_header)]['sequences'][-1] += sequence
                                outfile.write(current_header + '\n' + sequence + '\n')
                                sequences_found = True
                
                os.remove(temp_filepath)
                
                if not sequences_found:
                    logging.warning(f"No valid sequences found in {filename} after header processing")
                    empty_files.append(filename)
                    update_state(state, filename, "filtered", False, "No valid sequences found after header processing")
                else:
                    update_state(state, filename, "filtered", True)
                    
            except Exception as e:
                err = f"Error in header processing for {filename}: {e}"
                logging.error(err)
                update_state(state, filename, "filtered", False, err)
                continue

    # Report empty files
    if empty_files:
        logging.warning(f"The following files have no sequences: {', '.join(empty_files)}")

    repetitive_headers_report = []
    for (filename, header), info in headers.items():
        if info['count'] > 1:
            unique_headers = False
            repetitive_headers_report.append({
                'Filename': filename,
                'Header': header,
                'Repetitions': info['count'],
                'Line Numbers': ', '.join(map(str, info['lines']))
            })

    report_filepath = os.path.join(report_dir, 'repetitive_headers_report.xlsx')
    if repetitive_headers_report:
        df = pd.DataFrame(repetitive_headers_report)
        df.to_excel(report_filepath, index=False)

    return unique_headers, sequence_without_header, empty_files

def main():
    # Directory configuration
    original_input_dir = "."                        # Original input files (.txt or .fasta)
    converted_dir = "./converted_fasta"             # FASTA files converted from .txt files
    duplicate_removed_dir = "./duplicate_removed"   # Files after duplicate sequence removal
    final_output_dir = "./filtered_fasta"           # Final processed files after hyphen and header checks
    report_dir = "./filtered_fasta"                 # Directory to save report files

    # Load workflow state
    state = load_state()

    # Stage 1: Convert .txt files to .fasta using latin-1 encoding.
    logging.info("Starting conversion from .txt to .fasta...")
    convert_txt_to_fasta(original_input_dir, converted_dir)
    
    # New stage: Preprocess FASTA files to fix common format issues
    logging.info("Preprocessing FASTA files to fix formatting issues...")
    preprocess_fasta_files(converted_dir)

    # Stage 2: Remove duplicate sequences.
    logging.info("Starting duplicate removal (keeping longest sequence)...")
    process_duplicates_in_directory(converted_dir, duplicate_removed_dir, num_processes=8)

    # Stage 3: Remove hyphens, clean headers, and check duplicate headers.
    logging.info("Starting header and hyphen processing...")
    unique_headers, sequence_without_header, empty_files = remove_hyphens_and_check_headers(duplicate_removed_dir, final_output_dir, report_dir)

    # Report final statuses.
    if unique_headers:
        logging.info("All headers are unique.")
    else:
        logging.warning(f"Some files have repetitive headers. See {os.path.join(report_dir, 'repetitive_headers_report.xlsx')}")
    if sequence_without_header:
        logging.warning(f"Some sequences are missing headers. See {os.path.join(report_dir, 'sequences_without_headers.fasta')}")
    else:
        logging.info("All sequences have headers.")
        
    if empty_files:
        logging.warning(f"The following files contain no sequences after processing: {', '.join(empty_files)}")

    # Report files that failed any stage.
    state = load_state()  # Reload state to ensure it's updated.
    error_files = {fname: info["errors"] for fname, info in state.items() if info.get("errors")}
    if error_files:
        logging.error("Files with errors detected:")
        for fname, errors in error_files.items():
            logging.error(f"{fname}: {errors}")
    else:
        logging.info("All files processed successfully.")

if __name__ == "__main__":
    main()
