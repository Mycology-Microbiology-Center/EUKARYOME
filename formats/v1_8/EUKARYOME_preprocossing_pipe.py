"""
EUKARYOME Pre-Processing Pipeline - Revised with Cleanup, Multiprocessing, Duplicate Elimination,
Uppercase Sequence Output, Final Header Standardization, and Nucleotide Filtering

This pipeline converts TXT files to FASTA format, cleans headers (retaining only the accession ID and defined taxonomic ranks), 
removes hyphens, and filters sequences to include only valid nucleotide characters. It then generates duplicate reports 
based on cleaned full headers and accession IDs. After reporting, it eliminates duplicates by comparing sequence lengths—keeping only 
the longest sequence among duplicates (with sequences converted to uppercase). Finally, headers in the deduplicated FASTA files 
are standardized to retain only the accession ID plus 7 taxonomic fields (with a trailing semicolon). A cleanup routine removes 
previously generated files and state to ensure a fresh run.
"""

import os
import json
import logging
import csv
import shutil
from collections import defaultdict
import sys
import unicodedata
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
from functools import partial
import re

# ------------------------------
# Global Configuration
# ------------------------------

STATE_FILE = "workflow_state.json"
MAX_PROCESSES = 8

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("pipeline.log", mode='a', encoding='utf-8'),
        logging.StreamHandler()
    ]
)

# ------------------------------
# Cleanup Functions
# ------------------------------

def clear_directory(directory):
    """Delete the directory if it exists and recreate it."""
    if os.path.exists(directory):
        shutil.rmtree(directory)
    os.makedirs(directory)

def clear_previous_run(directories, state_file=STATE_FILE):
    """Clear all specified directories and remove the state file if it exists."""
    for directory in directories:
        clear_directory(directory)
    if os.path.exists(state_file):
        os.remove(state_file)
        logging.info(f"Removed previous state file: {state_file}")

# ------------------------------
# Utility Functions
# ------------------------------

def normalize_to_ascii(text):
    """Normalize Unicode characters to ASCII-compatible form."""
    normalized = unicodedata.normalize('NFKD', text)
    return ''.join(c for c in normalized if not unicodedata.combining(c))

def load_state():
    """Load the workflow state from a JSON file."""
    if os.path.exists(STATE_FILE):
        try:
            with open(STATE_FILE, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception as e:
            logging.error(f"Error loading state file: {e}")
            return {}
    return {}

def save_state(state):
    """Save the workflow state to a JSON file."""
    try:
        with open(STATE_FILE, "w", encoding="utf-8") as f:
            json.dump(state, f, indent=4)
    except Exception as e:
        logging.error(f"Error saving state file: {e}")

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

# ------------------------------
# Nucleotide Filtering Function
# ------------------------------

def filter_nucleotide_sequence(seq):
    """
    Filter a nucleotide sequence by removing any characters not in the allowed set.
    Allowed nucleotides: A, C, G, T, N, R, Y, S, W, K, M, B, D, H, V.
    """
    VALID_NUCLEOTIDES = "ACGTNRYSWKMBDHV"
    seq = seq.upper()
    return re.sub(f'[^{VALID_NUCLEOTIDES}]', '', seq)

# ------------------------------
# Header Cleaning Functions
# ------------------------------

def extract_accession_id(header):
    """
    Extract the accession ID from a FASTA header.
    Expected format: >EUK1703800;Fungi;Ascomycota;...
    Returns the accession ID.
    """
    if not header.startswith('>'):
        return None
    parts = header.split(';')
    if parts:
        return parts[0].lstrip('>')
    return None

def clean_taxonomic_header(header, expected_fields=9):
    """
    Clean the header by keeping only the accession ID and the defined taxonomic ranks.
    Splits the header by semicolon and returns only the first 'expected_fields' fields.
    """
    if not header.startswith('>'):
        return header
    parts = header.split(';')
    cleaned = ';'.join(parts[:expected_fields])
    return cleaned

# ------------------------------
# Multiprocessing Modules
# ------------------------------

def convert_single_txt(file, input_dir, output_dir):
    """
    Convert a single TXT file to FASTA format.
    """
    if not file.endswith('.txt'):
        return
    state = load_state()
    base = os.path.splitext(file)[0] + ".fasta"
    output_filepath = os.path.join(output_dir, base)
    if state.get(base, {}).get("converted", False):
        logging.info(f"Skipping conversion for {file} (already converted).")
        return
    input_filepath = os.path.join(input_dir, file)
    try:
        with open(input_filepath, 'r', encoding='utf-8', errors='replace') as infile:
            content = infile.read()
        with open(output_filepath, 'w', encoding='utf-8', errors='replace') as outfile:
            outfile.write(content)
        logging.info(f"Converted {file} to {base}")
        update_state(state, base, "converted", True)
    except Exception as e:
        err = f"Conversion error for {file}: {e}"
        logging.error(err)
        update_state(state, base, "converted", False, err)

def convert_txt_to_fasta(input_dir, output_dir):
    """
    Concurrently convert all TXT files in the input directory to FASTA format.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    files = [f for f in os.listdir(input_dir) if f.endswith('.txt')]
    with ProcessPoolExecutor(max_workers=MAX_PROCESSES) as executor:
        list(executor.map(partial(convert_single_txt, input_dir=input_dir, output_dir=output_dir), files))

def preprocess_single_fasta(filename, input_dir):
    """
    Preprocess a single FASTA file by cleaning headers.
    """
    filepath = os.path.join(input_dir, filename)
    tmp_filepath = os.path.join(input_dir, f"tmp_{filename}")
    try:
        with open(filepath, 'r', encoding='utf-8', errors='replace') as infile, \
             open(tmp_filepath, 'w', encoding='utf-8', errors='replace') as outfile:
            for line in infile:
                line = line.strip()
                if line.startswith('>') or line.startswith('">'):
                    if line.startswith('">'):
                        line = line.strip('"')
                    line = line.replace('"', '')
                    line = clean_taxonomic_header(line)
                    line = normalize_to_ascii(line)
                outfile.write(line + '\n')
        os.replace(tmp_filepath, filepath)
        logging.info(f"Preprocessed {filename} to standardize headers")
    except Exception as e:
        logging.error(f"Error preprocessing {filename}: {e}")
        if os.path.exists(tmp_filepath):
            os.remove(tmp_filepath)

def preprocess_fasta_files(input_dir):
    """
    Concurrently preprocess all FASTA files in the input directory to clean headers.
    """
    files = [f for f in os.listdir(input_dir) if f.endswith(".fasta")]
    with ProcessPoolExecutor(max_workers=MAX_PROCESSES) as executor:
        list(executor.map(partial(preprocess_single_fasta, input_dir=input_dir), files))

def remove_hyphens_single_file(filename, input_dir, output_dir):
    """
    Remove hyphens from sequences in a single FASTA file.
    """
    input_filepath = os.path.join(input_dir, filename)
    output_filepath = os.path.join(output_dir, filename)
    try:
        with open(input_filepath, 'r', encoding='utf-8', errors='replace') as infile, \
             open(output_filepath, 'w', encoding='utf-8', errors='replace') as outfile:
            for line in infile:
                if not line.startswith('>'):
                    line = line.replace('-', '')
                outfile.write(line)
        logging.info(f"Removed hyphens in {filename}")
    except Exception as e:
        logging.error(f"Error processing {filename} for hyphen removal: {e}")

def process_files_remove_hyphens(input_dir, output_dir):
    """
    Concurrently process FASTA files to remove hyphens.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    files = [f for f in os.listdir(input_dir) if f.endswith(".fasta")]
    with ProcessPoolExecutor(max_workers=MAX_PROCESSES) as executor:
        list(executor.map(partial(remove_hyphens_single_file, input_dir=input_dir, output_dir=output_dir), files))

def process_file_for_elimination(filename, input_dir, output_dir):
    """
    Process a single file for duplicate elimination.
    """
    in_path = os.path.join(input_dir, filename)
    out_path = os.path.join(output_dir, filename)
    eliminate_duplicates_in_file(in_path, out_path)

def eliminate_duplicates_in_directory(input_dir, output_dir):
    """
    Concurrently eliminate duplicate sequences in each FASTA file from input_dir.
    Writes deduplicated files to output_dir.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    files = [f for f in os.listdir(input_dir) if f.endswith(".fasta")]
    with ProcessPoolExecutor(max_workers=MAX_PROCESSES) as executor:
        list(executor.map(partial(process_file_for_elimination, input_dir=input_dir, output_dir=output_dir), files))

# ------------------------------
# Duplicate Reporting Module
# ------------------------------

def generate_duplicate_reports(input_dir, report_dir):
    """
    Analyze FASTA files to identify duplicates based on cleaned full headers and accession IDs.
    Generates two Excel reports: one for duplicate full headers and one for duplicate accession IDs.
    """
    if not os.path.exists(report_dir):
        os.makedirs(report_dir)
    
    duplicate_full_headers = []  # List of [Filename, Header, Count, Line Numbers]
    duplicate_accessions = []    # List of [Filename, Accession ID, Count, Headers]
    
    for filename in os.listdir(input_dir):
        if not filename.endswith('.fasta'):
            continue
        filepath = os.path.join(input_dir, filename)
        header_counts = defaultdict(lambda: {'count': 0, 'lines': []})
        accession_counts = defaultdict(lambda: {'count': 0, 'headers': []})
        
        try:
            with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
                line_num = 0
                for line in f:
                    line_num += 1
                    line = line.strip()
                    if line.startswith('>'):
                        header = line  # Already cleaned in preprocessing
                        header_counts[header]['count'] += 1
                        header_counts[header]['lines'].append(line_num)
                        acc_id = extract_accession_id(header)
                        if acc_id:
                            accession_counts[acc_id]['count'] += 1
                            accession_counts[acc_id]['headers'].append(header)
            for header, info in header_counts.items():
                if info['count'] > 1:
                    duplicate_full_headers.append([filename, header, info['count'], ', '.join(map(str, info['lines']))])
            for acc_id, info in accession_counts.items():
                if info['count'] > 1:
                    duplicate_accessions.append([filename, acc_id, info['count'], " | ".join(info['headers'])])
        except Exception as e:
            logging.error(f"Error analyzing duplicates in {filename}: {e}")
    
    headers_report_path = os.path.join(report_dir, 'duplicate_full_headers_report.xlsx')
    if duplicate_full_headers:
        df_full = pd.DataFrame(duplicate_full_headers, columns=['Filename', 'Full Header', 'Repetitions', 'Line Numbers'])
    else:
        df_full = pd.DataFrame([['No duplicates found', '', '', '']], columns=['Filename', 'Full Header', 'Repetitions', 'Line Numbers'])
    df_full.to_excel(headers_report_path, index=False)
    logging.info(f"Duplicate full headers report generated at {headers_report_path}")
    
    accessions_report_path = os.path.join(report_dir, 'duplicate_accession_ids_report.xlsx')
    if duplicate_accessions:
        df_acc = pd.DataFrame(duplicate_accessions, columns=['Filename', 'Accession ID', 'Repetitions', 'Headers'])
    else:
        df_acc = pd.DataFrame([['No duplicates found', '', '', '']], columns=['Filename', 'Accession ID', 'Repetitions', 'Headers'])
    df_acc.to_excel(accessions_report_path, index=False)
    logging.info(f"Duplicate accession IDs report generated at {accessions_report_path}")

# ------------------------------
# Duplicate Elimination Module
# ------------------------------

def eliminate_duplicates_in_file(input_filepath, output_filepath):
    """
    Eliminate duplicate sequences from a FASTA file by keeping only the longest sequence 
    for each accession ID (and hence full header). Reads sequences, compares lengths, 
    removes hyphens, filters out invalid nucleotide characters, converts to uppercase, and writes unique entries.
    """
    sequences = {}  # key: accession id, value: (header, sequence)
    order = []      # maintain order of first appearance
    current_header = None
    current_seq = []
    
    try:
        with open(input_filepath, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    # Process previous sequence
                    if current_header is not None:
                        # Remove hyphens, filter invalid characters, and convert to uppercase
                        seq_str = filter_nucleotide_sequence("".join(current_seq).replace('-', ''))
                        acc_id = extract_accession_id(current_header)
                        if acc_id:
                            if acc_id in sequences:
                                if len(seq_str) > len(sequences[acc_id][1]):
                                    sequences[acc_id] = (current_header, seq_str)
                            else:
                                sequences[acc_id] = (current_header, seq_str)
                                order.append(acc_id)
                    current_header = line
                    current_seq = []
                else:
                    current_seq.append(line)
            # Process final sequence
            if current_header is not None:
                seq_str = filter_nucleotide_sequence("".join(current_seq).replace('-', ''))
                acc_id = extract_accession_id(current_header)
                if acc_id:
                    if acc_id in sequences:
                        if len(seq_str) > len(sequences[acc_id][1]):
                            sequences[acc_id] = (current_header, seq_str)
                    else:
                        sequences[acc_id] = (current_header, seq_str)
                        order.append(acc_id)
    except Exception as e:
        logging.error(f"Error eliminating duplicates in {input_filepath}: {e}")
        return
    
    try:
        with open(output_filepath, "w", encoding="utf-8", errors="replace") as out:
            for acc in order:
                header, seq = sequences[acc]
                out.write(f"{header}\n{seq}\n")
        logging.info(f"Eliminated duplicates in {os.path.basename(input_filepath)}; kept {len(sequences)} unique sequences.")
    except Exception as e:
        logging.error(f"Error writing deduplicated file for {input_filepath}: {e}")

# ------------------------------
# Final Header Standardization Module
# ------------------------------

def standardize_dereplicated_headers_file(filename, input_dir, output_dir, expected_fields=8):
    """
    Standardize the header of a deduplicated FASTA file so that only the accession ID 
    plus the first 'expected_fields' taxonomic fields are kept.
    For each header, split by semicolon, take the first expected_fields, and append a trailing semicolon.
    """
    input_filepath = os.path.join(input_dir, filename)
    output_filepath = os.path.join(output_dir, filename)
    try:
        with open(input_filepath, "r", encoding="utf-8", errors="replace") as infile, \
             open(output_filepath, "w", encoding="utf-8", errors="replace") as outfile:
            for line in infile:
                if line.startswith('>'):
                    header = line[1:].strip()
                    fields = header.split(';')
                    fields = fields[:expected_fields]
                    new_header = '>' + ';'.join(fields) + ';\n'
                    outfile.write(new_header)
                else:
                    outfile.write(line)
        logging.info(f"Standardized header in {filename}")
    except Exception as e:
        logging.error(f"Error standardizing header in {filename}: {e}")

def standardize_dereplicated_headers_in_directory(input_dir, output_dir, expected_fields=8):
    """
    Concurrently standardize headers for all FASTA files in the deduplicated directory.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    files = [f for f in os.listdir(input_dir) if f.endswith(".fasta")]
    with ProcessPoolExecutor(max_workers=MAX_PROCESSES) as executor:
        list(executor.map(partial(standardize_dereplicated_headers_file, input_dir=input_dir, output_dir=output_dir, expected_fields=expected_fields), files))

# ------------------------------
# Main Function
# ------------------------------

def main():
    """
    Main function orchestrating the workflow:
      1. Concurrently convert TXT files to FASTA.
      2. Concurrently preprocess FASTA files to clean headers.
      3. Concurrently remove hyphens.
      4. Generate duplicate reports based on cleaned headers and accession IDs.
      5. Concurrently eliminate duplicates by retaining only the longest sequence.
      6. Final step - Standardize headers in the deduplicated FASTA files so that only the accession ID plus 7 taxonomic fields remain.
      7. Cleanup previous generated files before running the pipeline.
    """
    # Directory configuration
    original_input_dir = "."                        # Input files (.txt or .fasta)
    converted_dir = "./converted_fasta"             # Files converted from TXT
    preprocessed_dir = "./preprocessed_fasta"       # FASTA files after header cleaning
    final_output_dir = "./final_fasta"              # FASTA files after hyphen removal
    report_dir = "./reports"                        # Directory for duplicate reports
    dedup_dir = "./deduplicated_fasta"              # Directory for deduplicated FASTA files
    standardized_dir = "./final_standardized_fasta" # Directory for final standardized deduplicated FASTA files

    # Clear previous generated files and state
    directories_to_clear = [converted_dir, preprocessed_dir, final_output_dir, report_dir, dedup_dir, standardized_dir]
    clear_previous_run(directories_to_clear, STATE_FILE)
    
    # Step 1: Convert TXT files to FASTA concurrently
    logging.info("Starting concurrent conversion from TXT to FASTA...")
    convert_txt_to_fasta(original_input_dir, converted_dir)
    
    # Step 2: Copy converted FASTA files to preprocessed_dir and clean headers concurrently
    logging.info("Copying converted files and preprocessing to standardize headers...")
    for filename in os.listdir(converted_dir):
        if filename.endswith(".fasta"):
            src = os.path.join(converted_dir, filename)
            dst = os.path.join(preprocessed_dir, filename)
            with open(src, 'r', encoding='utf-8', errors='replace') as f_in, \
                 open(dst, 'w', encoding='utf-8', errors='replace') as f_out:
                f_out.write(f_in.read())
    preprocess_fasta_files(preprocessed_dir)
    
    # Step 3: Remove hyphens concurrently
    logging.info("Removing hyphens concurrently from preprocessed files...")
    process_files_remove_hyphens(preprocessed_dir, final_output_dir)
    
    # Step 4: Generate duplicate reports based on cleaned headers
    logging.info("Generating duplicate reports based on cleaned headers...")
    generate_duplicate_reports(final_output_dir, report_dir)
    
    # Step 5: Eliminate duplicates concurrently by retaining only the longest sequence
    logging.info("Eliminating duplicates by keeping only the longest sequence...")
    eliminate_duplicates_in_directory(final_output_dir, dedup_dir)
    
    # Step 6: Final step - Standardize headers in the deduplicated FASTA files
    logging.info("Standardizing headers in deduplicated FASTA files to retain only accession ID plus 7 taxonomic fields...")
    standardize_dereplicated_headers_in_directory(dedup_dir, standardized_dir, expected_fields=8)
    
    logging.info("===== PIPELINE COMPLETE =====")
    logging.info(f"Converted files: {converted_dir}")
    logging.info(f"Preprocessed files: {preprocessed_dir}")
    logging.info(f"Final output files (hyphen removed): {final_output_dir}")
    logging.info(f"Duplicate reports available in: {report_dir}")
    logging.info(f"Deduplicated FASTA files available in: {dedup_dir}")
    logging.info(f"Final standardized FASTA files available in: {standardized_dir}")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.critical(f"Pipeline failed with error: {e}")
        sys.exit(1)
