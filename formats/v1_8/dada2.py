import os
from multiprocessing import Pool

# Function to handle DADA2 Assign Taxonomy format conversion
def convert_to_dada2_format(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                # Remove '>' and split the header
                header_parts = line[1:].strip().split(';')
                
                # Initialize taxonomy levels list
                taxonomy = []
                
                # Process each part after the ID (first element)
                for part in header_parts[1:]:
                    # Handle parts with prefixes (like k__, p__, etc.)
                    if '__' in part:
                        taxonomy.append(part.split('__')[-1])
                    else:
                        taxonomy.append(part.strip())
                
                # Ensure exactly 6 levels by either truncating or padding with 'unclassified'
                if len(taxonomy) > 6:
                    taxonomy = taxonomy[:6]
                while len(taxonomy) < 6:
                    taxonomy.append('unclassified')
                
                # Create new DADA2 format header
                new_header = '>' + ';'.join(taxonomy) + ';'
                outfile.write(new_header + '\n')
            else:
                outfile.write(line)

# Function to handle DADA2 Assign Species format conversion
def convert_fasta_headers(input_fasta, output_fasta):
    with open(input_fasta, "r") as input_handle, open(output_fasta, "w") as output_handle:
        for line in input_handle:
            line = line.strip()
            if line.startswith(">"):
                # Split the header into parts
                header_parts = line[1:].split(";")
                
                # Extract the ID (first part before the first semicolon)
                seq_id = header_parts[0].split()[0]
                
                # Initialize genus and species as empty
                genus = ""
                species = ""
                
                # Loop through header parts to identify genus and species levels
                for part in header_parts:
                    if part.startswith("g__") and "unclassified" not in part:
                        genus = part.split("__")[1]
                    if part.startswith("s__") and "unclassified" not in part:
                        species = part.split("__")[1]
                
                # Construct the new header in the desired format
                new_header = f">{seq_id}"
                if genus:
                    new_header += f" {genus}"
                if species:
                    new_header += f" {species}"
                
                output_handle.write(f"{new_header}\n")
            else:
                output_handle.write(f"{line}\n")

# Wrapper function to handle both conversions for a single file
def process_file(file_path, output_dir):
    base_name = os.path.basename(file_path)
    # Generate the two output file names
    dada2_taxonomy_output = os.path.join(output_dir, f"{base_name}_dada2_taxonomy.fasta")
    dada2_species_output = os.path.join(output_dir, f"{base_name}_dada2_species.fasta")
    
    # Process both formats
    convert_to_dada2_format(file_path, dada2_taxonomy_output)
    convert_fasta_headers(file_path, dada2_species_output)

# Function to run multiprocessing for all files in a folder
def process_all_files_in_folder(input_folder, output_folder, num_processes=4):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    fasta_files = [os.path.join(input_folder, f) for f in os.listdir(input_folder) if f.endswith('.fasta')]

    # Use multiprocessing pool to process files in parallel
    with Pool(num_processes) as pool:
        pool.starmap(process_file, [(file, output_folder) for file in fasta_files])

# Usage Example
if __name__ == "__main__":
    input_folder = "input_fasta_files"  # Path to the folder containing the input FASTA files
    output_folder = "output_fasta_files"  # Path to the folder where output files will be stored
    process_all_files_in_folder(input_folder, output_folder, num_processes=4)
