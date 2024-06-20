import os
import subprocess
from collections import defaultdict

def remove_hyphens_from_fasta(input_dir, output_dir):
    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Dictionary to track headers and check uniqueness
    headers = defaultdict(int)
    unique_headers = True
    sequence_without_header = False

    # Process each file in the input directory
    for filename in os.listdir(input_dir):
        if filename.endswith(".fasta"):
            input_filepath = os.path.join(input_dir, filename)
            temp_filepath = os.path.join(input_dir, "temp_" + filename)
            output_filepath = os.path.join(output_dir, filename)
            
            # Run seqkit commands to handle lowercase sequences
            subprocess.run(["seqkit", "seq", "-u", input_filepath, "-o", temp_filepath])
            subprocess.run(["seqkit", "seq", "-w", "0", temp_filepath, "-o", input_filepath])
            
            with open(input_filepath, 'r') as infile, open(output_filepath, 'w') as outfile:
                current_header = None
                for line in infile:
                    if line.startswith('>'):
                        if current_header is None:
                            sequence_without_header = True
                        current_header = line.strip()
                        headers[current_header] += 1
                        outfile.write(current_header + '\n')
                    else:
                        if current_header is None:
                            sequence_without_header = True
                        sequence = line.replace('-', '').strip()
                        outfile.write(sequence + '\n')
                if current_header is None:
                    sequence_without_header = True

            # Clean up temporary file
            os.remove(temp_filepath)
    
    # Check for uniqueness of headers
    if any(count > 1 for count in headers.values()):
        unique_headers = False

    return unique_headers, sequence_without_header

# Usage
input_directory = 'input_fasta_files'  # Path to the directory containing input FASTA files
output_directory = 'output_fasta_files'  # Path to the directory to save modified FASTA files
unique_headers, sequence_without_header = remove_hyphens_from_fasta(input_directory, output_directory)

# Output the status
if unique_headers:
    print("All headers are unique.")
else:
    print("There are repetitive headers.")

if sequence_without_header:
    print("There are sequences without headers.")
else:
    print("All sequences have headers.")

