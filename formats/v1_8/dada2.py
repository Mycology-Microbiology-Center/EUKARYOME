import os
import shutil

def convert_to_dada2_format(input_file, output_file):
    """
    Convert taxonomy database headers to DADA2 format, handling various input formats.
    It will discard the remaining levels after detecting 'unclassified' in any level.
    
    Parameters:
    input_file (str): Path to input FASTA file
    output_file (str): Path to output FASTA file
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                # Remove '>' and split the header
                header_parts = line[1:].strip().split(';')
                
                # Initialize taxonomy levels list
                taxonomy = []

                # Process each part of the taxonomy, skipping the ID part (first element)
                for part in header_parts[1:]:
                    if '__' in part:
                        cleaned_part = part.split('__')[-1].lstrip('_')  # Remove leading underscores and prefix
                    else:
                        cleaned_part = part.strip().lstrip('_')  # Clean the part without prefixes
                    
                    # Stop processing further levels if 'unclassified' is found
                    if 'unclassified' in cleaned_part.lower():
                        break
                    
                    taxonomy.append(cleaned_part)

                # After filtering out 'unclassified', join the valid parts
                if len(taxonomy) == 0:
                    taxonomy.append('unclassified')  # In case no classification left
                
                # Join valid taxonomy levels
                final_taxonomy = ';'.join(taxonomy)

                # Write the new header in DADA2 format
                new_header = f">{final_taxonomy};"
                outfile.write(new_header + '\n')
            else:
                # Write sequence lines as they are
                outfile.write(line)

def process_fasta_files(input_dir, output_dir):
    """
    Process all FASTA files in a directory, converting each to DADA2 format
    and saving the results in the specified output directory.
    
    Parameters:
    input_dir (str): Path to the directory containing input FASTA files
    output_dir (str): Path to the directory where output files will be saved
    """
    # Check if output directory exists; if it does, delete and recreate it
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)
    
    # Loop through all files in the input directory
    for filename in os.listdir(input_dir):
        if filename.endswith(".fasta"):
            input_file = os.path.join(input_dir, filename)
            output_file = os.path.join(output_dir, f"dada2_{filename}")
            
            # Convert the FASTA file to DADA2 format
            convert_to_dada2_format(input_file, output_file)
            print(f"Processed {filename} and saved to {output_file}")

# Set input and output directory paths
input_directory = "."  # Change this to your input directory
output_directory = "./output"  # Change this to your desired output directory

# Process all FASTA files in the input directory
process_fasta_files(input_directory, output_directory)
