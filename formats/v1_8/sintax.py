##THE SCRIPT IS DESIGNED TO CHANGE THE BLAST FORMATTED FILE INTO THE SINTAX ONE####
#NOTE: the General release should be provided as an input#
import os
import re
from multiprocessing import Pool, cpu_count

def convert_blast_to_sintax(input_fasta):
    output_fasta = f'SINTAX_EUK_{os.path.basename(input_fasta)}'
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                # Remove '>' at the start and strip whitespace
                line = line[1:].strip()
                
                # Split the line into parts
                parts = line.split(';')
                
                # Extract the ID and taxonomy
                seq_id = parts[0]
                taxonomy = parts[1:]
                
                # Format the taxonomy for SINTAX
                sintax_taxonomy = []
                rank_dict = {
                    'k__': 'd:',  # Change kingdom to domain for SINTAX
                    'p__': 'p:',
                    'c__': 'c:',
                    'o__': 'o:',
                    'f__': 'f:',
                    'g__': 'g:',
                    's__': 's:'
                }
                
                for taxon in taxonomy:
                    for key, value in rank_dict.items():
                        if taxon.startswith(key):
                            taxon_name = taxon.split('__')[1]
                            if taxon_name != 'unclassified':
                                sintax_taxonomy.append(f"{value}{taxon_name}")
                            break
                
                # Join the SINTAX formatted taxonomy
                sintax_header = f">{seq_id};tax={','.join(sintax_taxonomy)};\n"
                outfile.write(sintax_header)
            else:
                outfile.write(line)

def process_fasta_files(fasta_file):
    convert_blast_to_sintax(fasta_file)

def main():
    # Get list of all FASTA files in the current directory
    fasta_files = [f for f in os.listdir('.') if f.endswith('.fasta')]
    
    # Use up to 8 CPUs for parallel processing
    num_cpus = min(8, cpu_count())
    with Pool(num_cpus) as pool:
        pool.map(process_fasta_files, fasta_files)

if __name__ == "__main__":
    main()
