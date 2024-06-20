##THE SCRIPT IS DESIGNED TO CHANGE THE BLAST FORMATTED FILE INTO THE SINTAX ONE####
import re

def convert_blast_to_sintax(input_fasta, output_fasta):
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                # Remove '>' at the start
                line = line[1:]
                
                # Split the line into parts
                parts = line.split('|')
                
                # Extract the ID and taxonomy
                seq_id = parts[0]
                taxonomy = parts[1:]
                
                # Format the taxonomy for SINTAX
                sintax_taxonomy = []
                rank_dict = {
                    'k__': 'd:',
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
                            taxon = taxon.replace(key, value)
                            if 'unclassified' not in taxon:
                                sintax_taxonomy.append(taxon)
                
                # Join the SINTAX formatted taxonomy
                sintax_header = f">{seq_id};tax=" + ",".join(sintax_taxonomy) + ";\n"
                outfile.write(sintax_header)
            else:
                outfile.write(line)

# Usage
input_fasta = 'BLAST_SSU.fasta'
output_fasta = 'SINTAX_SSU.fasta'
convert_blast_to_sintax(input_fasta, output_fasta)
