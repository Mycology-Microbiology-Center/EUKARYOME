##The script is changing the uchime format files into the Qiime2 
from Bio import SeqIO
import csv
import os

def generate_qiime2_files(input_fasta, output_fasta, output_tsv):
    with open(output_fasta, 'w') as fasta_file, open(output_tsv, 'w', newline='') as tsv_file:
        # Prepare the TSV writer
        tsv_writer = csv.writer(tsv_file, delimiter='\t')
        # Write the header of the TSV file
        tsv_writer.writerow(["Feature ID", "Taxonomy"])
        
        for record in SeqIO.parse(input_fasta, 'fasta'):
            # Extract the ID and taxonomy from the record description
            parts = record.description.split(';')
            feature_id = parts[0]  # The sequence ID
            taxonomy = ';'.join(parts[1:])  # The taxonomy string
            
            # Write the taxonomy to the TSV file
            tsv_writer.writerow([feature_id, taxonomy])
            
            # Write the sequence to the FASTA file with the modified ID
            fasta_file.write(f">{feature_id}\n{str(record.seq)}\n")

# List of FASTA filenames in the directory to process
fasta_files = ['ITS.fasta', 'LSU.fasta', 'SSU.fasta', 'longread.fasta']

for fasta_file in fasta_files:
    # Construct the output filenames based on the input filename
    base_name = os.path.splitext(fasta_file)[0]
    output_fasta = f'{base_name}_qiime2.fasta'
    output_tsv = f'{base_name}_qiime2.tsv'
    
    # Generate the QIIME2-compatible files
    generate_qiime2_files(fasta_file, output_fasta, output_tsv)
    print(f"Processed {fasta_file}: generated {output_fasta} and {output_tsv}")
