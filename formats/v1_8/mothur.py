from Bio import SeqIO
import csv
import os

def generate_mothur_files(input_fasta, output_fasta, output_tax):
    with open(output_fasta, 'w') as fasta_file, open(output_tax, 'w', newline='') as tax_file:
        # Prepare the .tax file writer
        tax_writer = csv.writer(tax_file, delimiter='\t')
        
        for record in SeqIO.parse(input_fasta, 'fasta'):
            # Extract the ID and taxonomy from the record description
            parts = record.description.split(';')
            feature_id = parts[0]  # The sequence ID
            taxonomy = ';'.join(parts[1:])  # The taxonomy string, joined with semicolons
            
            # Write the taxonomy to the .tax file
            tax_writer.writerow([feature_id, taxonomy])
            
            # Write the sequence to the FASTA file with the modified ID
            fasta_file.write(f">{feature_id}\n{str(record.seq)}\n")

# List of FASTA filenames in the directory to process
fasta_files = ['ITS.fasta', 'LSU.fasta', 'SSU.fasta', 'longread.fasta']

for fasta_file in fasta_files:
    # Construct the output filenames based on the input filename
    base_name = os.path.splitext(fasta_file)[0]
    output_fasta = f'{base_name}_mothur.fasta'
    output_tax = f'{base_name}_mothur.tax'
    
    # Generate the Mothur-compatible files
    generate_mothur_files(fasta_file, output_fasta, output_tax)
    print(f"Processed {fasta_file}: generated {output_fasta} and {output_tax}")
