##### duplicate remove and clean fasta file script #####
from Bio import SeqIO
import os
import subprocess
from collections import OrderedDict

def remove_duplicates(input_fasta, output_fasta):
    unique_sequences = OrderedDict()
    
    # Parse the input FASTA file and keep only one copy of each unique tax identifier
    for record in SeqIO.parse(input_fasta, 'fasta'):
        tax_id = record.description.split(';')[0]
        if tax_id not in unique_sequences:
            unique_sequences[tax_id] = record
    
    # Write unique sequences to the output FASTA file
    with open(output_fasta, 'w') as fasta_file:
        for record in unique_sequences.values():
            SeqIO.write(record, fasta_file, 'fasta')
    
    # Use seqkit to ensure proper formatting (sequence in one line)
    temp_filepath = 'temp.fasta'
    subprocess.run(["seqkit", "seq", "-u", output_fasta, "-o", temp_filepath])
    subprocess.run(["seqkit", "seq", "-w", "0", temp_filepath, "-o", output_fasta])
    os.remove(temp_filepath)

# Example usage
input_fasta = './original/longread.fasta'
output_fasta = './filtered_fasta/longread.fasta'
remove_duplicates(input_fasta, output_fasta)
print(f"Processed {input_fasta}: generated {output_fasta}")
