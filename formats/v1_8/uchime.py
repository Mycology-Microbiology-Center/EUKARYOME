from Bio import SeqIO

def transform_header(header):
    # Split the header, remove the last element
    parts = header.split(';')[:-1]  # Removes the 'unused' or any last word
    # Define taxonomic levels with their prefixes
    levels = ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
    # Replace "." with "unclassified" and prepend levels
    transformed_parts = []
    for i, part in enumerate(parts):
        # Check if the index exists in levels, if not break the loop
        if i >= len(levels):
            break
        transformed_part = levels[i] + (part if part != '.' else 'unclassified')
        transformed_parts.append(transformed_part)
    # Ensure all taxonomic levels are represented, even if not in original header
    for j in range(i, len(levels)):
        transformed_parts.append(levels[j] + "unclassified")
    # Rejoin the parts into a single string
    return ';'.join(transformed_parts)

def process_fasta(input_file, output_file):
    with open(output_file, 'w') as out_f:
        for record in SeqIO.parse(input_file, 'fasta'):
            # Transform the header
            new_header = transform_header(record.description)
            # Write the transformed header and unwrapped sequence to the output file
            out_f.write(f">{new_header}\n{str(record.seq)}\n")
            
input_file = 'LSU.fasta'
output_file = 'uchime_LSU.fasta'
process_fasta(input_file, output_file)
