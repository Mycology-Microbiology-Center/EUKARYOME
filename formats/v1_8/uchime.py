from Bio import SeqIO

def transform_header(header):
    # Split the header on semicolons and remove the last element, typically "unused"
    parts = header.split(';')
    # Initialize transformed_parts with the ID part, which doesn't change
    transformed_parts = [parts[0]]
    # Define the prefixes for taxonomic ranks
    taxonomic_ranks = ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
    # Process each part of the header, replacing '.' with 'unclassified', starting from the second element
    for rank, part in zip(taxonomic_ranks, parts[1:]):
        if part == '.' or part == 'unused':
            transformed_parts.append(rank + 'unclassified')
        else:
            transformed_parts.append(rank + part)
    # If there are fewer parts than ranks, fill the remaining ranks with 'unclassified'
    for i in range(len(parts) - 1, len(taxonomic_ranks)):
        transformed_parts.append(taxonomic_ranks[i] + 'unclassified')

    return ';'.join(transformed_parts)

def process_fasta(input_file, output_file):
    with open(output_file, 'w') as out_f:
        for record in SeqIO.parse(input_file, 'fasta'):
            # Transform the header using the first space as a delimiter between the ID and the rest of the header
            new_header = transform_header(record.description)
            # Write the transformed header and unwrapped sequence to the output file
            out_f.write(f">{new_header}\n{str(record.seq)}\n")

# Specify the paths to  input and output FASTA files
input_file = 'input.fasta'
output_file = 'uchime_input.fasta'
process_fasta(input_file, output_file)
