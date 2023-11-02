import pandas as pd

def clean_sequence(sequence):
    """Retain only IUPAC standard characters for DNA sequences and convert to uppercase."""
    iupac_chars = set('ACGTRYKMSWBDHVN')
    return ''.join(c for c in sequence.upper() if c in iupac_chars)

def generate_blast_fasta(excel_file_path, fasta_file_path, chunk_size=5000):
    # Convert Excel to temporary CSV for chunked reading
    temp_csv_path = 'temp_chunked_blast_file.csv'
    pd.read_excel(excel_file_path).to_csv(temp_csv_path, index=False)
    
    # Initialize output FASTA file
    with open(fasta_file_path, 'w') as f_fasta:
        pass
    
    # Read and process CSV file in chunks
    for chunk in pd.read_csv(temp_csv_path, chunksize=chunk_size):
        chunk['sequence'] = chunk['sequence'].apply(clean_sequence)
        
        fasta_data = []
        for idx, row in chunk.iterrows():
            # Skip the first column (presumably 'acc_new' is the second column)
            row_data = row.iloc[2:]  # Assuming 'acc_new' is at index 1
            # Replace "." with "unclassified" in the header elements
            header_elements = [
                f"{k}__{'unclassified' if v == '.' else v}" 
                for k, v in row_data.items() if k != 'sequence'
            ]
            header = f"{row['acc_new']}|{'|'.join(header_elements)}"
            fasta_data.append(f">{header}\n{row['sequence']}")
        
        # Write this chunk's data to FASTA file
        with open(fasta_file_path, 'a') as f_fasta:
            f_fasta.write('\n'.join(fasta_data) + '\n')

    # Clean up the temporary CSV file
    import os
    os.remove(temp_csv_path)
    
# File Paths
excel_file_path = 'path/to/your/excel/file.xlsx'
fasta_file_path = 'path/to/your/new/output/fasta/file.fasta'

# Generate FASTA File
generate_blast_fasta(excel_file_path, fasta_file_path)
