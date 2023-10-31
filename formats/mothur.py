import pandas as pd

def clean_sequence(sequence):
    """Retain only IUPAC standard characters for DNA sequences and convert to uppercase."""
    iupac_chars = set('ACGTRYKMSWBDHVN')
    return ''.join(c.upper() if c.upper() in iupac_chars else c for c in sequence)

def format_sequence(seq, line_length=80):
    """Break down long sequence lines into shorter segments."""
    return '\n'.join(seq[i:i+line_length] for i in range(0, len(seq), line_length))

def generate_mothur_files(excel_file_path, fasta_file_path, tax_file_path, chunk_size=5000):
    # Convert Excel to temporary CSV for chunked reading
    temp_csv_path = 'temp_chunked_file.csv'
    pd.read_excel(excel_file_path).to_csv(temp_csv_path, index=False)
    
    # Initialize output files
    with open(fasta_file_path, 'w') as f_fasta, open(tax_file_path, 'w') as f_tax:
        pass
    
    # Read and process CSV file in chunks
    for chunk in pd.read_csv(temp_csv_path, chunksize=chunk_size):
        # Clean sequences
        chunk['sequence'] = chunk['sequence'].apply(clean_sequence)
        
        # Generate FASTA and tax data for this chunk
        fasta_data = []
        tax_data = []
        for idx, row in chunk.iterrows():
            formatted_seq = format_sequence(row['sequence'])
            fasta_data.append(f">{row['acc_new']}_{idx}\n{formatted_seq}")
            taxon = '\t'.join([f"{row['acc_new']}_{idx}", ';'.join([f"{k}__{('unclassified' if v == '.' else v)}" 
                                                               for k, v in row[['k', 'p', 'c', 'o', 'f', 'g', 's']].items()]) + ';'])
            tax_data.append(taxon)
        
        # Write this chunk's data to FASTA and tax files
        with open(fasta_file_path, 'a') as f_fasta, open(tax_file_path, 'a') as f_tax:
            f_fasta.write('\n'.join(fasta_data) + '\n')
            f_tax.write('\n'.join(tax_data) + '\n')
        
        # Write this chunk's data to FASTA and tax files
        with open(fasta_file_path, 'a') as f_fasta, open(tax_file_path, 'a') as f_tax:
            f_fasta.write('\n'.join(fasta_data) + '\n')
            f_tax.write('\n'.join(tax_data) + '\n')

# Usage
excel_file_path = 'path/to/your/excel/file.xlsx'
fasta_file_path = 'path/to/your/new/output/fasta/file.fasta'
tax_file_path = 'path/to/your/new/output/tax/file.tax'
excel_to_fasta_and_tax(excel_file_path, fasta_file_path, tax_file_path)
