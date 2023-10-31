import pandas as pd

import pandas as pd

def clean_sequence(sequence):
    """Convert A, C, T, G to uppercase without altering other characters."""
    return ''.join(c.upper() if c.upper() in {'A', 'C', 'T', 'G'} else c for c in sequence)

def format_sequence(seq, line_length=80):
    """Break down long sequence lines into shorter segments."""
    return '\n'.join(seq[i:i+line_length] for i in range(0, len(seq), line_length))

def generate_files(excel_file_path, fasta_file_path, tsv_file_path):
    """Generate FASTA and TSV files from an Excel file."""
    # Read Excel file
    df = pd.read_excel(excel_file_path)
    
    # Clean sequences
    df['sequence'] = df['sequence'].apply(clean_sequence)
    
    # Generate FASTA file
    with open(fasta_file_path, 'w') as f_fasta:
        for idx, row in df.iterrows():
            unique_header = f"{row['acc_new']}_{idx}"  # Append index to make the header unique
            formatted_seq = format_sequence(row['sequence'])
            f_fasta.write(f">{unique_header}\n{formatted_seq}\n")
    
    # Generate TSV file
    with open(tsv_file_path, 'w') as f_tsv:
        f_tsv.write("Feature ID\tTaxon\n")
        for idx, row in df.iterrows():
            unique_header = f"{row['acc_new']}_{idx}"  # Append index to make the header unique
            taxon = ';'.join([f"{k}__{('unclassified' if v == '.' else v)}" for k, v in row[['k', 'p', 'c', 'o', 'f', 'g', 's']].items()])
            f_tsv.write(f"{unique_header}\t{taxon}\n")

# File paths
excel_file_path = '/media/ali/HDD1/reference/EUK_ITS_chim.xlsx'
fasta_file_path = '/media/ali/HDD1/reference/QIIME2_new/EUK_ITS_chim.fasta'
tsv_file_path = '/media/ali/HDD1/reference/QIIME2_new/EUK_ITS_chim.tsv'

# Generate files
generate_files(excel_file_path, fasta_file_path, tsv_file_path)


def format_sequence(seq, line_length=80):
    """
    Break down long sequence lines into shorter segments.
    """
    return '\n'.join(seq[i:i+line_length] for i in range(0, len(seq), line_length))

def generate_files(excel_file_path, fasta_file_path, tsv_file_path):
    """
    Generate FASTA and TSV files from an Excel file.
    """
    # Read Excel file
    df = pd.read_excel(excel_file_path)
    
    # Clean sequences
    df['sequence'] = df['sequence'].apply(clean_sequence)
    
    # Generate FASTA file
    with open(fasta_file_path, 'w') as f_fasta:
        for idx, row in df.iterrows():
            formatted_seq = format_sequence(row['sequence'])
            f_fasta.write(f">{row['acc_new']}\n{formatted_seq}\n")
    
    # Generate TSV file
    with open(tsv_file_path, 'w') as f_tsv:
        f_tsv.write("Feature ID\tTaxon\n")
        for idx, row in df.iterrows():
            taxon = ';'.join([f"{k}__{('unclassified' if v == '.' else v)}" for k, v in row[['k', 'p', 'c', 'o', 'f', 'g', 's']].items()])
            f_tsv.write(f"{row['acc_new']}\t{taxon}\n")

# File paths
excel_file_path = 'path/to/your/excel/file.xlsx'
fasta_file_path = 'path/to/save/fasta/file.fasta'
tsv_file_path = 'path/to/save/tsv/file.tsv'

# Generate files
generate_files(excel_file_path, fasta_file_path, tsv_file_path)
