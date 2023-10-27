import pandas as pd

def excel_to_new_fasta_format(excel_file_path, fasta_file_path):
    # Step 1: Read the Excel file
    df = pd.read_excel(excel_file_path)
    
    # Step 2: Create FASTA headers in the new format
    fasta_content = ""
    for index, row in df.iterrows():
        header = f"{row['acc_new']}    k__{row['kgd']};p__{row['ph']};c__{row['c']};o__{row['o']};f__{row['f']};g__{row['g']};s__{row['s']}"
        sequence = row['seq']
        
        # Step 3: Append to FASTA content
        fasta_content += f">{header}\n{sequence}\n"
        
    # Step 4: Write to new FASTA file
    with open(fasta_file_path, 'w') as f:
        f.write(fasta_content)

# Usage
excel_file_path = 'path/to/your/excel/file.xlsx'
fasta_file_path = 'path/to/your/new/output/fasta/file.fasta'
excel_to_new_fasta_format(excel_file_path, fasta_file_path)


def excel_to_qiime2_tsv(excel_file_path, tsv_file_path):
    """
    Convert an Excel file to a TSV file compatible with QIIME2's FeatureData[Taxonomy] format.
    
    Parameters:
    excel_file_path (str): Path to the input Excel file containing taxonomy and sequence information.
    tsv_file_path (str): Path to the output TSV file compatible with QIIME2.
    
    Returns:
    None
    """
    # Step 1: Read the Excel file
    df = pd.read_excel(excel_file_path)
    
    # Step 2: Create DataFrame for TSV format
    feature_ids = df['acc_new']
    taxons = df.apply(lambda row: f"k__{row['kgd']};p__{row['ph']};c__{row['c']};o__{row['o']};f__{row['f']};g__{row['g']};s__{row['s']}", axis=1)
    
    tsv_df = pd.DataFrame({
        'Feature ID': feature_ids,
        'Taxon': taxons
    })
    
    # Step 3: Write to TSV file
    tsv_df.to_csv(tsv_file_path, sep='\t', index=False)

# Usage example (Note: The actual file paths should be provided)
excel_file_path = 'path/to/your/excel/file.xlsx'
tsv_file_path = 'path/to/your/new/output/tsv/file.tsv'
excel_to_qiime2_tsv(excel_file_path, tsv_file_path)
