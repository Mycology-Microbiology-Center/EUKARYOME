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
