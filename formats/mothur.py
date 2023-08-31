import pandas as pd

def excel_to_fasta_and_tax(excel_file_path, fasta_file_path, tax_file_path):
    # Step 1: Read the Excel file
    df = pd.read_excel(excel_file_path)
    
    # Step 2: Generate FASTA and Taxonomic content
    fasta_content = ""
    tax_content = ""
    for index, row in df.iterrows():
        fasta_header = f"{row['acc_new']}    k__{row['kgd']};p__{row['ph']};c__{row['c']};o__{row['o']};f__{row['f']};g__{row['g']};s__{row['s']}"
        tax_line = f"{row['acc_new']}\tk__{row['kgd']};p__{row['ph']};c__{row['c']};o__{row['o']};f__{row['f']};g__{row['g']};s__{row['s']};"
        sequence = row['seq']
        
        fasta_content += f">{fasta_header}\n{sequence}\n"
        tax_content += f"{tax_line}\n"
        
    # Step 3: Write to new files
    with open(fasta_file_path, 'w') as f:
        f.write(fasta_content)
    
    with open(tax_file_path, 'w') as f:
        f.write(tax_content)

# Usage
excel_file_path = 'path/to/your/excel/file.xlsx'
fasta_file_path = 'path/to/your/new/output/fasta/file.fasta'
tax_file_path = 'path/to/your/new/output/tax/file.tax'
excel_to_fasta_and_tax(excel_file_path, fasta_file_path, tax_file_path)
