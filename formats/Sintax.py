import pandas as pd

def excel_to_sintax_fasta(excel_file_path, fasta_file_path):
    # Step 1: Read the Excel file
    df = pd.read_excel(excel_file_path)
    
    # Step 2: Create FASTA headers in the new SINTAX format
    fasta_content = ""
    for index, row in df.iterrows():
        header = f"{row['acc_new']};tax=d:{row['kgd']},p:{row['ph']},c:{row['c']},o:{row['o']},f:{row['f']},g:{row['g']};"
        sequence = row['seq']
        
        # Step 3: Append to FASTA content
        fasta_content += f">{header}\n{sequence}\n"
        
    # Step 4: Write to new FASTA file
    with open(fasta_file_path, 'w') as f:
        f.write(fasta_content)

# Usage
excel_file_path = 'path/to/your/excel/file.xlsx'
fasta_file_path = 'path/to/your/new/sintax/output/fasta/file.fasta'
excel_to_sintax_fasta(excel_file_path, fasta_file_path)
