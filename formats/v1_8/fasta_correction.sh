import os
import subprocess
from collections import defaultdict
import pandas as pd

def remove_hyphens_and_duplicates(input_dir, output_dir, report_dir):
    # Ensure the output and report directories exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if not os.path.exists(report_dir):
        os.makedirs(report_dir)

    # Dictionary to track headers and check uniqueness
    headers = defaultdict(lambda: {'count': 0, 'lines': [], 'sequences': []})
    unique_headers = True
    sequence_without_header = False

    # Process each file in the input directory
    for filename in os.listdir(input_dir):
        if filename.endswith(".fasta"):
            input_filepath = os.path.join(input_dir, filename)
            temp_filepath = os.path.join(input_dir, "temp_" + filename)
            output_filepath = os.path.join(output_dir, filename)
            
            # Run seqkit commands to handle lowercase sequences
            subprocess.run(["seqkit", "seq", "-u", input_filepath, "-o", temp_filepath])
            subprocess.run(["seqkit", "seq", "-w", "0", temp_filepath, "-o", input_filepath])
            
            with open(input_filepath, 'r') as infile, open(output_filepath, 'w') as outfile:
                current_header = None
                line_number = 0
                for line in infile:
                    line_number += 1
                    if line.startswith('>'):
                        current_header = line.strip()
                        headers[(filename, current_header)]['count'] += 1
                        headers[(filename, current_header)]['lines'].append(line_number)
                        if current_header in headers[(filename, current_header)]['sequences']:
                            continue
                        headers[(filename, current_header)]['sequences'].append('')
                    else:
                        if current_header is None:
                            sequence_without_header = True
                            with open(os.path.join(report_dir, 'sequences_without_headers.fasta'), 'a') as seq_no_header_file:
                                seq_no_header_file.write(f"{filename}: {line}")
                        else:
                            sequence = line.replace('-', '').strip()
                            if sequence in headers[(filename, current_header)]['sequences']:
                                continue
                            headers[(filename, current_header)]['sequences'][-1] += sequence
                            outfile.write(current_header + '\n' + sequence + '\n')

            # Clean up temporary file
            os.remove(temp_filepath)
    
    # Create a DataFrame for the report
    repetitive_headers_report = []
    for (filename, header), info in headers.items():
        if info['count'] > 1:
            unique_headers = False
            repetitive_headers_report.append({
                'Filename': filename,
                'Header': header,
                'Repetitions': info['count'],
                'Line Numbers': ', '.join(map(str, info['lines']))
            })

    # Write the report to an Excel file
    report_filepath = os.path.join(report_dir, 'repetitive_headers_report.xlsx')
    if repetitive_headers_report:
        df = pd.DataFrame(repetitive_headers_report)
        df.to_excel(report_filepath, index=False)

    return unique_headers, sequence_without_header

# Usage
input_directory = './original/duplicate_removed'  # Path to the directory containing input FASTA files
output_directory = './filtered_fasta'  # Path to the directory to save modified FASTA files
report_directory = './filtered_fasta'  # Path to the directory to save report files

unique_headers, sequence_without_header = remove_hyphens_and_duplicates(input_directory, output_directory, report_directory)

# Output the status
if unique_headers:
    print("All headers are unique.")
else:
    print(f"There are repetitive headers. Check {os.path.join(report_directory, 'repetitive_headers_report.xlsx')} for details.")

if sequence_without_header:
    print(f"There are sequences without headers. Check {os.path.join(report_directory, 'sequences_without_headers.fasta')} for details.")
else:
    print("All sequences have headers.")
