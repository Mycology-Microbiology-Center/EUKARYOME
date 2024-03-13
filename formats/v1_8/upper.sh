##The Script designed to tackle the problem of the sequences with lower case present in the original files. 

#!/bin/bash

# Create the input directory if it doesn't exist
mkdir -p input

# Loop over all .fasta files in the current directory
for fasta_file in *.fasta
do
    # Run the seqkit commands
    seqkit seq -u "$fasta_file" > up_"$fasta_file"
    seqkit seq -w 0 up_"$fasta_file" > "$fasta_file"

    # Move the final files to the input directory
    mv "$fasta_file" input/
done 
