## THE uchime files processed with this script to make them suitable for the blast
#!/bin/bash
for fasta_file in *.fasta
do
    awk '/^>/ { gsub(";", "|"); print; next } { print }' "$fasta_file" > "BLAST_$fasta_file"
done

