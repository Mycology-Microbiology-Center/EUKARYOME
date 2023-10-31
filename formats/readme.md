For running the scripts, you need to first install the Pandas library using the below command:

``` bash
pip install pandas
```

and later run the codes by the command `python script.py'. 
  
The Input excel file should have this headers format as an input for generating fasta files out of it for different software:
- `Acc_new`: This represents the sequence ID.
- `sequence`: This represents the sequence data.
- Taxonomic columns: `Kingdom`, `Phylum`, `class`, `order`, `family`, `genus`, `species`  as short forms in headers: k, p, c, o, f, g, s
