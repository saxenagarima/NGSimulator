# NGSimulator

NGSimulator is a tool to benchmark variant callers. It generates reads simulating sequencing errors such as insertion, deletion, and substitution. Users can define the position, type and length of variant.  

### Usage:
```
python ngsimulator.py
```

### Input:
ngsimulator takes two input files:

1) Reference genome to be sequenced in fasta format
2) A tsv file with following column:

Col1- Type of Varaiant (Str-"SNP,Ins,Del")

Col2- Position in genome (Int)

Col3- Frequency of variant in percentage (Int)

Col4- Length of the Variant (Int)

Col5- Allele (optional)

### Output:
ngsimulator makes two output files:

Simulated reads in fastq format
A text file with variant information
