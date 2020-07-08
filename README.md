# NGSimulator
##Usage:
python simseq.py
##Input:
Simseq takes two input files:

Reference genome to be sequenced in fasta format
A tsv file with following column:
Col1- Type of Varaiant (Str-"SNP,Ins,Del")
Col2- Position in genome (Int)
Col3- Frequency of variant in percentage (Int)
Col4- Length of the Variant (Int)
Col5- Allele (optional)

##Output:
Simseq makes two output files:

Simulated reads in fastq format
A text file with variant information
