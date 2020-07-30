#!/bin/sh

#module load blast-plus

# $1 query fasta $2 output name
# expects to be launched from non_coding_crypto with a directory called KN99_blastdb 

# for non tabular, remove -outfmt 6
#blastn -query $1  -db KN99_blastdb/KN99_blastdb -out $2 -outfmt "6 qseqid sseqid length pident evalue bitscore sstrand sstart send"

# for non tabular output
blastn -query $1 -db KN99_blastdb/KN99_blastdb -out $2 
