#!/bin/bash

# convert sam file to bam (specifically if you get the error "no SQ lines in header" when trying samtools -bS view input_sam.sam > output_bam.bam)
# first input is genome (.fasta) second is the sam file third is the name or path of the bam output

samtools view -bS -T $1 $2 > $3
