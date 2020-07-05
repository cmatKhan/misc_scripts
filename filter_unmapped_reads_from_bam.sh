#!/bin/bash

ml samtools

for f in experiments/unmapped_blast/*
do
   samtools view -f4 $f > ${f/_sorted_aligned_reads.bam/_umapped.bam}
done
