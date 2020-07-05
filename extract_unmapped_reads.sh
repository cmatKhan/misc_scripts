#!/bin/bash

# extract unmapped reads from a bam file
# usage: extract_unmapped_reads /path/to/sorted_aligned_reads.bam or /path/to/aligned_reads.bam

# dependency: samtools
# see http://www.novocraft.com/documentation/novoalign-2/novoalign-ngs-quick-start-tutorial/1040-2/
#     https://www.biostars.org/p/56246/#56251

samtools view -f 4 $1
