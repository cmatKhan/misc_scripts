#!/bin/bash

# use samtools index to index all files in the input directory. must have sorted alignment files with suffix _sorted_aligned_reads.bam
# usage: samtools_index_all_files.sh path_to_directory

for file in $1/*_sorted_aligned_reads.bam;
do
    samtools index $file
done
