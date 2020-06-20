#!/bin/bash

# from a directory with multi runs, loop through run directories and extract unmapped reads as sam files for each
# result will be <library_name>_unmapped.sam for each library in run
# required that novoalign was run with -r All (or at least -r and some threshold of multimaps) see novoalign docs

# $1 is the directory $2 is the gtf

for run_dir in $1/*;
do
    for bam_file in $run_dir/*_sorted_aligned_reads.bam;
    do
        unmapped_suffix=_unmapped.bam
        unmapped_path=${bam_file}/_sorted_aligned_reads.bam/${unmapped_suffix}
        # see http://www.novocraft.com/userfiles/file/Novocraft.pdf last page of sam format where it says ZS:Z:R not present for unique alignments
        echo "working on ${unmapped_path}"
        #samtools view $bam_file | grep ZS:Z:R | samtools view -T $2 -S -b > $unmapped_path
    done
done
