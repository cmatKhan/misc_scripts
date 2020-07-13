#!/bin/bash

# subsample bam file. Note: for saturation curve for protein_coding_count
# this assumes that the bamfile has already been filtered for only protein protein_coding_count
# unique alignment

# INPUT $1: total number of rows in bamfile
# INPUT $2: the bamfile to subset

input_bam=$1
no_path_input_bam=${input_bam##*/}
basename_no_ext=${no_path_input_bam%_sorted_aligned_reads_with_annote.bam}

for power in {11..24};do
  for seed in {1..10};do

    # set power given power of 2 to variable n
    power_of_two=$(echo "2^${power}"|bc)
 
   # get the fraction of the input bam file that the given power_of_two corresponds to
    frac=$( samtools idxstats $input_bam | cut -f3 | awk -v "n=$power_of_two" 'BEGIN {total=0} {total += $1} END {frac=n/total; if (frac > 1) {print 1} else {print frac}}' )

    # drop leading zero from $frac and replace with $seed
    if (( $(echo "$frac < 1" |bc -l) )); then
      nozero=$(echo $frac | sed 's/^0*//')
      seed_frac=${seed}${nozero}
    else
      seed_frac=1.1
  fi
  echo "${seed}, ${frac}, ${seed_frac}"
      
    # use sametools with the given seed (see outer loop) and the fraction of the library to subset thereads
    #echo "working on 2 to the ${power} reads, replicate ${seed}"
    #samtools view -bs $seed.$frac $input_bam > ${basename_no_ext}_${seed}_${power}_protein_coding_subset.bam
  done
done
