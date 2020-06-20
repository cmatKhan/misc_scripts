#!/bin/sh

#SBATCH -N 1
#SBATCH --mem=12G


# $1 sam output $2 is strandedness -- either yes, no or reverse -- $3 is the sorted_bam file $4 is the gtf $5 is the name of the output $6 is the name of the log

module load htseq

htseq-count -o $1 -s $2 -f bam -t exon $3 $4> $5 2> $6
