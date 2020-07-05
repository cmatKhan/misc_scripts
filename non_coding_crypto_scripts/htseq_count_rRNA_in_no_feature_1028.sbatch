#!/bin/sh

#SBATCH -N 1
#SBATCH --mem=12G


module load htseq

htseq-count -f bam -t exon rnaseq_pipeline/run_1028_no_feature/run_1028_AGACTGA_no_feature_sorted.bam non_coding_crypto/KN99_ribo.gtf > non_coding_crypto/run_1028_no_feature_rRNA_count.tsv 2> non_coding_crypto/run_1028_no_feature_rRNA_count.log
