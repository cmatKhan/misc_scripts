#!/bin/bash


cat $1 | while read fastq_file run_dir; do
   fastq_path=/lts/mblab/Crypto/rnaseq_data/lts_sequence/${run_dir}_samples/${fastq_file}

    echo "moving $fastq_path"
    rsync -aHv $fastq_path /scratch/mblab/chasem/rnaseq_pipeline/scratch_sequence/unident_fastq_all

done
