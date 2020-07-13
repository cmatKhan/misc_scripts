#!/bin/bash

# written by chase.mateusiak@gmail.com
# date: 07/05/2020

# loop over all files in align_count_results and create a look_up.txt sheet for
# unindexed bam files (see wustl htcf 'getting started' about array jobs with slurm)

# INPUT $1: name of the lookup file (_lookup.txt will be appended), and the file
# will be placed in rnaseq_pipeline/job_scripts

# TODO: check that the timestamp for the .bai and .bam are the same somehow --
# if the .bam file was created after hte .bai, then it probably won't work (it turns out)

# loop over all .bam files in align_count_results/*/align
for align_file in /scratch/mblab/${USER}/rnaseq_pipeline/align_count_results/*/align/*.bam; do
    # if the .bam.bai file does not exist, then
    if [ ! -f ${align_file}.bai ]; then
        # add the file to the look_up list
        echo ${align_file} >> /scratch/mblab/${USER}/rnaseq_pipeline/job_scripts/${1}_lookup.txt
    fi
done

echo "Look up file at: /scratch/mblab/${USER}/rnaseq_pipeline/job_scripts/${1}_lookup.txt"
