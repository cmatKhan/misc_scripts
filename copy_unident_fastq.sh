#!/bin/bash

for file in /scratch/mblab/chasem/misc_scripts/unident_by_run/*million*; do
    cat $file | while read line; do
        if [[! -e scratch_sequence/unident_fastq/$line]]; then
             echo "$line does not exist in scratch yet"
            #rsync -aHv $line scratch_sequence/unident_fastq/
        fi
    done
done
