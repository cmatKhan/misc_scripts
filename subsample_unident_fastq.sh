#!/bin/bash

# subsample fastq files. launch from rnaseq_pipeline dir

ml seqtk

for file in scratch_sequence/unident_fastq/*.fastq.gz; do
    sample_basename=$(basename $file)
    subsample_path = scratch_sequence/unident_subsample/$sample_basename
    if [[! -e $subsample_path]]; then
        seqtk sample $file 50000 > $subsample_path
    fi
  done
