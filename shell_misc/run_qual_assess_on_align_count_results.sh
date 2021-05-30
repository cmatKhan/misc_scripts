#/bin/bash

# intended to be run from rnaseq_pipeline.
# usage: run_qual_assess_on_align_count_results.sh align_count_results query_sheet.csv

for run_directory in $1/*; do
    quality_assess_1.py -r ${run_directory} -o ${run_directory} -cc -qs $2 --interactive
done
