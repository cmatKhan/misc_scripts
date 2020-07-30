#!/bin/bash


# usage: total_cnag_counts.sh /path/to/some_library_read_count.tsv  <-- read counts totaled by htseq counts

cat $1 | grep CNAG | cut -f2 | paste -sd+ | bc
