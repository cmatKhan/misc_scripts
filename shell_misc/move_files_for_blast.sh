#!/bin/bash

cat ../count_files_for_blast.txt | while read line; do
    rsync -aHv $line experiments/unmapped_blast/
done
