#!/usr/bin/env python3

"""
    remove the sequences identified by their header in the filter_file.txt
 
        eg) if there are chromosomes 1 to 3 in fasta in format >chr1
                                                               ATAGCA...
                                                               >chr2
                                                               ATAGCCC...
                                                               >chr2
                                                               AAAATCC...

            and the filter_file.txt has the following: chr1
                                                       chr3

            the output will be >chr2
                               ATAGCCC...

    usage: remove_scaffold_from_genome_fasta.py -g genome.fasta -f filter_file.txt

    written_by: Kamil S Jaron (https://bioinformatics.stackexchange.com/a/3940)
    included_in_misc_scripts: 5/21/20 by chase.mateusiak@gmail.com
    dependencies: Biopython, argparse
"""
from Bio import SeqIO
import argparse
import os

def main(argv):
    """
        main method
    """
    print('...parsing input')
    args = parse_args(argv)
    try:
        if not os.path.isfile(args.genome)
            raise FileNotFoundError('GenomeFileNotFound')
    except FileNotFoundError:
        print('The path to the genome file is not valid. Check and try again')
    else:
        genome_fasta_path = args.genome
    try:
        if not os.path.isfile(filter_file_path = args.filter_file


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="filter scaffold from genome fasta")
    parser.add_argument("-g", "--genome", required=True,
                        help="[REQUIRED] genome in fasta format"
    parser.add_argument("-f", "--filter_file", required=True,
                        help="[REQUIRED] file with the headers to remove (see doc string at top of file, use head, to see example/description)")
    args = parser.parse_args(argv[1:])
    return args


if __name__ == "__main__":
    main(sys.argv)
