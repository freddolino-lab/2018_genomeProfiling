#!/usr/bin/env python2 

"""
This script is to convert SAM files to FASTA files, keeping read names and read sequences, print to screen

USAGE: sam2fasta.py <input_sam_file> 

ARGS:
input_sam_file: a SAM file to be converted to FASTA file

INPUT: a SAM file
OUTPUT: converted FASTA	file, printed to standard output (screen)
"""

import sys

# global constants
## naming conventions: Section 1.4 https://samtools.github.io/hts-specs/SAMv1.pdf
QNAME_COL = 1 - 1
SEQ_COL = 10 - 1


def format_sam_line_to_fasta_line(line):
    """
    Format each line of SAM file to a FASTA file entry
    """
    fields = line.split("\t")
    cat_string = ">" + fields[QNAME_COL] + "\n" + fields[SEQ_COL]
    return(cat_string)


def main():
    """ 
    Main function: parse arguments and handles file I/O
    """
    in_filename = sys.argv[1]
    # I decided to remove the output file sys argv because
    # this can be easily redirected using bash `>`
    # out_filename = sys.argv[]

    with open(in_filename) as in_file:
        for sam_line in in_file:
            if(not sam_line.startswith("@")):
                fasta_line = format_sam_line_to_fasta_line(sam_line)
                print(fasta_line)


if __name__ == "__main__": 
    main()
