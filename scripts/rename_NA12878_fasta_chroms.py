#!/usr/bin/env python3

import sys, os, re
from random import random, choice
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Seq import Seq
import pysam

__version__ = "0.0.1"

"""
Given a fasta with non-chromosome IDs and a SAM/BAM alignment
to a reference genome with chromosome names, rename contigs with
the chromosome belonging to their primary alignment.
"""

def parse_bam(bam_f):
    """
    Parse a bam file into a dict, keyed on read name.
    """
    ret = {}
    bam = pysam.AlignmentFile(bam_f, "rb")
    for read in bam:
        if not read.is_supplementary and not read.is_secondary:
            ret[read.query_name] = read
    return ret

def main():
    parser = ArgumentParser(description='Given a fasta with non-chromosome IDs and a SAM/BAM alignment to a reference genome with chromosome names, rename contigs with the chromosome belonging to their primary alignment.')
    parser.add_argument('-f', '--fasta', metavar='FASTA', type=str,
                        required=True, help='Path to genome assembly fasta with non-chromosome names.')
    parser.add_argument('-b', '--bam', metavar='SAM/BAM', type=str,
                        required=True, help='Path to SAM/BAM alignment of the input fasta to a reference genome with chromosome names.')
    
    args = parser.parse_args()

    bam = parse_bam(args.bam)
    
    for fasta in SeqIO.parse(args.fasta, "fasta"):
        #sys.stderr.write("%s\n" % (fasta.name))
        if fasta.name in bam:
            #sys.stderr.write("%s\n" % (bam[fasta.name].reference_name))
            fasta.id = bam[fasta.name].reference_name
        else:
            sys.stderr.write("WARNING: %s not found in BAM. Not renaming!\n" % (fasta.name))
        sys.stdout.write("%s\n" % fasta.format("fasta"))
            
    sys.stderr.write("Done!\n")
        
if __name__ == '__main__':
    main()
    exit(0)
