#!/usr/bin/env python3

import sys, os, re
from random import random, choice
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Seq import Seq

__version__ = "0.0.1"

"""
Add back fasta sequences to SURVIVOR output that were omitted due to
the length cutoff.
"""

def parse_fasta(fasta_f):
    """
    Parse a fasta sequence into a dict keyed on name.
    """
    ret = {}
    for fasta in SeqIO.parse(fasta_f, "fasta"):
        ret[fasta.name] = fasta
    return ret


def main():
    parser = ArgumentParser(description='Add back fasta sequences to SURVIVOR output that were omitted due to the length cutoff.')
    parser.add_argument('-r', '--reference_fasta', metavar='FASTA', type=str,
                        required=True, help='Path to genome assembly fasta containing reference contigs supplied to SURVIVOR.')
    parser.add_argument('-s', '--survivor_fasta', metavar='FASTA', type=str,
                        required=True, help='Path to fasta results from SURVIVOR.')
    
    
    args = parser.parse_args()

    survivor_fasta = parse_fasta(args.survivor_fasta)

    for ref_fasta in SeqIO.parse(args.reference_fasta, "fasta"):
        fasta = ref_fasta
        if ref_fasta.name in survivor_fasta:
            fasta = survivor_fasta[ref_fasta.name]        
        sys.stdout.write("%s\n" % fasta.format("fasta"))
        
    sys.stderr.write("Done!\n")
        
if __name__ == '__main__':
    main()
    exit(0)
