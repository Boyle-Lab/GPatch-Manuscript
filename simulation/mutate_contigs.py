#!/usr/bin/env python3

import sys, os, re
from random import random, choice
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Seq import Seq

__version__ = "0.0.1"

"""
Introduce mutations/errors at a given, fixed rate into fasta records.
Assumes equal chance of switching between all nucleotides and all other
nucleotides.
"""


def main():
    parser = ArgumentParser(description='Introduce mutations/errors at a given, fixed rate into fasta records. Assumes equal chance of switching between all nucleotides and all other nucleotides.')
    parser.add_argument('-f', '--fasta', metavar='FASTA', type=str,
                        required=True, help='Path to genome assembly fasta containing chromosomes to break into pseudocontigs.')
    parser.add_argument('-r', '--rate', metavar='x', type=float,
                        required=False, default=0.01,
                        help='Error/mutation rate: the per-base chance of switching from one nucleotide to any other nucleotide. Default = 0.01 (1%)')
    
    args = parser.parse_args()

    for fasta in SeqIO.parse(args.fasta, "fasta"):
        seq = list(fasta.seq)
        # Traverse the sequence, randomly choosing to mutate or
        # not mutate, at the given rate.
        for i, s in enumerate(seq):
            val = random()
            if val < args.rate:
                # choose a random nucleotide that's different.
                seq[i] = choice([x for x in "ACTG" if x != s.upper()])
        # Print the mutated sequence in fasta format
        fasta.seq = "".join(seq)
        sys.stdout.write("%s\n" % fasta.format("fasta"))
        
    sys.stderr.write("Done!\n")
        
if __name__ == '__main__':
    main()
    exit(0)
