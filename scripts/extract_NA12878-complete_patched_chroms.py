#!/usr/bin/env python3

import sys, os, re
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Seq import Seq

__version__ = "0.0.1"

"""
Given a patched assembly, xtract only chromosomes present 
in the given fasta assembly. Used to filter out patched 
chromosomes corresponding to complete chromosomes in the
NA12878 target assembly.
"""

def parse_fasta(fasta_f):
    """
    Parse a fasta file into a dict, keyed on ID.
    """
    ret = {}
    for fasta in SeqIO.parse(fasta_f, "fasta"):
        ret[fasta.id] = fasta.name
    return ret


def main():
    parser = ArgumentParser(description='Given a patched assembly, xtract only chromosomes present in the given fasta assembly. Used to filter out patched chromosomes corresponding to complete chromosomes in the NA12878 target assembly.')
    parser.add_argument('-p', '--patched_fasta', metavar='FASTA', type=str,
                        required=True, help='Path to patched genome assembly fasta.')
    parser.add_argument('-t', '--target_fasta', metavar='FASTA', type=str,
                        required=True, help='Path to target genome assembly fasta.')
    
    args = parser.parse_args()

    fasta_dict = parse_fasta(args.target_fasta)

    for fasta in SeqIO.parse(args.patched_fasta, "fasta"):
        #sys.stderr.write("%s\n" % (fasta.id))
        if fasta.id in fasta_dict:
            sys.stdout.write("%s\n" % fasta.format("fasta"))
            
    sys.stderr.write("Done!\n")
        
if __name__ == '__main__':
    main()
    exit(0)
