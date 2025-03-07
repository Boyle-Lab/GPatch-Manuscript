#!/usr/bin/env python3

import sys, os, re
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Seq import Seq

__version__ = "0.0.1"

"""
Extract only complete chromosomes (those with only a single contig) from
the given fasta assembly.
"""

def parse_fasta(fasta_f):
    """
    Parse a fasta file into a dict, keyed on ID. Duplicated
    IDs will be stored in lists.
    """
    ret = {}
    for fasta in SeqIO.parse(fasta_f, "fasta"):
        if fasta.id in ret:
            ret[fasta.id].append(fasta.name)
        else:
            ret[fasta.id] = [fasta.name]
    return ret


def main():
    parser = ArgumentParser(description='Extract only complete chromosomes (those with only a single contig) from the given fasta assembly. ')
    parser.add_argument('-f', '--fasta', metavar='FASTA', type=str,
                        required=True, help='Path to genome assembly fasta with complete and partial chromosomes. Partial chromsomes are assumed to share the same ID.')
    
    args = parser.parse_args()

    fasta_dict = parse_fasta(args.fasta)

    for fasta in SeqIO.parse(args.fasta, "fasta"):
        #sys.stderr.write("%s\n" % (fasta.id))
        if fasta.name in fasta_dict:
            if len(fasta_dict[fasta.id]) == 1:
                # Single contig. Print fasta to output.
                sys.stdout.write("%s\n" % fasta.format("fasta"))
        else:
            sys.stderr.write("WARNING: %s not found in FASTA ID dict! This should not happen!\n" % (fasta.id))
            
    sys.stderr.write("Done!\n")
        
if __name__ == '__main__':
    main()
    exit(0)
