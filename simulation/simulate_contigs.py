#!/usr/bin/env python3

import sys, os, re
from random import randrange
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Seq import Seq

__version__ = "0.0.1"

"""
Simulate a contig-based genome from a chromosome-scale genome by 
breaking chromosomes into pseudo-contigs. Splitting is based on
a set of contig lengths extracted from a contig-level "model" genome
assembly.
"""

def read_lengths(lengths_f):
    """
    Read the lengths vector used to specify pseudo-contig lengths.
    """
    ret = []
    with open(lengths_f, "r") as lengths:
        for rec in lengths:
            rec.strip()
            ret.append(int(rec))
    return ret


def main():
    parser = ArgumentParser(description='Simulate a contig-based genome from a chromosome-scale genome by breaking chromosomes into pseudo-contigs. Splitting is based on a set of contig lengths extracted from a contig-level "model" genome assembly.')
    parser.add_argument('-f', '--fasta', metavar='FASTA', type=str,
                        required=True, help='Path to genome assembly fasta containing chromosomes to break into pseudocontigs.')
    parser.add_argument('-l', '--lengths', metavar='STR', type=str,
                        required=True, help='File containing contig lengths extracted from a "model" genome assembly. These will be used as the lengths of pseudo-contigs.')
    parser.add_argument('-p', '--prefix', metavar='STR', type=str, default="",
                        required=False, help='Prefix for output fasta and bed file names.')
    
    args = parser.parse_args()

    contig_lengths = read_lengths(args.lengths)

    contigs_bed_f = "contigs.bed"
    if args.prefix != "":
        contigs_bed_f = args.prefix + "." + contigs_bed_f

    contigs_fasta_f = "contigs.fa"
    if args.prefix != "":
        contigs_fasta_f = args.prefix + "." + contigs_fasta_f

    contigs_bed = open(contigs_bed_f, "w")
    contigs_fasta = open(contigs_fasta_f, "w")
    
    for fasta in SeqIO.parse(args.fasta, "fasta"):
        root_id = fasta.id
        seq = fasta.seq
        pos = 0    # Start position in fasta sequence
        frag = 0   # Fragment number

        # From here, we need to break the chromosome into some number of
        # fragments, with random lengths selected based on a set of observed
        # contig lengths from a model assembly.
        while pos < len(seq):
            # Select a random index from the sequence lengths. Note we
            # are not tracking which indeces are used already, so lengths
            # can be reused. This way we ensure we don't run out of
            # indeces, which would be problematic.
            idx = randrange(len(contig_lengths))
            contig_len = contig_lengths[idx]
            cend = pos + contig_len
            if cend > len(seq):
                # Truncate at sequence end if needed.
                cend = len(seq)
            cseq = seq[pos:cend]
            cname = root_id + '_' + str(frag)
            fasta.seq = seq[pos:cend]
            fasta.id = cname
            contigs_fasta.write("%s\n" % fasta.format("fasta"))
            contigs_bed.write("%s\t%d\t%d\t%s\t%d\n" % (root_id, pos, cend, cname, cend-pos))
            pos = cend
            frag += 1

    contigs_bed.close()
    contigs_fasta.close()
    sys.stderr.write("Done!\n")
        
if __name__ == '__main__':
    main()
    exit(0)
