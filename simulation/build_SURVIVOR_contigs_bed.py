#!/usr/bin/env python3

import sys, os, re
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

__version__ = "0.0.1"

"""
Given a set of contigs containing indels introduced by
SURVIVOR simSV, build a contigs.bed file with correct
genomic breakpoints for the pseudo-assembly corresponding
to the concatenation of all contigs according to their
known order in the reference genome source.
"""

def parse_fasta(fasta_f):
    """
    Parse a fasta sequence into a dict keyed on fasta name.
    """
    ret = {}
    for fasta in SeqIO.parse(fasta_f, "fasta"):
        ret[fasta.name] = fasta
    return ret


def parse_bed(bed_f):
    """
    Parse the given BED file into a dict keyed on chrom, where
    each key contains a list of BED records.
    """
    ret = {}
    with open(bed_f, "r") as bed:
        for line in bed:
            rec = line.strip().split()
            if rec[0] in ret:
                ret[rec[0]].append(rec)
            else:
                ret[rec[0]] = [rec]
    return ret


def main():
    parser = ArgumentParser(description='Given a set of contigs containing indels introduced by SURVIVOR simSV, build a contigs.bed file with correct genomic breakpoints for the pseudo-assembly corresponding to the concatenation of all contigs according to their known order in the reference genome source.')
    parser.add_argument('-s', '--survivor_fasta', metavar='FASTA', type=str,
                        required=True, help='Path to fasta results from SURVIVOR.')
    parser.add_argument('-r', '--reference_bed', metavar='BED', type=str,
                        required=True, help='BED file containing entries for all reference elements supplied to SURVIVOR to produce the mutated fasta. This determines contig order in the simulated fasta seq and must be sorted by position or the contig order will be wrong!')
    parser.add_argument('-f', '--fasta_out', metavar='FASTA', type=str,
                        required=False, default="", help='Store concatenated fasta sequence to given file.')
    
    
    args = parser.parse_args()

    survivor_fasta = parse_fasta(args.survivor_fasta)
    reference_bed = parse_bed(args.reference_bed)

    # Loop over chromosomes, using the inherent order of the BED records
    # to produce a contig ordering. We will build the actual sequence and
    # keep track of position in the sequence to get the contig coordinates
    # in the mutated genome frame.
    if args.fasta_out != "":
        # We will write concatenated fasta sequences to a file.
        fasta_out = open(args.fasta_out, "w")
    for chrom in reference_bed.keys():
        cseq = ""
        pos = 0
        for contig in reference_bed[chrom]:
            seq = survivor_fasta[contig[3]].seq
            if contig[5] == "-":
                seq = seq.reverse_complement()
            cseq = cseq + str(seq)
            sys.stdout.write("%s\t%d\t%d\t%s\t%d\t%s\n" % (contig[0], pos, len(cseq), contig[3], len(cseq)-pos, contig[5]))
            pos = len(cseq)
        if args.fasta_out != "":
            # Write the fasta sequence to output.            
            rec = SeqRecord(Seq(cseq), id = chrom)
            fasta_out.write("%s\n" % (rec.format("fasta")))

    if args.fasta_out != "":
        fasta_out.close()
        
    sys.stderr.write("Done!\n")
        
if __name__ == '__main__':
    main()
    exit(0)
