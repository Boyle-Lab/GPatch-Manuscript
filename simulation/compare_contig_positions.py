#!/usr/bin/env python3

import sys, os, re
from argparse import ArgumentParser
import pysam

__version__ = "0.0.1"

"""
Compare the locations of contigs between simulated contigs in the
reference genome and their positions in the final patched genome.
Reports whether each contig is (correctly) mapped according to its
chromosomal coordinates, order, and orientation in the patched genome.
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


def parse_contigs(contigs_f):
    """
    Parse contigs in a BED format into a dict keyed on BED Name.
    """
    ret = {}
    with open(contigs_f, "r") as contigs:
        for line in contigs:
            contig = line.strip().split()
            ret[contig[3]] = contig
    return ret


def main():
    parser = ArgumentParser(description='Simulate a contig-based genome from a chromosome-scale genome by breaking chromosomes into pseudo-contigs. Splitting is based on a set of contig lengths extracted from a contig-level "model" genome assembly.')
    parser.add_argument('-r', '--reference_contigs', metavar='FASTA', type=str,
                        required=True, help='BED file containing locations of contigs in the reference sequence.')
    parser.add_argument('-p', '--patched_contigs', metavar='STR', type=str,
                        required=True, help='BED file containing locations of contigs in the patched sequence.')
    parser.add_argument('-b', '--bam', metavar='STR', type=str, default="",
                        required=False, help='BAM file originally used to generate the patched genome.')
    
    args = parser.parse_args()

    patched_contigs = parse_contigs(args.patched_contigs)
    contigs_bam = None
    if args.bam != "":
        contigs_bam = parse_bam(args.bam)

    # Loop through reference contig records, checking each one
    # for correct placement in the patched genome.
    sys.stdout.write("ref_chrom\tref_start\tref_end\tref_contig\tref_length\tmapped_chrom\tmapped_start\tmapped_end\tmapped_length\tis_mapped\tis_mapped_correctly\tis_correct_orientation\tis_correct_order\tmapping_quality\n")
    previous_reference_contig = []
    with open(args.reference_contigs) as reference_contigs:
        for line in reference_contigs:
            reference_contig = line.strip().split()
            mapq = "."
            if contigs_bam is not None:
                aligned_contig = contigs_bam[reference_contig[3]]
                mapq = aligned_contig.mapping_quality

            # First see if there is a mapped contig in patched_contigs
            is_mapped = False
            is_mapped_correctly = False
            is_correct_orientation = False
            is_correct_order = False
            if reference_contig[3] not in patched_contigs:
                sys.stdout.write("%s\t%d\t%d\t%s\t%d\t.\t.\t.\t.\tFalse\tFalse\tFalse\tFalse\t%s\n" % (reference_contig[0], int(reference_contig[1]), int(reference_contig[2]), reference_contig[3], int(reference_contig[4]), mapq))
            else:
                patched_contig = patched_contigs[reference_contig[3]]
                is_mapped = True
                # Check chromosome and position (note we require exact match)
                if patched_contig[0] == reference_contig[0]:
                    if int(patched_contig[1]) == int(reference_contig[1]):
                        if int(patched_contig[2]) == int(reference_contig[2]):
                            is_mapped_correctly = True
                else:
                    is_mapped_correctly = False
                # Check for proper orientation. All should map to the + strand.
                if patched_contig[5] == reference_contig[5]:
                    is_correct_orientation = True
                # Check for correct order. We will do this by comparing to the
                # mapped position of the previous reference contig.
                if re.search("_0", reference_contig[3]) and patched_contig[3] == reference_contig[3]:
                    # Starting a new chromosome with the correct contig
                    is_correct_order = True
                elif patched_contig[0] == previous_reference_contig[0]:
                    if int(patched_contig[1]) >= int(previous_reference_contig[2]):
                        is_correct_order = True
                sys.stdout.write("%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\n" % (reference_contig[0], int(reference_contig[1]), int(reference_contig[2]), reference_contig[3], int(reference_contig[4]), patched_contig[0], int(patched_contig[1]), int(patched_contig[2]), int(patched_contig[7]), is_mapped, is_mapped_correctly, is_correct_orientation, is_correct_order, mapq))
                        
            previous_reference_contig = reference_contig

    sys.stderr.write("Done!\n")
        
if __name__ == '__main__':
    main()
    exit(0)
