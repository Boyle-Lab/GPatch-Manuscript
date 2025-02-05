#!/usr/bin/env python3

import sys, os, re
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Seq import Seq
#import editdistance
import textdistance

__version__ = "0.0.1"

"""
Compare common sequences between to fasta files by edit distance.
"""

def parse_fasta(fasta_f):
    """
    Parse a fasta sequence into a dict keyed on id.
    """
    ret = {}
    for fasta in SeqIO.parse(fasta_f, "fasta"):
        ret[fasta.id] = fasta
    return ret


def main():
    parser = ArgumentParser(description='Compare common sequences between to fasta files by edit distance.')
    parser.add_argument('-r', '--reference_fasta', metavar='FASTA', type=str,
                        required=True, help='Path to genome assembly fasta containing reference sequences used to produce the patched genome.')
    parser.add_argument('-p', '--patched_fasta', metavar='FASTA', type=str,
                        required=True, help='Path to patched fasta sequence assembly.')
    
    
    args = parser.parse_args()

    reference_fasta = parse_fasta(args.reference_fasta)

    cum_rdist = 0
    cum_sdist = 0
    cum_rlen = 0
    cum_reflen = 0
    cum_pchlen = 0
    sys.stdout.write("chrom\tref_len\tpatched_len\tpct_coverage\traw_bag_dist\tnorm_bag_dist\tscaled_bag_dist\tnorm_scaled_bag_dist\n")
    for patched_fasta in SeqIO.parse(args.patched_fasta, "fasta"):
        if patched_fasta.id in reference_fasta:
            ID = patched_fasta.id
            if re.search("scf", patched_fasta.description):                
                ID = patched_fasta.id + "_" + patched_fasta.description.split()[1]
            ref_seq = str(reference_fasta[patched_fasta.id].seq).lower()
            patched_seq = str(patched_fasta.seq).lower()
            rlen = max(len(ref_seq), len(patched_seq))
            # Raw bag distance
            raw_dist = textdistance.bag(ref_seq, patched_seq)
            # Bag distance minus difference in length between patched and ref sequence
            scaled_dist = raw_dist - abs(len(ref_seq)-len(patched_seq))
            norm_rdist = raw_dist / rlen
            norm_sdist = scaled_dist / rlen
            sys.stdout.write("%s\t%d\t%d\t%f\t%d\t%f\t%d\t%f\n" % (ID, len(ref_seq), len(patched_seq), len(patched_seq)/len(ref_seq), raw_dist, norm_rdist, scaled_dist, norm_sdist))
            cum_rdist += raw_dist
            cum_sdist += scaled_dist
            cum_rlen += rlen
            cum_reflen += len(ref_seq)
            cum_pchlen += len(patched_seq)
    sys.stdout.write("%s\t%d\t%d\t%f\t%d\t%f\t%d\t%f\n" % ("average", cum_reflen, cum_pchlen, cum_pchlen/cum_reflen, cum_rdist, cum_rdist/cum_rlen, cum_sdist, cum_sdist/cum_rlen))
    sys.stderr.write("Done!\n")
        
if __name__ == '__main__':
    main()
    exit(0)
