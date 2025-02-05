#!/usr/bin/env python3

import sys, os, re
from argparse import ArgumentParser
import editdistance

__version__ = "0.0.1"

"""
Compare the order of contigs between simulated contigs in the
reference genome and their positions in the final patched genome.
Reports edit distance between simulated and patched ordering vectors
and the percent of adjacent contig pairs that are correct in the
patched versus simulated genomes.
"""

def parse_contigs(contigs_f):
    """
    Parse contigs in a BED format into a dict keyed on chromosome. Each
    key in ret will contain a list of BED entries corresponding to 
    contigs mapped to that chromosome.
    """
    ret = {}
    with open(contigs_f, "r") as contigs:
        for line in contigs:
            contig = line.strip().split()
            ID = contig[0]
            if len(contig) > 8:
                # Scaffold IDs are present (RagTag Results). Since RagTag
                # often splits chromosomes into separat scaffolds, we need
                # to keep track so we don't make invalid coparisons!
                ID = contig[0] + "_" + contig[8]
            if ID in ret:
                ret[ID].append(contig)
            else:
                ret[ID] = [contig]
    return ret


def extract_contig_ordering_vector(contigs_list):
    """
    Extract the list of contig names in the current set.
    """
    ret = []
    for contig in contigs_list:
        ret.append(contig[3])
    return ret


def extract_pairs(contigs_list):
    """
    Extract pairs of adjacent contig names from the given set
    and use as dict keys in ret.
    """
    ret = {}
    i = 0
    while i < len(contigs_list)-1:
        key = (contigs_list[i][3], contigs_list[i+1][3])
        ret[key] = key
        i += 1
    return ret


def count_correct_pairs(reference_pairs, patched_pairs):
    """
    Count the number of correct pairs in patched_pairs
    compared to reference_pairs.
    """
    ret = 0
    for pair in patched_pairs.keys():
        if pair in reference_pairs:
            ret += 1
    return ret


def names_vector_to_rank_vector(names_list):
    """
    Extract ranks from sequence names and return in a list.
    """
    ret = []
    for name in names_list:
        rank = name.split("_")[1]
        ret.append(int(rank))
    return ret


def main():
    parser = ArgumentParser(description='Compare the order of contigs between simulated contigs in the reference genome and their positions in the final patched genome. Reports edit distance between simulated and patched ordering vectors and the percent of adjacent contig pairs that are correct in the patched versus simulated genomes.')
    parser.add_argument('-r', '--reference_contigs', metavar='BED', type=str,
                        required=True, help='BED file containing locations of contigs in the reference sequence.')
    parser.add_argument('-p', '--patched_contigs', metavar='BED', type=str,
                        required=True, help='BED file containing locations of contigs in the patched sequence.')
    
    args = parser.parse_args()

    reference_contigs = parse_contigs(args.reference_contigs)
    patched_contigs = parse_contigs(args.patched_contigs)
    
    # Loop over reference chromosomes, calculating edit distance and pairwise
    # accuracy between reference and patched chromosomes.
    cum_edit_dist = 0
    cum_total_len = 0
    cum_correct_pairs = 0
    cum_patched_pairs = 0
    cum_n_pairs = 0
    cum_n_switches = 0  # Use cum_patched_pairs as denominator
    sys.stdout.write("Chrom\tedit_dist\tnorm_edit_dist\tpct_pairs_recovered\tpct_pairs_correct\tswitch_err_count\tswitch_err_rate\n")
    for scf in patched_contigs.keys():
        #sys.stderr.write("%s\n" % (scf))
        chrom = scf
        if scf not in reference_contigs:
            chrom = scf.split("_")[0]
            if chrom not in reference_contigs:
                sys.stderr.write("%s is not in reference contigs -- skipping!\n" % (chrom))
                continue
        reference_ordering = extract_contig_ordering_vector(reference_contigs[chrom])
        patched_ordering = extract_contig_ordering_vector(patched_contigs[scf])
        #sys.stderr.write("%s\n%s\n" % (len(reference_ordering), len(patched_ordering)))
        edit_dist = editdistance.eval(reference_ordering, patched_ordering)
        norm_edit_dist = edit_dist / len(reference_ordering)
        cum_edit_dist += edit_dist
        cum_total_len += len(reference_ordering)
        reference_pairs = extract_pairs(reference_contigs[chrom])
        patched_pairs = extract_pairs(patched_contigs[scf])
        n_correct_pairs = count_correct_pairs(reference_pairs, patched_pairs)
        cum_correct_pairs += n_correct_pairs
        cum_patched_pairs += len(patched_pairs)
        cum_n_pairs += len(reference_pairs)
        patched_ranks = names_vector_to_rank_vector(patched_ordering)
        sorted_ranks = sorted(patched_ranks)
        # Switch error count is half the edit distance between the observed and
        # sorted ranks. (One switch equals two edit operations.)
        switch_err_count = editdistance.eval(patched_ranks, sorted_ranks) / 2
        cum_n_switches += switch_err_count

        rlen = len(reference_pairs.keys())
        if rlen < 1:
            rlen = 1

        plen = len(patched_pairs.keys())
        if plen < 1:
            plen = 1
                    
        sys.stdout.write("%s\t%d\t%f\t%f\t%f\t%d\t%f\n" % (scf, edit_dist, norm_edit_dist, n_correct_pairs/rlen, n_correct_pairs/plen, switch_err_count, switch_err_count/plen))

    sys.stdout.write("%s\t%d\t%f\t%f\t%f\t%d\t%f\n" % ("Average", cum_edit_dist, cum_edit_dist/cum_total_len, cum_correct_pairs/cum_n_pairs, cum_correct_pairs/cum_patched_pairs, cum_n_switches, cum_n_switches/cum_patched_pairs))


    sys.stderr.write("Done!\n")
        
if __name__ == '__main__':
    main()
    exit(0)
