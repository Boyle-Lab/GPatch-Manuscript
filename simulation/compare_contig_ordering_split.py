#!/usr/bin/env python3

import sys
from argparse import ArgumentParser
import editdistance

import split_benchmark as sb

__version__ = "0.1.0"

"""
Version of compare_contig_ordering.py for benchmark assemblies whose
chromosomes are split into several contigs (see split_benchmark.py).

Comparisons are restricted to one benchmark contig at a time, and each
benchmark contig may lie in either orientation relative to the patched genome:

  reference vector  simulated contigs of the benchmark contig in benchmark order
  patched vector    contigs on the target patched sequence (the one holding most
                    of the benchmark contig's placed sequence) from its first to
                    its last contig, in patched order, reversed when the
                    benchmark contig's orientation is - (majority of placed
                    contigs on the opposite strand to their simulated strand).
                    Contigs from other benchmark contigs placed between them are
                    kept, so interleaving counts as an error; contigs placed on
                    other sequences are missing from the vector.

Edit distance, pair recovery/correctness and switch errors are then computed
as in compare_contig_ordering.py, with one output row per benchmark contig
that has at least one placed contig (benchmark contigs without placed contigs
are skipped, as unplaced chromosomes are in the original). The output columns
are identical to the original; --details writes the orientation and target
sequence of every benchmark contig.
"""


def pairs(names):
    return {(a, b) for a, b in zip(names, names[1:])}


def main():
    parser = ArgumentParser(description='Compare contig order between simulated and patched genomes within each benchmark contig, allowing each benchmark contig to be reversed.')
    parser.add_argument('-r', '--reference_contigs', metavar='BED', type=str, required=True,
                        help='BED file containing locations of contigs in the benchmark contigs.')
    parser.add_argument('-p', '--patched_contigs', metavar='BED', type=str, required=True,
                        help='BED file containing locations of contigs in the patched sequence.')
    parser.add_argument('-d', '--details', metavar='TSV', type=str, default=None,
                        help='Optional per-benchmark-contig details (orientation, target sequence, counts).')
    parser.add_argument('--sep', default=".", help='Separator between chromosome and descriptor in benchmark contig IDs (default ".").')
    args = parser.parse_args()

    sim_by_name, sim_by_bench = sb.parse_sim_bed(args.reference_contigs)
    pat_by_name, pat_by_seq = sb.parse_patched_bed(args.patched_contigs)
    bcs = sb.bench_contigs(sim_by_bench, pat_by_name, args.sep)

    cum_edit_dist = cum_total_len = cum_correct_pairs = cum_patched_pairs = cum_n_pairs = 0
    cum_n_switches = 0
    extra = {}
    sys.stdout.write("Chrom\tedit_dist\tnorm_edit_dist\tpct_pairs_recovered\tpct_pairs_correct\tswitch_err_count\tswitch_err_rate\n")
    for b in sorted(bcs, key=lambda x: (x.split(args.sep)[0], x)):
        c = bcs[b]
        if not c.placed:
            sys.stderr.write("%s has no placed contigs -- skipping!\n" % b)
            continue
        reference_ordering = [r[3] for r in c.sim]
        patched_ordering = [p[3] for p in c.patched_slice(pat_by_seq)]
        edit_dist = editdistance.eval(reference_ordering, patched_ordering)
        cum_edit_dist += edit_dist
        cum_total_len += len(reference_ordering)
        reference_pairs = pairs(reference_ordering)
        patched_pairs = pairs(patched_ordering)
        n_correct_pairs = len(patched_pairs & reference_pairs)
        cum_correct_pairs += n_correct_pairs
        cum_patched_pairs += len(patched_pairs)
        cum_n_pairs += len(reference_pairs)
        own = {r[3] for r in c.sim}
        patched_ranks = [sb.rank(n) for n in patched_ordering if n in own]
        switch_err_count = editdistance.eval(patched_ranks, sorted(patched_ranks)) / 2
        cum_n_switches += switch_err_count
        rlen = max(len(reference_pairs), 1)
        plen = max(len(patched_pairs), 1)
        sys.stdout.write("%s\t%d\t%f\t%f\t%f\t%d\t%f\n" % (
            b, edit_dist, edit_dist / len(reference_ordering), n_correct_pairs / rlen,
            n_correct_pairs / plen, switch_err_count, switch_err_count / plen))
        extra[b] = {"foreign_in_slice": len(patched_ordering) - len(patched_ranks)}

    sys.stdout.write("%s\t%d\t%f\t%f\t%f\t%d\t%f\n" % (
        "Average", cum_edit_dist, cum_edit_dist / cum_total_len, cum_correct_pairs / max(cum_n_pairs, 1),
        cum_correct_pairs / max(cum_patched_pairs, 1), cum_n_switches, cum_n_switches / max(cum_patched_pairs, 1)))
    if args.details:
        sb.write_bench_summary(args.details, bcs, extra)
    sys.stderr.write("Done!\n")


if __name__ == '__main__':
    main()
    exit(0)
