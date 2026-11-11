#!/usr/bin/env python3

import sys
from argparse import ArgumentParser
from collections import defaultdict
import pysam
import textdistance

import split_benchmark as sb

__version__ = "0.1.0"

"""
Version of compare_fasta_by_editdist.py for benchmark assemblies whose
chromosomes are split into several contigs (see split_benchmark.py).

Each benchmark contig is compared with the part of the patched genome that
corresponds to it, oriented to match it (reverse-complemented when the
benchmark contig's orientation is -):

  chrom scope  the benchmark contig is the only one for its chromosome (a
               complete chromosome) and its target patched sequence is named
               after that chromosome: the whole patched sequence is compared
               with the whole benchmark contig, as in the original script.
  span scope   otherwise (and for all benchmark contigs with --scope span):
               the patched sequence from the start of the first to the end of
               the last of the benchmark contig's contigs on the target
               sequence, compared with the benchmark contig from the start to
               the end of the same simulated contigs. Contigs from other
               benchmark contigs placed in between are included.

Benchmark contigs with no placed contigs are skipped (as unplaced chromosomes
are in the original). Distances are bag distances, as in the original, which
are insensitive to order; output columns are identical to the original
(chrom is the benchmark contig ID). --details writes the scope, orientation
and compared intervals.
"""

COMP = str.maketrans("acgtnACGTN", "tgcanTGCAN")


def main():
    parser = ArgumentParser(description='Compare each benchmark contig with the corresponding part of the patched genome by bag distance, allowing each benchmark contig to be reversed.')
    parser.add_argument('-r', '--reference_fasta', metavar='FASTA', type=str, required=True,
                        help='Benchmark assembly FASTA with unique benchmark contig IDs (indexed or indexable).')
    parser.add_argument('-p', '--patched_fasta', metavar='FASTA', type=str, required=True,
                        help='Patched FASTA sequence assembly (indexed or indexable).')
    parser.add_argument('-s', '--sim_contigs', metavar='BED', type=str, required=True,
                        help='BED file of simulated contig locations in the benchmark contigs.')
    parser.add_argument('-g', '--patched_contigs', metavar='BED', type=str, required=True,
                        help='BED file of contig locations in the patched sequence.')
    parser.add_argument('--scope', choices=("auto", "span"), default="auto",
                        help='auto: whole-sequence comparison for complete chromosomes, span otherwise (default); span: span comparison for all.')
    parser.add_argument('--fasta-id-col', type=int, default=1,
                        help='Column of the patched BED holding the patched FASTA sequence ID (default 1; e.g. 9 for RagTag-derived BEDs).')
    parser.add_argument('-d', '--details', metavar='TSV', type=str, default=None,
                        help='Optional per-benchmark-contig details (scope, orientation, compared intervals).')
    parser.add_argument('--sep', default=".", help='Separator between chromosome and descriptor in benchmark contig IDs (default ".").')
    args = parser.parse_args()

    sim_by_name, sim_by_bench = sb.parse_sim_bed(args.sim_contigs)
    pat_by_name, pat_by_seq = sb.parse_patched_bed(args.patched_contigs)
    bcs = sb.bench_contigs(sim_by_bench, pat_by_name, args.sep)
    ref = pysam.FastaFile(args.reference_fasta)
    pat = pysam.FastaFile(args.patched_fasta)
    n_per_chrom = defaultdict(int)
    for b in ref.references:
        n_per_chrom[sb.bench_chrom(b, args.sep)] += 1

    cum_rdist = cum_sdist = cum_rlen = cum_reflen = cum_pchlen = 0
    extra = {}
    sys.stdout.write("chrom\tref_len\tpatched_len\tpct_coverage\traw_bag_dist\tnorm_bag_dist\tscaled_bag_dist\tnorm_scaled_bag_dist\n")
    for b in sorted(bcs, key=lambda x: (x.split(args.sep)[0], x)):
        c = bcs[b]
        if not c.on_target:
            sys.stderr.write("%s has no placed contigs -- skipping!\n" % b)
            continue
        fasta_id = c.on_target[0][args.fasta_id_col - 1]
        complete = n_per_chrom[c.chrom] == 1 and c.target.split("_")[0] == c.chrom
        if args.scope == "auto" and complete:
            scope = "chrom"
            rs, re_ = 0, ref.get_reference_length(b)
            ps, pe = 0, pat.get_reference_length(fasta_id)
        else:
            scope = "span"
            sims = [sim_by_name[p[3]] for p in c.on_target]
            rs, re_ = min(int(r[1]) for r in sims), max(int(r[2]) for r in sims)
            ps, pe = min(int(p[1]) for p in c.on_target), max(int(p[2]) for p in c.on_target)
        ref_seq = ref.fetch(b, rs, re_).lower()
        patched_seq = pat.fetch(fasta_id, ps, pe).lower()
        if c.orientation == "-":
            patched_seq = patched_seq.translate(COMP)[::-1]
        rlen = max(len(ref_seq), len(patched_seq))
        raw_dist = textdistance.bag(ref_seq, patched_seq)
        scaled_dist = raw_dist - abs(len(ref_seq) - len(patched_seq))
        sys.stdout.write("%s\t%d\t%d\t%f\t%d\t%f\t%d\t%f\n" % (
            b, len(ref_seq), len(patched_seq), len(patched_seq) / len(ref_seq), raw_dist,
            raw_dist / rlen, scaled_dist, scaled_dist / rlen))
        cum_rdist += raw_dist
        cum_sdist += scaled_dist
        cum_rlen += rlen
        cum_reflen += len(ref_seq)
        cum_pchlen += len(patched_seq)
        extra[b] = {"scope": scope, "bench_interval": "%d-%d" % (rs, re_),
                    "patched_interval": "%s:%d-%d" % (fasta_id, ps, pe),
                    "bench_contig_len": ref.get_reference_length(b)}
    sys.stdout.write("%s\t%d\t%d\t%f\t%d\t%f\t%d\t%f\n" % (
        "average", cum_reflen, cum_pchlen, cum_pchlen / cum_reflen, cum_rdist, cum_rdist / cum_rlen,
        cum_sdist, cum_sdist / cum_rlen))
    if args.details:
        sb.write_bench_summary(args.details, bcs, extra)
    sys.stderr.write("Done!\n")


if __name__ == '__main__':
    main()
    exit(0)
