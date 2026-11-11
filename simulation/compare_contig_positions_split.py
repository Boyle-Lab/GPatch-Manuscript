#!/usr/bin/env python3

import sys
from argparse import ArgumentParser
from collections import Counter
import pysam

import split_benchmark as sb

__version__ = "0.1.0"

"""
Version of compare_contig_positions.py for benchmark assemblies whose
chromosomes are split into several contigs (see split_benchmark.py).

Benchmark-contig coordinates cannot be compared directly with patched-genome
coordinates, and each benchmark contig may lie in either orientation, so each
simulated contig is evaluated relative to the other contigs from its own
benchmark contig. With orientation o of the benchmark contig and simulated
interval [s, e):

  is_mapped              the contig is placed in the patched genome
  is_mapped_correctly    placed on the target patched sequence (the one holding
                         most of the benchmark contig's placed sequence), whose
                         name matches the benchmark chromosome, at the
                         benchmark contig's common offset: patched start - s
                         and patched end - e (o = +), or patched start + e and
                         patched end + s (o = -), both equal the most frequent
                         such offset among the benchmark contig's contigs
                         (ties: the offset covering most sequence). This is
                         the relative analogue of the original exact-coordinate
                         test; a benchmark contig with one placed contig is
                         trivially at its own offset.
  is_correct_orientation patched strand = simulated strand (o = +) or the
                         opposite strand (o = -)
  is_correct_order       taking the benchmark contig's placed contigs in
                         benchmark order (reversed for o = -): the first is in
                         correct order if no other contig of its benchmark
                         contig precedes it on the target sequence; each later
                         one if it lies on the same patched sequence as the
                         previous placed one and after it (by patched start).

Output columns are identical to the original (ref_chrom is the benchmark
contig ID); --details writes per-benchmark-contig orientation and offset.
"""


def parse_bam(bam_f):
    ret = {}
    for read in pysam.AlignmentFile(bam_f, "rb"):
        if not read.is_supplementary and not read.is_secondary:
            ret[read.query_name] = read
    return ret


def offsets(c, sim_rec, p):
    s, e = int(sim_rec[1]), int(sim_rec[2])
    ps, pe = int(p[1]), int(p[2])
    return (ps - s, pe - e) if c.orientation == "+" else (ps + e, pe + s)


def main():
    parser = ArgumentParser(description='Check placement, orientation and order of simulated contigs in the patched genome within each benchmark contig, allowing each benchmark contig to be reversed.')
    parser.add_argument('-r', '--reference_contigs', metavar='BED', type=str, required=True,
                        help='BED file containing locations of contigs in the benchmark contigs.')
    parser.add_argument('-p', '--patched_contigs', metavar='BED', type=str, required=True,
                        help='BED file containing locations of contigs in the patched sequence.')
    parser.add_argument('-b', '--bam', metavar='BAM', type=str, default="",
                        help='BAM file originally used to generate the patched genome.')
    parser.add_argument('-d', '--details', metavar='TSV', type=str, default=None,
                        help='Optional per-benchmark-contig details (orientation, target sequence, offset).')
    parser.add_argument('--sep', default=".", help='Separator between chromosome and descriptor in benchmark contig IDs (default ".").')
    args = parser.parse_args()

    sim_by_name, sim_by_bench = sb.parse_sim_bed(args.reference_contigs)
    pat_by_name, pat_by_seq = sb.parse_patched_bed(args.patched_contigs)
    bcs = sb.bench_contigs(sim_by_bench, pat_by_name, args.sep)
    contigs_bam = parse_bam(args.bam) if args.bam else None

    # Per benchmark contig: common offset, and per-contig order calls.
    offset, order_ok, extra = {}, {}, {}
    for b, c in bcs.items():
        votes, bp = Counter(), Counter()
        for p in c.on_target:
            ds, de = offsets(c, sim_by_name[p[3]], p)
            if ds == de:
                votes[ds] += 1
                bp[ds] += int(p[2]) - int(p[1])
        offset[b] = max(votes, key=lambda k: (votes[k], bp[k])) if votes else None
        placed_norm = [r for r in c.normalized_sim() if r[3] in pat_by_name]
        # In normalized (orientation-corrected) benchmark order, patched starts
        # increase when the order is correct.
        prev = None
        for r in placed_norm:
            p = pat_by_name[r[3]]
            if prev is None:
                ok = p[-1] == c.target and all(int(q[1]) >= int(p[1]) for q in c.on_target)
            else:
                ok = p[-1] == prev[-1] and int(p[1]) > int(prev[1])
            order_ok[r[3]] = ok
            prev = p
        extra[b] = {"offset": offset[b] if offset[b] is not None else "NA",
                    "n_at_offset": votes[offset[b]] if offset[b] is not None else 0}

    sys.stdout.write("ref_chrom\tref_start\tref_end\tref_contig\tref_length\tmapped_chrom\tmapped_start\tmapped_end\tmapped_length\tis_mapped\tis_mapped_correctly\tis_correct_orientation\tis_correct_order\tmapping_quality\n")
    for b in sorted(bcs, key=lambda x: (x.split(args.sep)[0], x)):
        c = bcs[b]
        for r in c.sim:
            mapq = "."
            if contigs_bam is not None and r[3] in contigs_bam:
                mapq = contigs_bam[r[3]].mapping_quality
            if r[3] not in pat_by_name:
                sys.stdout.write("%s\t%d\t%d\t%s\t%d\t.\t.\t.\t.\tFalse\tFalse\tFalse\tFalse\t%s\n" % (
                    r[0], int(r[1]), int(r[2]), r[3], int(r[4]), mapq))
                continue
            p = pat_by_name[r[3]]
            ds, de = offsets(c, r, p)
            is_mapped_correctly = (p[-1] == c.target and p[0] == c.chrom
                                   and ds == de == offset[b])
            is_correct_orientation = p[5] == c.expected_strand(r)
            sys.stdout.write("%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\n" % (
                r[0], int(r[1]), int(r[2]), r[3], int(r[4]), p[0], int(p[1]), int(p[2]), int(p[7]),
                True, is_mapped_correctly, is_correct_orientation, order_ok[r[3]], mapq))

    if args.details:
        sb.write_bench_summary(args.details, bcs, extra)
    sys.stderr.write("Done!\n")


if __name__ == '__main__':
    main()
    exit(0)
