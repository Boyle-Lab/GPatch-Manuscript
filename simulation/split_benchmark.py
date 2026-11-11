#!/usr/bin/env python3
"""
Shared helpers for comparing a patched genome against a benchmark assembly
whose chromosomes may be split into several contigs (e.g. NA12878.1, where
10 of 23 chromosomes are fragmented).

The simulated contig BED is in benchmark-contig coordinates: column 1 is the
benchmark contig ID (made unique as <chrom>.<descriptor>, e.g.
chr1.mat-0000003, by T2T_NA12878/rename_benchmark_contigs.sh) and simulated
contigs are named <benchmark contig>_<rank>. Because the order and orientation
of benchmark contigs from the same chromosome are unknown, all comparisons are
made within one benchmark contig at a time, and each benchmark contig may lie
in either orientation relative to the patched genome (whose coordinates follow
the patching reference). The orientation of a benchmark contig is the
length-weighted majority vote of its placed simulated contigs: + when the
patched strand equals the simulated strand, - when it is the opposite.
"""

from collections import defaultdict, Counter

__version__ = "0.1.0"

FLIP = {"+": "-", "-": "+"}


def bench_chrom(bench_id, sep="."):
    """Chromosome label of a benchmark contig ID (part before the first sep)."""
    return bench_id.split(sep)[0]


def rank(name):
    """Rank of a simulated contig named <benchmark contig>_<rank>."""
    return int(name.rsplit("_", 1)[1])


def parse_sim_bed(bed_f):
    """
    Simulated contig BED (benchmark coordinates).
    Returns (by_name: name -> record, by_bench: bench_id -> [records sorted by start]).
    Records are the split BED fields.
    """
    by_name, by_bench = {}, defaultdict(list)
    with open(bed_f) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            r = line.rstrip("\n").split("\t")
            by_name[r[3]] = r
            by_bench[r[0]].append(r)
    for v in by_bench.values():
        v.sort(key=lambda r: (int(r[1]), -int(r[2])))
    return by_name, by_bench


def parse_patched_bed(bed_f):
    """
    Patched contig BED (GPatch contigs.bed, or RagTag-derived BED with a 9th
    scaffold column). The patched sequence key is column 1, or column 1 +
    "_" + column 9 when more than 8 columns are present (as in
    compare_contig_ordering.py), since RagTag can split a chromosome into
    several scaffolds.
    Returns (by_name: name -> record, by_seq: key -> [records sorted by start]).
    Each record gets an extra trailing element: the sequence key.
    """
    by_name, by_seq = {}, defaultdict(list)
    with open(bed_f) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            r = line.rstrip("\n").split("\t")
            key = r[0] + "_" + r[8] if len(r) > 8 else r[0]
            r.append(key)
            by_name[r[3]] = r
            by_seq[key].append(r)
    for v in by_seq.values():
        v.sort(key=lambda r: (int(r[1]), int(r[2])))
    return by_name, by_seq


class BenchContig:
    """Placement summary of the simulated contigs from one benchmark contig."""

    def __init__(self, bench_id, sim_recs, patched_by_name, sep="."):
        self.id = bench_id
        self.chrom = bench_chrom(bench_id, sep)
        self.sim = sim_recs
        self.placed = [patched_by_name[r[3]] for r in sim_recs if r[3] in patched_by_name]
        sim = {r[3]: r for r in sim_recs}
        # Target patched sequence: the one holding most of the placed contigs'
        # (simulated) length; ties prefer the sequence named after the
        # benchmark chromosome.
        bp = Counter()
        for p in self.placed:
            bp[p[-1]] += int(sim[p[3]][2]) - int(sim[p[3]][1])
        self.target = (max(bp, key=lambda k: (bp[k], k.split("_")[0] == self.chrom)) if bp else None)
        # Orientation: length-weighted strand agreement over all placed contigs.
        agree = dis = 0
        for p in self.placed:
            L = int(sim[p[3]][2]) - int(sim[p[3]][1])
            if p[5] == sim[p[3]][5]:
                agree += L
            else:
                dis += L
        self.orientation = "-" if dis > agree else "+"
        tot = agree + dis
        self.strand_agreement = (max(agree, dis) / tot) if tot else float("nan")
        # Contigs on the target sequence, in patched order.
        self.on_target = [p for p in self.placed if p[-1] == self.target]
        self.on_target.sort(key=lambda p: (int(p[1]), int(p[2])))

    def normalized_sim(self):
        """Simulated contigs in benchmark order, reversed when orientation is -."""
        return self.sim if self.orientation == "+" else self.sim[::-1]

    def expected_strand(self, sim_rec):
        return sim_rec[5] if self.orientation == "+" else FLIP[sim_rec[5]]

    def patched_slice(self, patched_by_seq):
        """
        Contigs on the target sequence from the first to the last contig of this
        benchmark contig (patched order), including contigs from other benchmark
        contigs placed between them; reversed when orientation is -.
        """
        if not self.target:
            return []
        seq = patched_by_seq[self.target]
        idx = [i for i, p in enumerate(seq) if p[3] in {r[3] for r in self.sim}]
        sl = seq[idx[0]:idx[-1] + 1]
        return sl if self.orientation == "+" else sl[::-1]


def bench_contigs(sim_by_bench, patched_by_name, sep="."):
    return {b: BenchContig(b, recs, patched_by_name, sep) for b, recs in sim_by_bench.items()}


def write_bench_summary(out_f, bcs, extra=None):
    """Per-benchmark-contig placement summary (orientation, target, counts)."""
    with open(out_f, "w") as o:
        cols = ["bench_contig", "chrom", "n_sim", "n_placed", "n_on_target", "target_seq",
                "target_matches_chrom", "orientation", "strand_agreement"]
        ek = sorted({k for d in (extra or {}).values() for k in d})
        o.write("\t".join(cols + ek) + "\n")
        for b in sorted(bcs):
            c = bcs[b]
            row = [b, c.chrom, len(c.sim), len(c.placed), len(c.on_target), c.target or ".",
                   (c.target or "").split("_")[0] == c.chrom if c.target else ".", c.orientation,
                   "%.4f" % c.strand_agreement if c.placed else "NA"]
            row += [(extra or {}).get(b, {}).get(k, "NA") for k in ek]
            o.write("\t".join(map(str, row)) + "\n")
