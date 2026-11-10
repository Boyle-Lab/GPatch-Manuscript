#!/usr/bin/env python3
"""
Compare GPatch junction handling (trims and reference patches) against the
true junctions of the simulated contig assembly.

GPatch contigs.bed (patched-genome frame) columns:
  chrom start end contig . strand qstart qend
qstart > 0: qstart bp were trimmed from the 5' end of the placed contig. A gap
between the previous contig's end and this contig's start (patched frame) is a
reference patch of that length; otherwise the contigs are bookended.

Truth comes from the simulation BED (benchmark coordinates) and junctions.tsv
(distance < 0 = overlap of -distance bp, > 0 = gap).

Each consecutive pair of placed contigs is classified:
  true_adjacent  a true junction in junctions.tsv, same left/right order
  true_reversed  a true junction, placed in swapped order
  not_adjacent   same benchmark chromosome, not adjacent (intervening contigs
                 not placed); true distance from the BED intervals
  diff_chrom     contigs from different benchmark chromosomes (no truth)
Chromosome ends (START/END rows) compare the leading/terminal patch (or trim)
with the benchmark sequence before the first / after the last placed contig,
when the contig comes from the same-named benchmark chromosome.

For every unit, gpatch_dist = patch length (>0), -trim (<0) or 0 (bookend);
true_dist = true gap (>0) or -overlap (<0); diff = gpatch_dist - true_dist:
> 0 means extra sequence (duplicated overlap or excess patch), < 0 means
missing sequence (over-trim or short/absent patch).
"""

import sys
from argparse import ArgumentParser
from collections import Counter

__version__ = "0.2.0"


def stats(d):
    d = sorted(d)
    n = len(d)
    if not n:
        return "n=0"
    ad = sorted(abs(x) for x in d)
    return ("n=%d exact=%d |d|<=10=%d |d|<=1kb=%d |d|<=10kb=%d median|d|=%d max|d|=%d extra=%d missing=%d"
            % (n, sum(x == 0 for x in d), sum(x <= 10 for x in ad), sum(x <= 1000 for x in ad),
               sum(x <= 10000 for x in ad), ad[n // 2], ad[-1],
               sum(x for x in d if x > 0), -sum(x for x in d if x < 0)))


def main():
    p = ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("-g", "--gpatch-bed", required=True)
    p.add_argument("-P", "--patches-bed", required=True, help="patches.bed (reference frame).")
    p.add_argument("-s", "--sim-bed", required=True)
    p.add_argument("-j", "--junctions", required=True)
    p.add_argument("--ref-fai", required=True, help="Reference .fai (chromosome lengths).")
    p.add_argument("--bench-paf", required=True, help="Simulation PAF (benchmark chromosome lengths).")
    p.add_argument("-l", "--label", default="")
    p.add_argument("-o", "--out", required=True, help="Per-junction TSV.")
    a = p.parse_args()

    rlen = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(a.ref_fai)}
    blen = {}
    with open(a.bench_paf) as f:
        for line in f:
            x = line.split("\t", 8)
            blen[x[5]] = int(x[6])
    truth = {}
    with open(a.sim_bed) as f:
        for line in f:
            x = line.rstrip("\n").split("\t")
            truth[x[3]] = (x[0], int(x[1]), int(x[2]), x[5])
    junc = {}
    with open(a.junctions) as f:
        next(f)
        for line in f:
            x = line.rstrip("\n").split("\t")
            junc[(x[1], x[2])] = int(x[5])
    gp = []
    with open(a.gpatch_bed) as f:
        for line in f:
            x = line.rstrip("\n").split("\t")
            gp.append((x[0], int(x[1]), int(x[2]), x[3], x[5], int(x[6]), int(x[7])))
    term = {}   # terminal patch length per chrom
    with open(a.patches_bed) as f:
        for line in f:
            c, s, e = line.split()[:3]
            if int(e) == rlen[c]:
                term[c] = int(e) - int(s)

    def call(d):
        return "trim" if d < 0 else "patch" if d > 0 else "bookend"

    def ttype(d):
        return None if d is None else "overlap" if d < 0 else "gap" if d > 0 else "abut"

    rows = []
    for i, (c, s, e, name, strand, qs, ql) in enumerate(gp):
        first = i == 0 or gp[i - 1][0] != c
        last = i == len(gp) - 1 or gp[i + 1][0] != c
        t = truth[name]
        if first:
            gd = -qs if qs > 0 else s
            td = t[1] if t[0] == c else None
            rows.append([c, "START", name, strand, "chrom_end" if td is not None else "diff_chrom",
                         ttype(td), td, call(gd), gd])
        else:
            pname = gp[i - 1][3]
            gd = -qs if qs > 0 else s - gp[i - 1][2]
            tl = truth[pname]
            if (pname, name) in junc:
                cls, td = "true_adjacent", junc[(pname, name)]
            elif (name, pname) in junc:
                cls, td = "true_reversed", junc[(name, pname)]
            elif tl[0] == t[0]:
                l, r = (tl, t) if tl[1] <= t[1] else (t, tl)
                cls, td = "not_adjacent", r[1] - l[2]
            else:
                cls, td = "diff_chrom", None
            rows.append([c, pname, name, strand, cls, ttype(td), td, call(gd), gd])
        if last:
            gd = term.get(c, 0)
            td = blen[t[0]] - t[2] if t[0] == c else None
            rows.append([c, name, "END", strand, "chrom_end" if td is not None else "diff_chrom",
                         ttype(td), td, call(gd), gd])
    for r in rows:
        r.append(None if r[6] is None else r[8] - r[6])

    with open(a.out, "w") as o:
        o.write("# compare_overlap_trimming.py %s %s\n" % (__version__, a.label))
        o.write("chrom\tprev_contig\tcontig\tstrand\tpair_class\ttrue_type\ttrue_dist"
                "\tgpatch_call\tgpatch_dist\tdiff\n")
        for r in rows:
            o.write("\t".join("NA" if v is None else str(v) for v in r) + "\n")

    w = sys.stdout.write
    w("==== %s ====\n" % a.label)
    w("placed contigs: %d; units: %s\n" % (len(gp), dict(Counter(r[4] for r in rows))))
    flip = [g[3] for g in gp if g[4] != truth[g[3]][3]]
    w("placed in the opposite orientation to the benchmark: %d %s\n" % (len(flip), ",".join(flip)))
    for scope, keep in (("all chromosomes", lambda r: True), ("excluding chrY", lambda r: r[0] != "chrY")):
        R = [r for r in rows if keep(r)]
        w("\n-- %s --\n" % scope)
        w("true type x call (true_adjacent | not_adjacent | chrom_end):\n")
        for cls in ("true_adjacent", "not_adjacent", "chrom_end"):
            cnt = Counter((r[5], r[7]) for r in R if r[4] == cls)
            w("  %-14s %s\n" % (cls, "  ".join("%s->%s:%d" % (k[0], k[1], v) for k, v in sorted(cnt.items()))))
        adj = [r for r in R if r[4] == "true_adjacent"]
        w("TRIMS   true overlap, trimmed      (diff = overlap - trim):   %s\n"
          % stats([r[9] for r in adj if r[5] == "overlap" and r[7] == "trim"]))
        w("        true overlap, any call     (diff = net bp error):     %s\n"
          % stats([r[9] for r in adj if r[5] == "overlap"]))
        w("PATCHES true gap, patched          (diff = patch - gap):      %s\n"
          % stats([r[9] for r in adj if r[5] == "gap" and r[7] == "patch"]))
        w("        true gap, any call         (diff = net bp error):     %s\n"
          % stats([r[9] for r in adj if r[5] == "gap"]))
        w("        not-adjacent pairs, any    (diff = net bp error):     %s\n"
          % stats([r[9] for r in R if r[4] == "not_adjacent"]))
        w("        chromosome ends, any       (diff = net bp error):     %s\n"
          % stats([r[9] for r in R if r[4] == "chrom_end"]))
        w("ALL units with truth (net bp error):                           %s\n"
          % stats([r[9] for r in R if r[9] is not None]))
    w("\n")


if __name__ == "__main__":
    main()
