#!/usr/bin/env python3
"""
Why GPatch overlap extension did not place more unplaced simulated contigs.

1. Walk stops (extension.tsv): for each walk that stopped with "no candidate",
   the true benchmark neighbour of the walk's last contig, in the walk
   direction, is classified (already placed, unmapped, separated by a true gap
   or < --min-overlap overlap, or overlapping but excluded by GPatch's
   candidate rules: other chromosome, outside the gap region, opposite
   orientation).
2. Reachability: unplaced contigs reachable from a placed flanking contig of
   the same benchmark contig through a chain of true overlaps >= --min-overlap,
   versus contigs cut off by gaps or short overlaps.
3. Islands: runs of unplaced contigs joined by true overlaps >= --min-overlap.
"""
import sys
import statistics as st
from argparse import ArgumentParser
from collections import defaultdict, Counter
import pysam

__version__ = "0.1.0"


def main():
    p = ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("-s", "--sim-bed", required=True)
    p.add_argument("-j", "--junctions", required=True)
    p.add_argument("-g", "--gpatch-bed", required=True)
    p.add_argument("-u", "--unplaced", required=True)
    p.add_argument("-e", "--extension", required=True)
    p.add_argument("-b", "--bam", required=True)
    p.add_argument("--region-pad", type=int, default=1000000, help="GPatch --ext_region_pad (default 1 Mb).")
    p.add_argument("--min-overlap", type=int, default=5000, help="GPatch --ext_min_overlap (default 5 kb).")
    a = p.parse_args()
    flip = {"+": "-", "-": "+"}

    truth, byb = {}, defaultdict(list)
    for l in open(a.sim_bed):
        x = l.split("\t")
        truth[x[3]] = (x[0], int(x[1]), int(x[2]), x[5])
        byb[x[0]].append(x[3])
    for b in byb:
        byb[b].sort(key=lambda n: (truth[n][1], -truth[n][2]))
    pos = {n: (b, i) for b, L in byb.items() for i, n in enumerate(L)}
    gp = {l.split("\t")[3]: l.split("\t")[5] for l in open(a.gpatch_bed)}
    unpl = {l.split("\t")[0]: l.split("\t")[2] for l in list(open(a.unplaced))[1:]}
    J = {}
    for l in open(a.junctions):
        if not l.startswith("chrom"):
            x = l.split("\t")
            J[(x[1], x[2])] = int(x[5])
    prim = {}
    for r in pysam.AlignmentFile(a.bam):
        if not (r.is_unmapped or r.is_secondary or r.is_supplementary):
            prim[r.query_name] = (r.reference_name, r.reference_start, r.reference_end, "-" if r.is_reverse else "+")

    # 1. walk stops
    res = Counter()
    for l in open(a.extension):
        if l.startswith(("#", "chrom\t")):
            continue
        x = l.rstrip("\n").split("\t")
        if x[5] != "stop":
            continue
        chrom, A, B, d, X, why = x[0], x[2], x[3], x[4], x[6], x[9]
        if why.startswith("reached"):
            res["reached the flanking contig"] += 1
            continue
        if why.startswith("ambiguous"):
            res["ambiguous (several candidates)"] += 1
            continue
        if X not in gp:
            res["no candidate: walk end rejected later (claimed by several walks)"] += 1
            continue
        b, i = pos[X]
        same = gp[X] == truth[X][3]
        j = i + 1 if (same if d == "fwd" else not same) else i - 1
        if not 0 <= j < len(byb[b]):
            res["no candidate: end of benchmark contig"] += 1
            continue
        Y = byb[b][j]
        ov = -J[(X, Y)] if (X, Y) in J else -J[(Y, X)]
        if Y in gp:
            c = "true neighbour already placed (%s)" % ("gap or overlap < min" if ov < a.min_overlap else "overlap >= min")
        elif unpl.get(Y) == "unmapped":
            c = "true neighbour unmapped"
        elif ov < a.min_overlap:
            c = "true neighbour unplaced, separated by a gap or overlap < min"
        else:
            q = prim[Y]
            lo = (prim[A][1] - a.region_pad) if A != "START" else 0
            hi = (prim[B][2] + a.region_pad) if B != "END" else float("inf")
            need = truth[Y][3] if same else flip[truth[Y][3]]
            if q[0] != chrom:
                c = "true neighbour overlaps, but aligned to another chromosome"
            elif not (q[1] >= lo and q[2] <= hi):
                c = "true neighbour overlaps, but aligned outside the gap region"
            elif q[3] != need:
                c = "true neighbour overlaps, but aligned in the opposite orientation"
            else:
                c = "true neighbour overlaps and is eligible (overlap not detected)"
        res["no candidate: " + c] += 1
    print("== Extension walk stops (%d)" % sum(res.values()))
    for k, v in res.most_common():
        print("%6d  %s" % (v, k))

    # 2. reachability and 3. islands
    cat, bp, isl = Counter(), Counter(), []
    for b, L in byb.items():
        i = 0
        while i < len(L):
            if L[i] in gp:
                i += 1
                continue
            j = i
            while j < len(L) and L[j] not in gp:
                j += 1
            run = L[i:j]
            left, right = (L[i - 1] if i > 0 else None), (L[j] if j < len(L) else None)
            reach, prev = set(), left
            for n in run:
                if prev is None or -J[(prev, n)] < a.min_overlap:
                    break
                reach.add(n)
                prev = n
            nxt = right
            for n in reversed(run):
                if nxt is None or -J[(n, nxt)] < a.min_overlap:
                    break
                reach.add(n)
                nxt = n
            for n in run:
                k = ("reachable from a placed contig by an overlap chain" if n in reach else
                     "cut off by gaps or short overlaps" if (left or right) else
                     "in a benchmark contig with no placed contig")
                cat[k] += 1
                bp[k] += truth[n][2] - truth[n][1]
            cur = [run[0]]
            for p_, n in zip(run, run[1:]):
                if -J[(p_, n)] >= a.min_overlap:
                    cur.append(n)
                else:
                    isl.append(cur)
                    cur = [n]
            isl.append(cur)
            i = j
    print("\n== Unplaced contigs by reachability through true overlaps >= %d bp" % a.min_overlap)
    for k in cat:
        print("%6d contigs %7.2f Mb  %s" % (cat[k], bp[k] / 1e6, k))
    sz = [len(i) for i in isl]
    span = [truth[i[-1]][2] - truth[i[0]][1] for i in isl]
    print("\n== Islands of unplaced contigs joined by true overlaps >= %d bp: %d (singletons %d; max %d contigs)"
          % (a.min_overlap, len(isl), sum(s == 1 for s in sz), max(sz)))
    print("   span median %.1f kb; >= 50 kb: %d; >= 100 kb: %d (%.1f Mb)"
          % (st.median(span) / 1e3, sum(s >= 5e4 for s in span), sum(s >= 1e5 for s in span),
             sum(s for s in span if s >= 1e5) / 1e6))
    c = Counter()
    for b, L in byb.items():
        for x_, y_ in zip(L, L[1:]):
            if x_ not in gp and y_ not in gp:
                d = J[(x_, y_)]
                c["overlap >= min" if d <= -a.min_overlap else "overlap < min" if d < 0 else "gap/abut"] += 1
    print("   junctions between consecutive unplaced contigs: %s" % dict(c))


if __name__ == "__main__":
    main()
