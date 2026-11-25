#!/usr/bin/env python3
"""Loop length distributions by T2T-CHM13 sharing category, plus
scaffold-only loops absent from the patched assembly.

Usage: loop_length_dists.py <Hi-C/NA12878 dir> <output.txt>
Requires bedtools on PATH and loop_intervals.{5kb,10kb}.bed / loop_anchors.*.bed
files in each *_loops directory (see README.md)."""
import subprocess, statistics, bisect, sys
from collections import defaultdict

ROOT = sys.argv[1]
ASM = {"patched": f"{ROOT}/patched/juicer/NA12878_patched_loops",
       "scaffold-only": f"{ROOT}/scaffold-only/juicer/NA12878_patched_scaffold-only_loops"}
T2T = f"{ROOT}/T2T-CHM13/juicer/T2T-CHM13_loops"
BOUND = f"{ASM['scaffold-only']}/NA12878.2.cbreak_2.contig-boundaries.bed"
RES = {"5kb": 5000, "10kb": 10000}
OUT = sys.argv[2]

def read_intervals(path):
    d = {}
    for l in open(path):
        c, s, e, i = l.split()[:4]
        d[i] = (c, int(s), int(e))
    return d

def read_bedpe(path):
    """Native loops: id -> (chr, x1, x2, y1, y2), id = data-line order."""
    d, n = {}, 0
    for l in open(path):
        f = l.split("\t")
        if f[0] == "chr1" or l.startswith("#"):
            continue
        n += 1
        d[f"loop{n}"] = ("chr" + f[0], int(f[1]), int(f[2]), int(f[4]), int(f[5]))
    return d

def summarize(lengths):
    n = len(lengths)
    if n == 0:
        return dict(n=0)
    s = sorted(lengths)
    q = statistics.quantiles(s, n=4) if n > 1 else [s[0]] * 3
    return dict(n=n, min=s[0], q1=q[0], median=statistics.median(s), mean=statistics.mean(s),
                q3=q[2], max=s[-1],
                ge1M=sum(x >= 1_000_000 for x in s), ge2M=sum(x >= 2_000_000 for x in s))

def fmt(name, st):
    if st["n"] == 0:
        return f"{name}\t0" + "\tNA" * 11
    return "\t".join(map(str, [name, st["n"], st["min"], int(st["q1"]), int(st["median"]),
                               round(st["mean"]), int(st["q3"]), st["max"],
                               st["ge1M"], f"{st['ge1M']/st['n']:.4f}",
                               st["ge2M"], f"{st['ge2M']/st['n']:.4f}"]))

HDR = "category\tn\tmin\tq1\tmedian\tmean\tq3\tmax\tn_ge1Mb\tfrac_ge1Mb\tn_ge2Mb\tfrac_ge2Mb"

def same_loop_full(adir, res):
    """Loops whose two lifted anchors each overlap an anchor of the SAME T2T loop."""
    cmd = f"bedtools intersect -a {adir}/loop_anchors.{res}.T2T-CHM13.bed -b {T2T}/loop_anchors.{res}.bed -wa -wb"
    hits = defaultdict(lambda: [set(), set()])
    for l in subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True).stdout.splitlines():
        f = l.split("\t")
        lid, a = f[3].split("_")
        hits[lid][int(a) - 1].add(f[7].split("_")[0])
    return {lid for lid, (h1, h2) in hits.items() if h1 & h2}

# contig boundaries per chrom (patched == scaffold coordinates)
bounds = defaultdict(list)
for l in open(BOUND):
    c, s = l.split()[:2]
    bounds[c].append(int(s))
for c in bounds:
    bounds[c].sort()

def spans_boundary(c, s, e):
    b = bounds.get(c, [])
    i = bisect.bisect_left(b, s)
    return i < len(b) and b[i] < e

out = open(OUT, "w")
def p(*a):
    print(*a); print(*a, file=out)

# Reference: all T2T-CHM13 loops
for res in RES:
    p(f"\n### T2T-CHM13, {res}: all loops")
    p(HDR)
    p(fmt("all", summarize([e - s for c, s, e in read_intervals(f"{T2T}/loop_intervals.{res}.bed").values()])))

native = {}
for aname, adir in ASM.items():
    for res in RES:
        iv = read_intervals(f"{adir}/loop_intervals.{res}.bed")
        native[(aname, res)] = read_bedpe(f"{adir}/enriched_pixels_{RES[res]}.bedpe")
        cat = {i: "unique" for i in iv}
        for l in open(f"{adir}/loops.{res}.T2T-CHM13.shared.bed"):
            f = l.split()
            cat[f[6]] = "partially_shared" if f[0] == "." else "fully_shared"
        strict = same_loop_full(adir, res)
        L = lambda ids: [iv[i][2] - iv[i][1] for i in ids]
        p(f"\n### {aname}, {res}: loops by T2T-CHM13 sharing (length = x1..y2 in native coords)")
        p(HDR)
        p(fmt("all", summarize(L(iv))))
        for k in ("fully_shared", "partially_shared", "unique"):
            p(fmt(k, summarize(L([i for i in iv if cat[i] == k]))))
        p(fmt("fully_shared_same_T2T_loop", summarize(L([i for i in iv if cat[i] == "fully_shared" and i in strict]))))
        p(fmt("fully_shared_diff_T2T_loops", summarize(L([i for i in iv if cat[i] == "fully_shared" and i not in strict]))))
        for k in ("fully_shared", "partially_shared", "unique"):
            ids = [i for i in iv if cat[i] == k]
            nb = sum(spans_boundary(*iv[i]) for i in ids)
            p(f"#   {k}: {nb}/{len(ids)} ({nb/len(ids):.3f}) span a contig junction")
        native[(aname, res, "cat")] = cat

# Scaffold-only loops not in patched (and the reverse), native coords
def match_set(qry, ref, tol):
    idx = defaultdict(list)
    for i, (c, x1, x2, y1, y2) in ref.items():
        idx[c].append((x1, y1))
    hit = set()
    for i, (c, x1, x2, y1, y2) in qry.items():
        for rx, ry in idx[c]:
            if abs(rx - x1) <= tol and abs(ry - y1) <= tol:
                hit.add(i); break
    return hit

for res, bs in RES.items():
    for q, r in (("scaffold-only", "patched"), ("patched", "scaffold-only")):
        Q, R = native[(q, res)], native[(r, res)]
        cat = native[(q, res, "cat")]
        p(f"\n### {res}: {q} loops NOT found in {r} (both anchors matched within tolerance)")
        p(HDR)
        for tol_bins in (0, 1, 2):
            hit = match_set(Q, R, tol_bins * bs)
            miss = [i for i in Q if i not in hit]
            lens = [Q[i][4] - Q[i][1] for i in miss]
            p(fmt(f"absent_tol{tol_bins}bin", summarize(lens)))
            if tol_bins == 0:
                p(fmt("present_tol0bin", summarize([Q[i][4] - Q[i][1] for i in hit])))
                for k in ("fully_shared", "partially_shared", "unique"):
                    p(fmt(f"absent_tol0bin_{k}", summarize([Q[i][4] - Q[i][1] for i in miss if cat[i] == k])))
                nb = sum(spans_boundary(Q[i][0], Q[i][1], Q[i][4]) for i in miss)
                nh = sum(spans_boundary(Q[i][0], Q[i][1], Q[i][4]) for i in hit)
                p(f"#   absent spanning contig junction: {nb}/{len(miss)} ({nb/len(miss):.3f}); present: {nh}/{len(hit)} ({nh/len(hit):.3f})")
out.close()
