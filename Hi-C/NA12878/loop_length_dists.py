#!/usr/bin/env python3
"""Loop length distributions by T2T-CHM13 sharing category, plus
scaffold-only loops absent from the patched assembly (and the reverse).

Usage: loop_length_dists.py <Hi-C/NA12878 dir> <output prefix>

Writes two tab-delimited tables:
  <prefix>.by_T2T_sharing.tsv       lengths by T2T-CHM13 sharing category
  <prefix>.scaffold_vs_patched.tsv  loops absent from / present in the other
                                    GPatch assembly (native coordinates)

Loop length is x1..y2 in native assembly coordinates. Requires bedtools on
PATH and loop_intervals.{5kb,10kb}.bed / loop_anchors.*.bed files in each
*_loops directory (see README.md)."""
import subprocess, statistics, bisect, sys
from collections import defaultdict

ROOT, PREFIX = sys.argv[1], sys.argv[2]
ASM = {"patched": f"{ROOT}/patched/juicer/NA12878_patched_loops",
       "scaffold-only": f"{ROOT}/scaffold-only/juicer/NA12878_patched_scaffold-only_loops"}
T2T = f"{ROOT}/T2T-CHM13/juicer/T2T-CHM13_loops"
BOUND = f"{ASM['scaffold-only']}/NA12878.2.cbreak_2.contig-boundaries.bed"
RES = {"5kb": 5000, "10kb": 10000}
CATS = ("fully_shared", "partially_shared", "unique")

STAT_COLS = ["n", "min", "q1", "median", "mean", "q3", "max",
             "n_ge1Mb", "frac_ge1Mb", "n_ge2Mb", "frac_ge2Mb",
             "n_span_contig_junction", "frac_span_contig_junction"]


def read_intervals(path):
    """id -> (chr, start, end); start = x1, end = y2."""
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


# contig boundaries per chrom (patched == scaffold-only coordinates)
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


def stats(spans, junctions=True):
    """spans: list of (chr, start, end). Returns values for STAT_COLS."""
    n = len(spans)
    if n == 0:
        return [0] + ["NA"] * (len(STAT_COLS) - 1)
    s = sorted(e - b for c, b, e in spans)
    q = statistics.quantiles(s, n=4) if n > 1 else [s[0]] * 3
    g1 = sum(x >= 1_000_000 for x in s)
    g2 = sum(x >= 2_000_000 for x in s)
    row = [n, s[0], round(q[0]), round(statistics.median(s)), round(statistics.mean(s)),
           round(q[2]), s[-1], g1, f"{g1/n:.4f}", g2, f"{g2/n:.4f}"]
    if junctions:
        nj = sum(spans_boundary(*x) for x in spans)
        row += [nj, f"{nj/n:.4f}"]
    else:
        row += ["NA", "NA"]
    return row


def same_loop_full(adir, res):
    """Loops whose two lifted anchors each overlap an anchor of the SAME T2T loop."""
    cmd = (f"bedtools intersect -a {adir}/loop_anchors.{res}.T2T-CHM13.bed "
           f"-b {T2T}/loop_anchors.{res}.bed -wa -wb")
    hits = defaultdict(lambda: [set(), set()])
    out = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True).stdout
    for l in out.splitlines():
        f = l.split("\t")
        lid, a = f[3].split("_")
        hits[lid][int(a) - 1].add(f[7].split("_")[0])
    return {lid for lid, (h1, h2) in hits.items() if h1 & h2}


def match_set(qry, ref, tol):
    """Query loops with a reference loop on the same chrom whose x1 and y1
    are each within tol bp."""
    idx = defaultdict(list)
    for c, x1, x2, y1, y2 in ref.values():
        idx[c].append((x1, y1))
    hit = set()
    for i, (c, x1, x2, y1, y2) in qry.items():
        if any(abs(rx - x1) <= tol and abs(ry - y1) <= tol for rx, ry in idx[c]):
            hit.add(i)
    return hit


# ---- Table 1: lengths by T2T-CHM13 sharing category ----
t1 = [["assembly", "resolution", "category"] + STAT_COLS]
for res in RES:
    iv = read_intervals(f"{T2T}/loop_intervals.{res}.bed")
    # contig junctions are defined in GPatch coordinates; not applicable to T2T-CHM13
    t1.append(["T2T-CHM13", res, "all"] + stats(list(iv.values()), junctions=False))

native, cats = {}, {}
for aname, adir in ASM.items():
    for res in RES:
        iv = read_intervals(f"{adir}/loop_intervals.{res}.bed")
        native[(aname, res)] = read_bedpe(f"{adir}/enriched_pixels_{RES[res]}.bedpe")
        cat = {i: "unique" for i in iv}
        for l in open(f"{adir}/loops.{res}.T2T-CHM13.shared.bed"):
            f = l.split()
            cat[f[6]] = "partially_shared" if f[0] == "." else "fully_shared"
        cats[(aname, res)] = cat
        strict = same_loop_full(adir, res)
        sel = lambda pred: [iv[i] for i in iv if pred(i)]
        t1.append([aname, res, "all"] + stats(list(iv.values())))
        for k in CATS:
            t1.append([aname, res, k] + stats(sel(lambda i: cat[i] == k)))
        t1.append([aname, res, "fully_shared_same_T2T_loop"] +
                  stats(sel(lambda i: cat[i] == "fully_shared" and i in strict)))
        t1.append([aname, res, "fully_shared_diff_T2T_loops"] +
                  stats(sel(lambda i: cat[i] == "fully_shared" and i not in strict)))

# ---- Table 2: scaffold-only vs patched, native coordinates ----
t2 = [["resolution", "query_assembly", "reference_assembly", "match_tolerance_bins",
       "group", "T2T_sharing_category"] + STAT_COLS]
for res, bs in RES.items():
    for q, r in (("scaffold-only", "patched"), ("patched", "scaffold-only")):
        Q, R, cat = native[(q, res)], native[(r, res)], cats[(q, res)]
        span = lambda i: (Q[i][0], Q[i][1], Q[i][4])
        for tol in (0, 1, 2):
            hit = match_set(Q, R, tol * bs)
            groups = {"absent_from_reference": [i for i in Q if i not in hit],
                      "present_in_reference": [i for i in Q if i in hit]}
            for g, ids in groups.items():
                t2.append([res, q, r, tol, g, "all"] + stats([span(i) for i in ids]))
                for k in CATS:
                    t2.append([res, q, r, tol, g, k] +
                              stats([span(i) for i in ids if cat[i] == k]))

for suffix, rows in (("by_T2T_sharing", t1), ("scaffold_vs_patched", t2)):
    path = f"{PREFIX}.{suffix}.tsv"
    with open(path, "w") as fh:
        for row in rows:
            fh.write("\t".join(map(str, row)) + "\n")
    sys.stderr.write(f"Wrote {path}\n")
