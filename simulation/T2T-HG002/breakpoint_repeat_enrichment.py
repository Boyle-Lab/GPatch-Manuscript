#!/usr/bin/env python3
"""
Test whether simulated contig breakpoints are enriched in repeats and
segmental duplications of the benchmark genome (T2T-HG002 v1.1 paternal).

Two units of analysis:
  * endpoints  - every contig end that is not a chromosome end, classified as a
                 true model-contig terminus (the alignment reaches the end of the
                 HG002.1 contig, within --terminus-tol bp) or an internal
                 alignment break (the HG002.1 contig continues past the end).
                 Left ends are placed at `start`, right ends at `end - 1`.
  * junctions  - the interval between adjacent contigs (the gap, or the overlap
                 for overlapping pairs). One unit per junction, so the two
                 correlated contig ends at a junction are counted once.

For each annotation category and window size w, the statistic is the fraction
of units whose interval, padded by w on each side, intersects the category.
The null is obtained by placing units uniformly at random on the same
chromosome (preserving per-chromosome counts and, for junctions, interval
lengths) --n-perm times. Reports observed and expected fractions, fold
enrichment, two-sided empirical p-values and Benjamini-Hochberg q-values
(computed within each unit/window).
"""

import sys
import gzip
from argparse import ArgumentParser
from collections import defaultdict

import numpy as np

__version__ = "0.1.0"

SUFFIX = "_PATERNAL"


def open_any(f):
    return gzip.open(f, "rt") if f.endswith(".gz") else open(f)


def strip(chrom):
    return chrom[:-len(SUFFIX)] if chrom.endswith(SUFFIX) else chrom


def load_categories(a, chrom_lens):
    """Return {category: {chrom: list of (start, end)}}."""
    cats = defaultdict(lambda: defaultdict(list))

    def add(cat, c, s, e):
        if c in chrom_lens:
            cats[cat][c].append((s, e))

    # Segmental duplications (SEDEF/UW). The published set is already
    # filtered (>= 1 kb, >= 90% identity), so all records are used; high-
    # identity SDs are those with fracMatch (SEDEF identity) >= 0.98.
    with open_any(a.sd) as f:
        for line in f:
            if line.startswith("#"):
                ifr = line[1:].rstrip("\n").split("\t").index("fracMatch")
                continue
            x = line.rstrip("\n").split("\t")
            c, s, e = strip(x[0]), int(x[1]), int(x[2])
            add("SD (all)", c, s, e)
            if float(x[ifr]) >= 0.98:
                add("SD (>=98% identity)", c, s, e)

    # Centromeric satellite annotation (UCSC CenSat v2.0).
    alpha = {"active", "hor", "dhor", "mon", "mixedAlpha"}
    hsat = {"HSat1A", "HSat1B", "HSat2", "HSat3"}
    with open_any(a.censat) as f:
        for line in f:
            if line.startswith(("track", "#")):
                continue
            x = line.split("\t")
            lab = x[3].replace("(", "_").split("_")[0]
            if lab in ("ct", "GAP"):
                continue
            c, s, e = strip(x[0]), int(x[1]), int(x[2])
            add("CenSat (any)", c, s, e)
            if lab in alpha:
                add("CenSat alpha-satellite", c, s, e)
            elif lab in hsat:
                add("CenSat HSat1-3", c, s, e)
            elif lab == "rDNA":
                add("CenSat rDNA", c, s, e)
            else:
                add("CenSat other (beta/gamma/other)", c, s, e)

    # RepeatMasker (UConn) in UCSC "joined" bigRmsk format: name field is
    # "<repeat>#<class>/<family>"; chromStart/chromEnd span joined fragments
    # plus unaligned consensus flanks, so only the aligned blocks (relative
    # start != -1, size > 0) are used as genomic intervals.
    keep = {"SINE", "LINE", "LTR", "DNA", "Satellite", "Simple_repeat", "Low_complexity"}
    with open_any(a.rmsk) as f:
        for line in f:
            x = line.split("\t", 13)
            cls = x[3].split("#", 1)[-1].split("/", 1)[0]
            c, base = strip(x[0]), int(x[1])
            sizes = x[10].rstrip(",").split(",")
            rels = x[11].rstrip(",").split(",")
            for sz, rel in zip(sizes, rels):
                sz, rel = int(sz), int(rel)
                if rel < 0 or sz <= 0:
                    continue
                s, e = base + rel, base + rel + sz
                add("RepeatMasker (any)", c, s, e)
                if cls in keep:
                    add("RM " + cls, c, s, e)

    # TRF simple repeats (both haplotypes in file; keep paternal).
    with open_any(a.trf) as f:
        for line in f:
            x = line.split("\t", 3)
            if x[0].endswith(SUFFIX):
                add("TRF tandem repeats", strip(x[0]), int(x[1]), int(x[2]))

    # Microsatellites / homopolymer-rich (mnr10).
    with open_any(a.mnr) as f:
        for line in f:
            if line.startswith(("track", "#")):
                continue
            x = line.split("\t", 3)
            add("Microsatellites (mnr10)", strip(x[0]), int(x[1]), int(x[2]))

    # Union of everything.
    union = defaultdict(list)
    for cat in list(cats):
        for c, iv in cats[cat].items():
            union[c].extend(iv)
    cats["Any repeat or SD"] = union
    return cats


def merge(ivs):
    """Merge intervals; return (starts, ends) numpy arrays, sorted, disjoint."""
    if not ivs:
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    a = np.array(sorted(ivs), dtype=np.int64)
    s, e = a[:, 0], a[:, 1]
    run_end = np.maximum.accumulate(e)
    new = np.ones(len(s), bool)
    new[1:] = s[1:] > run_end[:-1]
    idx = np.flatnonzero(new)
    starts = s[idx]
    ends = np.append(run_end[idx[1:] - 1], run_end[-1])
    return starts, ends


def hits(starts, ends, qs, qe, w):
    """Boolean array: does [qs-w, qe+w) intersect any merged interval?"""
    if len(starts) == 0:
        return np.zeros(len(qs), bool)
    lo, hi = qs - w, qe + w
    j = np.searchsorted(ends, lo, side="right")
    ok = j < len(starts)
    out = np.zeros(len(qs), bool)
    out[ok] = starts[j[ok]] < hi[ok]
    return out


def bh(p):
    p = np.asarray(p)
    n = len(p)
    o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1)
    return out


def main():
    p = ArgumentParser(description="Breakpoint repeat/SD enrichment test.")
    p.add_argument("-b", "--bed", required=True, help="Simulated contigs BED (paf_to_simulated_contigs.py).")
    p.add_argument("-p", "--paf", required=True, help="PAF used to build the BED (for chromosome lengths).")
    p.add_argument("--sd", required=True)
    p.add_argument("--censat", required=True)
    p.add_argument("--rmsk", required=True)
    p.add_argument("--trf", required=True)
    p.add_argument("--mnr", required=True)
    p.add_argument("-w", "--windows", default="0,1000,10000",
                   help="Comma-separated padding sizes in bp (default 0,1000,10000).")
    p.add_argument("-n", "--n-perm", type=int, default=1000)
    p.add_argument("--terminus-tol", type=int, default=1000,
                   help="Max distance (bp) from a model-contig end for a breakpoint to count as a terminus (default 1000).")
    p.add_argument("--seed", type=int, default=20261109)
    p.add_argument("-o", "--out", required=True, help="Output TSV.")
    a = p.parse_args()
    windows = [int(x) for x in a.windows.split(",")]
    rng = np.random.default_rng(a.seed)

    chrom_lens = {}
    with open(a.paf) as f:
        for line in f:
            x = line.split("\t", 8)
            chrom_lens[x[5]] = int(x[6])

    # ---- Units -------------------------------------------------------------
    units = {k: defaultdict(list) for k in ("endpoints: all", "endpoints: model-contig terminus",
                                             "endpoints: internal alignment break", "junctions")}
    contigs = defaultdict(list)
    with open(a.bed) as f:
        for line in f:
            x = line.rstrip("\n").split("\t")
            c, s, e, strand = x[0], int(x[1]), int(x[2]), x[5]
            attr = dict(kv.split("=") for kv in x[7].split(";"))
            qs, qe, ql = int(attr["qs"]), int(attr["qe"]), int(attr["qlen"])
            contigs[c].append((s, e))
            left_q_end = qs <= a.terminus_tol if strand == "+" else qe >= ql - a.terminus_tol
            right_q_end = qe >= ql - a.terminus_tol if strand == "+" else qs <= a.terminus_tol
            for pos, is_term, at_chrom_end in ((s, left_q_end, s == 0),
                                               (e - 1, right_q_end, e == chrom_lens[c])):
                if at_chrom_end:
                    continue
                units["endpoints: all"][c].append((pos, pos + 1))
                key = "endpoints: model-contig terminus" if is_term else "endpoints: internal alignment break"
                units[key][c].append((pos, pos + 1))
    for c, ivs in contigs.items():
        for (s1, e1), (s2, e2) in zip(ivs, ivs[1:]):
            units["junctions"][c].append((min(e1, s2), max(e1, s2)))   # gap or overlap interval

    sys.stderr.write("units: %s\n" % {k: sum(len(v) for v in d.values()) for k, d in units.items()})

    # ---- Annotations -------------------------------------------------------
    cats = load_categories(a, chrom_lens)
    merged = {cat: {c: merge(iv) for c, iv in d.items()} for cat, d in cats.items()}
    genome = sum(chrom_lens.values())
    cov = {cat: sum(int((e - s).sum()) for s, e in d.values()) / genome for cat, d in merged.items()}
    empty = (np.zeros(0, np.int64), np.zeros(0, np.int64))
    cat_order = sorted(merged, key=lambda k: (k != "Any repeat or SD", k))

    # ---- Observed & permuted ------------------------------------------------
    rows = []
    for unit, d in units.items():
        chroms = sorted(d)
        qs = {c: np.array([s for s, _ in d[c]], np.int64) for c in chroms}
        ln = {c: np.array([e - s for s, e in d[c]], np.int64) for c in chroms}
        n = sum(len(qs[c]) for c in chroms)
        # Pre-draw permuted starts: per chrom, n_perm x n_c uniform in [0, L - len].
        perm = {c: (rng.random((a.n_perm, len(qs[c]))) * (chrom_lens[c] - ln[c] + 1)).astype(np.int64)
                for c in chroms}
        for w in windows:
            block = []
            for cat in cat_order:
                obs = 0
                null = np.zeros(a.n_perm)
                for c in chroms:
                    st, en = merged[cat].get(c, empty)
                    obs += hits(st, en, qs[c], qs[c] + ln[c], w).sum()
                    ps = perm[c]
                    pe = ps + ln[c][None, :]
                    h = hits(st, en, ps.ravel(), pe.ravel(), w).reshape(ps.shape)
                    null += h.sum(axis=1)
                obs_f, null_f = obs / n, null / n
                exp_f = null_f.mean()
                p_hi = (1 + (null_f >= obs_f).sum()) / (a.n_perm + 1)
                p_lo = (1 + (null_f <= obs_f).sum()) / (a.n_perm + 1)
                block.append([unit, w, cat, n, int(obs), obs_f, exp_f,
                              obs_f / exp_f if exp_f > 0 else float("nan"),
                              min(1.0, 2 * min(p_hi, p_lo)), cov[cat]])
            q = bh([r[8] for r in block])
            for r, qq in zip(block, q):
                rows.append(r + [qq])
            sys.stderr.write("done %s w=%d\n" % (unit, w))

    with open(a.out, "w") as o:
        o.write("# breakpoint_repeat_enrichment.py %s  n_perm=%d seed=%d terminus_tol=%d bed=%s\n"
                % (__version__, a.n_perm, a.seed, a.terminus_tol, a.bed))
        o.write("unit\twindow_bp\tcategory\tn_units\tn_hit\tobs_frac\texp_frac\tfold\tp_two_sided\tgenome_frac\tq_BH\n")
        for r in rows:
            o.write("%s\t%d\t%s\t%d\t%d\t%.4f\t%.4f\t%.3f\t%.4g\t%.4f\t%.4g\n" % tuple(r))


if __name__ == "__main__":
    main()
