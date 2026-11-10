#!/usr/bin/env python3
"""
Prototype of revised GPatch junction handling (options 1 + 2).

Takes GPatch's contig placement (order, strand, chromosome) as given and
re-decides only how adjacent placed contigs are joined. No N characters are
ever introduced: every junction is a trim, a reference patch, or a bookend.

Option 1 - contig extent on the reference:
  The primary alignment is chained to supplementary alignments on the same
  chromosome and strand that are colinear with it in both query and
  reference (MAPQ >= --chain-mapq, aligned query length >= --chain-minlen,
  gaps <= --chain-maxgap). Each contig end is then extrapolated 1:1 only by the
  residual unaligned query beyond the outermost chained alignment (GPatch
  extrapolates the full soft clip of the primary alignment).

Option 2 - direct overlap measurement (anchored dovetail test):
  For every junction, candidate overlaps between the last L bp of the previous
  contig (as placed) and the first L bp of the next contig come from exact
  matches of two --anchor bp anchors: the next contig's first bases found in
  the previous contig's tail, and the previous contig's last bases found in
  the next contig's head. Candidates supported by both anchors at the same
  offset (or within --indel-tol, verified with minimap2) are accepted if the
  whole overlap has identity >= --min-id. The trim is the position in the next
  contig that corresponds to the end of the previous contig. L starts at
  max(--win-min, 2 x estimated overlap + residual clips) and is doubled (up to
  --win-max) if an anchor suggests the overlap extends beyond the window.
  (A general aligner is not used to find dovetails: in palindromes and tandem
  arrays the short true overlap is reported as a low-scoring secondary of a
  longer paralogous alignment and filtered out.)

Decision per junction. The option 1 estimate is "trusted" when the residual
clips at both contig ends of the junction are <= --max-resid, "one-sided" when
only one end exceeds it, and "untrusted" when both do.
  trim_measured  one dovetail overlap found -> trim it
  trim_ambiguous several distinct dovetail overlaps (e.g. tandem repeats)
                 -> if the estimate is trusted or one-sided, trim the
                 candidate closest to the estimated overlap; otherwise trim
                 the smallest (never over-trims if the true one is among them)
  patch          no overlap; estimate >= 0 and trusted -> reference patch
  patch_one_sided no overlap; estimate >= 0 and one-sided -> reference patch
  bookend        no overlap and the estimate is negative or untrusted
                 -> join contigs directly
Chromosome ends are patched from the reference if the option 1 estimate is
trustworthy (residual <= --max-resid and positive extent); contigs are never
trimmed against a chromosome end.

Outputs (<prefix>.*): contigs.bed (GPatch format, patched-genome frame),
patches.bed (reference frame), junctions.tsv (decision and reasons),
patched.fasta (only with --write-fasta).
"""

import os
import bisect
import sys
import subprocess
import tempfile
from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor
from collections import defaultdict

import numpy as np
import pysam

__version__ = "0.3.0"

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s):
    return s.translate(COMP)[::-1]


class Aln:
    __slots__ = ("chrom", "rs", "re", "qs", "qe", "qlen", "strand", "mapq", "primary")

    def __init__(self, r):
        c = r.cigartuples
        lead = 0
        for op, n in c:
            if op in (4, 5):
                lead += n
            else:
                break
        self.chrom = r.reference_name
        self.rs, self.re = r.reference_start, r.reference_end
        self.qlen = r.infer_read_length()
        self.qs = lead          # placement-orientation query coords (BAM orientation)
        self.qe = lead + r.query_alignment_length
        self.strand = "-" if r.is_reverse else "+"
        self.mapq = r.mapping_quality
        self.primary = not (r.is_secondary or r.is_supplementary)


def chain_extent(p, supps, a):
    """Option 1: chain colinear supplementaries; return extent info."""
    cands = [s for s in supps if s.chrom == p.chrom and s.strand == p.strand
             and s.mapq >= a.chain_mapq and s.qe - s.qs >= a.chain_minlen]
    left = right = p
    nl = nr = 0
    while True:
        best = None
        for s in cands:
            dq, dr = left.qs - s.qe, left.rs - s.re
            if s.qs < left.qs and s.rs < left.rs and -a.chain_tol <= dq <= a.chain_maxgap \
                    and -a.chain_tol <= dr <= a.chain_maxgap:
                if best is None or s.qs < best.qs:
                    best = s
        if best is None:
            break
        left, nl = best, nl + 1
    while True:
        best = None
        for s in cands:
            dq, dr = s.qs - right.qe, s.rs - right.re
            if s.qe > right.qe and s.re > right.re and -a.chain_tol <= dq <= a.chain_maxgap \
                    and -a.chain_tol <= dr <= a.chain_maxgap:
                if best is None or s.qe > best.qe:
                    best = s
        if best is None:
            break
        right, nr = best, nr + 1
    return dict(ext_start=left.rs - left.qs, ext_end=right.re + (p.qlen - right.qe),
                resid5=left.qs, resid3=p.qlen - right.qe,
                gp_start=p.rs - p.qs, gp_end=p.re + (p.qlen - p.qe),
                clip5=p.qs, clip3=p.qlen - p.qe, nchain=1 + nl + nr)


def find_all(hay, needle):
    out, i = [], hay.find(needle)
    while i != -1:
        out.append(i)
        i = hay.find(needle, i + 1)
    return out


def identity_equal(x, y):
    """Fraction of identical bases between equal-length strings."""
    ax = np.frombuffer(x.encode("ascii"), np.uint8)
    ay = np.frombuffer(y.encode("ascii"), np.uint8)
    return float((ax == ay).mean())


def identity_mm2(x, y, a, tag):
    """Identity of a near-equal-length pair via minimap2 (indel-containing overlaps);
    requires the alignment to cover >= 99% of both sequences."""
    with open(tag + ".q.fa", "w") as o:
        o.write(">q\n%s\n" % x)
    with open(tag + ".t.fa", "w") as o:
        o.write(">t\n%s\n" % y)
    out = subprocess.run([a.minimap2, "-x", "asm5", "-c", "-t", "1", tag + ".t.fa", tag + ".q.fa"],
                         capture_output=True, text=True, check=True).stdout
    best = 0.0
    for line in out.splitlines():
        f = line.split("\t")
        qlen, qs, qe, st, tlen, ts, te, nm, bl = (int(f[1]), int(f[2]), int(f[3]), f[4],
                                                  int(f[6]), int(f[7]), int(f[8]), int(f[9]), int(f[10]))
        if st == "+" and (qe - qs) >= 0.99 * qlen and (te - ts) >= 0.99 * tlen:
            best = max(best, nm / bl)
    return best


def dovetail(fa, prev, cur, L0, a, tmp):
    """Option 2: anchored suffix(prev)/prefix(cur) overlap detection.

    Exact matches of the first --anchor bp of `cur` in the tail window of
    `prev` give candidate overlaps ov1 = lp - i; exact matches of the last
    --anchor bp of `prev` in the head window of `cur` give ov2 = j + anchor.
    Candidates where ov1 == ov2 are verified by direct base comparison;
    pairs with |ov1 - ov2| <= --indel-tol are verified with minimap2. A
    candidate is accepted if identity over the whole overlap >= --min-id.
    Returns (sorted trims, window, note)."""
    L, k = L0, a.anchor
    note = ""
    tag = os.path.join(tmp, "%s__%s" % (prev["name"], cur["name"]))
    while True:
        lp, lc = min(L, prev["len"]), min(L, cur["len"])
        ps = fetch(fa, prev, prev["len"] - lp, prev["len"])
        cs = fetch(fa, cur, 0, lc)
        ov1 = sorted({lp - i for i in find_all(ps, cs[:k])})
        ov2 = sorted({j + k for j in find_all(cs, ps[-k:])})
        s2 = set(ov2)
        trims, n_tested = [], 0
        for o in ov1:
            if o in s2:
                n_tested += 1
                if identity_equal(ps[lp - o:], cs[:o]) >= a.min_id:
                    trims.append(o)
        if not trims and a.indel_tol > 0:
            # near-consistent anchor pairs (indels inside the overlap)
            near = []
            for o1 in ov1:
                lo = bisect.bisect_left(ov2, o1 - a.indel_tol)
                hi = bisect.bisect_right(ov2, o1 + a.indel_tol)
                near.extend((abs(o1 - o2), o1, o2) for o2 in ov2[lo:hi] if o2 != o1)
            for _, o1, o2 in sorted(near)[:a.max_mm2_tests]:
                n_tested += 1
                if identity_mm2(ps[lp - o1:], cs[:o2], a, tag) >= a.min_id:
                    trims.append(o2)   # trim = position in cur of prev's end
        # overlap may extend beyond the window: cur's prefix anchors at the very
        # start of the prev window, or prev's suffix anchors at the end of cur's
        if not trims and L < a.win_max and ((ov1 and ov1[-1] > lp - k and lp < prev["len"]) or
                                            (ov2 and ov2[-1] > lc - k and lc < cur["len"])):
            L = min(2 * L, a.win_max)
            note = "window doubled"
            continue
        trims = [t for t in sorted(set(trims)) if t < cur["len"]]
        note = ("%s; " % note if note else "") + "anchors %d/%d, tested %d" % (len(ov1), len(ov2), n_tested)
        return trims, L, note


def dovetail_job(assembly, prev, cur, L0, a, tmp):
    return dovetail(pysam.FastaFile(assembly), prev, cur, L0, a, tmp)


def fetch(fa, c, s, e):
    """Sequence [s, e) of contig c in placement orientation."""
    if c["strand"] == "+":
        return fa.fetch(c["name"], s, e).upper()
    n = c["len"]
    return rc(fa.fetch(c["name"], n - e, n - s).upper())


def main():
    p = ArgumentParser(description="Prototype GPatch junction handling (options 1+2).")
    p.add_argument("-b", "--bam", required=True, help="GPatch input BAM (assembly vs reference).")
    p.add_argument("-c", "--contigs-bed", required=True, help="GPatch contigs.bed (placement order/strand).")
    p.add_argument("-a", "--assembly", required=True, help="Assembly FASTA (faidx-indexed).")
    p.add_argument("-r", "--reference", required=True, help="Reference FASTA (faidx-indexed).")
    p.add_argument("-x", "--prefix", required=True)
    p.add_argument("--write-fasta", action="store_true")
    p.add_argument("--chain-mapq", type=int, default=20)
    p.add_argument("--chain-minlen", type=int, default=10000)
    p.add_argument("--chain-maxgap", type=int, default=1000000)
    p.add_argument("--chain-tol", type=int, default=10000,
                   help="Allowed overlap between chained alignments (bp).")
    p.add_argument("--max-resid", type=int, default=10000,
                   help="Max residual clip (bp, per contig end) for a trustworthy gap estimate.")
    p.add_argument("--win-min", type=int, default=1000000)
    p.add_argument("--win-max", type=int, default=4000000)
    p.add_argument("--anchor", type=int, default=500,
                   help="Exact anchor length (bp); also the minimum detectable overlap.")
    p.add_argument("--indel-tol", type=int, default=50,
                   help="Max |ov1-ov2| for anchor pairs verified with minimap2 (0 disables).")
    p.add_argument("--max-mm2-tests", type=int, default=20,
                   help="Max near-consistent anchor pairs verified with minimap2 per junction.")
    p.add_argument("--min-id", type=float, default=0.999,
                   help="Min identity over the whole candidate overlap.")
    p.add_argument("--minimap2", default="minimap2")
    p.add_argument("-j", "--jobs", type=int, default=12)
    p.add_argument("--tmpdir", default=None)
    a = p.parse_args()

    order = defaultdict(list)
    names = set()
    with open(a.contigs_bed) as f:
        for line in f:
            x = line.rstrip("\n").split("\t")
            order[x[0]].append((x[3], x[5]))
            names.add(x[3])

    prim, supps = {}, defaultdict(list)
    for r in pysam.AlignmentFile(a.bam):
        if r.is_unmapped or r.is_secondary or r.query_name not in names:
            continue
        al = Aln(r)
        (prim.__setitem__(r.query_name, al) if al.primary else supps[r.query_name].append(al))

    fa = pysam.FastaFile(a.assembly)
    ref = pysam.FastaFile(a.reference)
    rlen = dict(zip(ref.references, ref.lengths))

    info = {}
    for chrom, lst in order.items():
        for name, strand in lst:
            pr = prim[name]
            assert pr.chrom == chrom and pr.strand == strand, name
            d = chain_extent(pr, supps[name], a)
            d.update(name=name, strand=strand, len=fa.get_reference_length(name))
            info[name] = d

    # ---- Option 2 for every adjacent pair (in parallel) ---------------------
    pairs = []
    for chrom, lst in order.items():
        for (pn, _), (cn, _) in zip(lst, lst[1:]):
            pairs.append((chrom, pn, cn))
    # Running position (GPatch's `pos`) only depends on extents, not decisions.
    tmp = tempfile.mkdtemp(prefix="gpj_", dir=a.tmpdir)
    jobs = {}
    with ProcessPoolExecutor(a.jobs) as ex:
        for chrom, pn, cn in pairs:
            P, C = info[pn], info[cn]
            est_ov = P["ext_end"] - C["ext_start"]
            L0 = min(a.win_max, max(a.win_min, 2 * max(est_ov, 0) + P["resid3"] + C["resid5"]))
            jobs[(pn, cn)] = ex.submit(dovetail_job, a.assembly, P, C, L0, a, tmp)
        meas = {k: v.result() for k, v in jobs.items()}
    for fn in os.listdir(tmp):
        os.remove(os.path.join(tmp, fn))
    os.rmdir(tmp)

    # ---- Build patched genome ----------------------------------------------
    cb = open(a.prefix + ".contigs.bed", "w")
    pb = open(a.prefix + ".patches.bed", "w")
    jt = open(a.prefix + ".junctions.tsv", "w")
    pf = open(a.prefix + ".patched.fasta", "w") if a.write_fasta else None
    jt.write("# gpatch_junction_prototype.py %s %s\n" % (__version__, " ".join(sys.argv[1:])))
    jt.write("chrom\tprev_contig\tcontig\tprev_ext_end\tcontig_ext_start\tprev_resid3\tcontig_resid5"
             "\tprev_chain_n\tcontig_chain_n\test_dist_opt1\test_dist_gpatch\twindow_bp\tdovetail_trims"
             "\tdecision\ttrim_bp\tpatch_ref_start\tpatch_ref_end\tpatch_bp\treason\n")

    def jrow(*v):
        jt.write("\t".join(str(x) for x in v) + "\n")

    for chrom in ref.references:
        if chrom not in order:
            continue
        L = rlen[chrom]
        seq = []
        plen = 0
        pos = None   # reference position of the end of the last placed contig
        prev = None
        for name, strand in order[chrom]:
            C = info[name]
            trim, patch = 0, None
            if prev is None:
                s = C["ext_start"]
                if s > 0 and C["resid5"] <= a.max_resid:
                    patch, dec, why = (0, min(s, L)), "patch", "chromosome start; option 1 extent trusted"
                elif s > 0:
                    dec, why = "bookend", "chromosome start; residual 5' clip %d > max" % C["resid5"]
                else:
                    dec, why = "bookend", "chromosome start; extrapolated start %d <= 0 (no trim at chrom end)" % s
                jrow(chrom, "START", name, 0, s, ".", C["resid5"], ".", C["nchain"], s, C["gp_start"],
                     ".", ".", dec, 0, *(patch if patch else (".", ".")),
                     (patch[1] - patch[0]) if patch else 0, why)
            else:
                P = info[prev]
                trims, win, note = meas[(prev, name)]
                est = C["ext_start"] - pos
                gp_est = C["gp_start"] - P["gp_end"]
                why_c = ("; " + note) if note else ""
                big = [x for x, v in (("prev residual 3' clip %d" % P["resid3"], P["resid3"]),
                                      ("residual 5' clip %d" % C["resid5"], C["resid5"])) if v > a.max_resid]
                trust = ("trusted", "one-sided", "untrusted")[len(big)]
                if len(trims) == 1:
                    dec, trim, why = "trim_measured", trims[0], "single dovetail overlap"
                elif len(trims) > 1 and trust != "untrusted":
                    trim = min(trims, key=lambda t: (abs(t + est), t))
                    dec, why = "trim_ambiguous", ("%d dovetail offsets; closest to %s estimated overlap %d used"
                                                  % (len(trims), trust, -est))
                elif len(trims) > 1:
                    dec, trim, why = "trim_ambiguous", trims[0], ("%d dovetail offsets; estimate untrusted, smallest used"
                                                                 % len(trims))
                elif est >= 0 and trust == "trusted":
                    dec, why = "patch", "no overlap; option 1 gap estimate trusted"
                    patch = (min(pos, L), min(C["ext_start"], L))
                elif est >= 0 and trust == "one-sided":
                    dec, why = "patch_one_sided", "no overlap; one-sided gap estimate used (%s)" % big[0]
                    patch = (min(pos, L), min(C["ext_start"], L))
                else:
                    dec = "bookend"
                    rs = (["estimate predicts overlap (%d) but none found" % -est] if est < 0 else []) + big
                    why = "no overlap; unusable gap estimate: " + ", ".join(rs)
                why += why_c
                jrow(chrom, prev, name, pos, C["ext_start"], P["resid3"], C["resid5"], P["nchain"], C["nchain"],
                     est, gp_est, win, ",".join(map(str, trims)) or ".", dec, trim,
                     *(patch if patch else (".", ".")), (patch[1] - patch[0]) if patch else 0, why)
            if patch and patch[1] > patch[0]:
                pb.write("%s\t%d\t%d\n" % (chrom, patch[0], patch[1]))
                plen += patch[1] - patch[0]
                if pf:
                    seq.append(ref.fetch(chrom, patch[0], patch[1]).lower())
            start = plen
            plen += C["len"] - trim
            if pf:
                seq.append(fetch(fa, C, trim, C["len"]))
            cb.write("%s\t%d\t%d\t%s\t.\t%s\t%d\t%d\n" % (chrom, start, plen, name, strand, trim, C["len"]))
            e = C["ext_end"]
            pos = e if pos is None else max(pos, e)
            prev = name
        # chromosome end
        P = info[prev]
        if pos < L and P["resid3"] <= a.max_resid:
            dec, why, patch = "patch", "chromosome end; option 1 extent trusted", (max(pos, 0), L)
        elif pos < L:
            dec, why, patch = "bookend", "chromosome end; residual 3' clip %d > max" % P["resid3"], None
        else:
            dec, why, patch = "bookend", "chromosome end; extrapolated end %d >= chrom length" % pos, None
        jrow(chrom, prev, "END", pos, L, P["resid3"], ".", P["nchain"], ".", L - pos, L - P["gp_end"],
             ".", ".", dec, 0, *(patch if patch else (".", ".")), (patch[1] - patch[0]) if patch else 0, why)
        if patch:
            pb.write("%s\t%d\t%d\n" % (chrom, patch[0], patch[1]))
            if pf:
                seq.append(ref.fetch(chrom, patch[0], patch[1]).lower())
        if pf:
            s = "".join(seq)
            pf.write(">%s\n" % chrom)
            for i in range(0, len(s), 80):
                pf.write(s[i:i + 80] + "\n")
    for fh in (cb, pb, jt) + ((pf,) if pf else ()):
        fh.close()


if __name__ == "__main__":
    main()
