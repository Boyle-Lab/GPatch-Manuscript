#!/usr/bin/env python3
"""
Prototype: place unplaced contigs by extending from placed contigs through
unambiguous sequence overlaps.

Starting from a GPatch placement (contigs.bed: order and strand along each
reference chromosome) and its unplaced-contig report, each gap between
adjacent placed contigs A and B (and each chromosome end) is extended:

  forward from A:  find unplaced contigs X whose 5' end overlaps the 3' end of
                   the current contig (anchored suffix/prefix test, as used for
                   junction handling); if exactly one such contig qualifies,
                   append it and continue from it.
  backward from B: likewise, unplaced contigs whose 3' end overlaps the 5' end
                   of the current contig, prepended.

A candidate qualifies only if its primary alignment is on the same reference
chromosome, it is used in the orientation of that alignment, and the
alignment lies within the gap region (from --region-pad bp before A's
alignment to --region-pad bp after B's). A walk stops when it reaches the
other flanking contig (an overlap with it is found), when no candidate or
more than one candidate qualifies, or after --max-steps steps. Contigs claimed
by more than one walk are rejected and walks are truncated before them.

Output: an augmented placement BED (same format as GPatch contigs.bed; only
chrom, name and strand are meaningful) for gpatch_junction_prototype.py, and
a walk log describing every step and stop reason.
"""

import sys
from argparse import ArgumentParser
from collections import defaultdict, Counter

import pysam

__version__ = "0.1.0"

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s):
    return s.translate(COMP)[::-1]


def find_all(hay, needle):
    out, i = [], hay.find(needle)
    while i != -1:
        out.append(i)
        i = hay.find(needle, i + 1)
    return out


def identity_at_least(x, y, min_id, chunk=1024):
    max_mm = int((1.0 - min_id) * len(x))
    mm = 0
    for i in range(0, len(x), chunk):
        a, b = x[i:i + chunk], y[i:i + chunk]
        if a != b:
            mm += sum(1 for c, d in zip(a, b) if c != d)
            if mm > max_mm:
                return False
    return True


class Seqs:
    """Head/tail windows of contigs in placement orientation (+ = as stored)."""

    def __init__(self, fasta, window):
        self.fa = pysam.FastaFile(fasta)
        self.w = window
        self.cache = {}

    def ends(self, name, strand):
        k = (name, strand)
        if k not in self.cache:
            L = self.fa.get_reference_length(name)
            h = self.fa.fetch(name, 0, min(self.w, L)).upper()
            t = self.fa.fetch(name, max(0, L - self.w), L).upper()
            self.cache[k] = (h, t, L) if strand == "+" else (rc(t), rc(h), L)
        return self.cache[k]


def overlaps(seqs, a, b, k, min_id, min_ovl=0):
    """Overlap lengths (>= min_ovl) where the 3' end of oriented contig a overlaps the 5' end of b."""
    _, ta, la = seqs.ends(*a)
    hb, _, lb = seqs.ends(*b)
    lp, lc = len(ta), len(hb)
    if lp < k or lc < k:
        return []
    ov1 = set(lp - i for i in find_all(ta, hb[:k]))
    if not ov1:
        return []
    ov2 = set(j + k for j in find_all(hb, ta[-k:]))
    return [o for o in sorted(ov1 & ov2)
            if min_ovl <= o < min(la, lb) and identity_at_least(ta[lp - o:], hb[:o], min_id)]


def main():
    p = ArgumentParser(description="Prototype overlap-based extension of GPatch placements.")
    p.add_argument("-c", "--contigs-bed", required=True, help="GPatch contigs.bed (placement).")
    p.add_argument("-u", "--unplaced", required=True, help="GPatch unplaced.tsv.")
    p.add_argument("-b", "--bam", required=True, help="GPatch input BAM (primary alignments).")
    p.add_argument("-a", "--assembly", required=True, help="Assembly FASTA (faidx-indexed).")
    p.add_argument("-r", "--ref-fai", required=True, help="Reference .fai (chromosome lengths).")
    p.add_argument("-o", "--out", required=True, help="Augmented placement BED.")
    p.add_argument("-l", "--log", required=True, help="Walk log TSV.")
    p.add_argument("-e", "--extents", required=True,
                   help="Output TSV of inserted contigs with the contig and overlap they were reached from, "
                        "for gpatch_junction_prototype.py --propagated (reference extents are carried through "
                        "overlaps instead of taken from the inserted contigs' own alignments).")
    p.add_argument("--anchor", type=int, default=500)
    p.add_argument("--min-id", type=float, default=0.999)
    p.add_argument("--window", type=int, default=2000000)
    p.add_argument("--region-pad", type=int, default=1000000)
    p.add_argument("--max-steps", type=int, default=50)
    p.add_argument("--min-overlap", type=int, default=5000,
                   help="Minimum overlap (bp) for an extension step or for closing a gap; short "
                        "overlaps in tandem/palindromic repeats are often spurious. Default 5000.")
    a = p.parse_args()

    rlen = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(a.ref_fai)}
    order = defaultdict(list)
    for line in open(a.contigs_bed):
        x = line.rstrip("\n").split("\t")
        order[x[0]].append((x[3], x[5]))
    unplaced = {}
    with open(a.unplaced) as f:
        next(f)
        for line in f:
            x = line.rstrip("\n").split("\t")
            unplaced[x[0]] = x[2]
    aln = {}
    for r in pysam.AlignmentFile(a.bam):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        aln[r.query_name] = (r.reference_name, r.reference_start, r.reference_end, "-" if r.is_reverse else "+")
    seqs = Seqs(a.assembly, a.window)

    cands = defaultdict(list)   # chrom -> [(name, strand, rs, re)]
    for n in unplaced:
        if n in aln:
            c, s, e, st = aln[n]
            cands[c].append((n, st, s, e))

    log = []
    ovl = {}     # (chrom, gap, direction, name) -> (anchor contig, overlap used)
    walks = {}   # (chrom, gap index) -> {"fwd": [...], "bwd": [...]}
    for chrom, placed in order.items():
        L = rlen[chrom]
        # gap i lies between placed[i-1] and placed[i]; gap 0 / len = chromosome ends
        for gi in range(len(placed) + 1):
            A = placed[gi - 1] if gi > 0 else None
            B = placed[gi] if gi < len(placed) else None
            lo = (aln[A[0]][1] - a.region_pad) if A else 0
            hi = (aln[B[0]][2] + a.region_pad) if B else L
            pool = [(n, st) for n, st, s, e in cands[chrom] if s >= lo and e <= hi]
            res = {"fwd": [], "bwd": []}
            fwd_closed = False
            for direction in ("fwd", "bwd"):
                if direction == "fwd":
                    start, target = A, B
                else:
                    # the backward walk aims for the end of the forward walk
                    start, target = B, (res["fwd"][-1] if res["fwd"] else A)
                if start is None or (direction == "bwd" and fwd_closed):
                    continue
                cur, used, why = start, set(), "max steps"
                for step in range(a.max_steps):
                    if target is not None:
                        closed = (overlaps(seqs, cur, target, a.anchor, a.min_id, a.min_overlap) if direction == "fwd"
                                  else overlaps(seqs, target, cur, a.anchor, a.min_id, a.min_overlap))
                        if closed:
                            why = "reached %s" % target[0]
                            fwd_closed = fwd_closed or direction == "fwd"
                            break
                    hits = []
                    for x in pool:
                        if x[0] in used or x in res["fwd"]:
                            continue
                        o = (overlaps(seqs, cur, x, a.anchor, a.min_id, a.min_overlap) if direction == "fwd"
                             else overlaps(seqs, x, cur, a.anchor, a.min_id, a.min_overlap))
                        if o:
                            hits.append((x, o))
                    if len(hits) != 1:
                        why = "no candidate" if not hits else "ambiguous: %s" % ",".join(h[0][0] for h in hits)
                        break
                    x, o = hits[0]
                    res[direction].append(x)
                    ovl[(chrom, gi, direction, x[0])] = (cur[0], min(o))
                    used.add(x[0])
                    log.append([chrom, gi, A[0] if A else "START", B[0] if B else "END", direction, step,
                                cur[0], x[0], x[1], ",".join(map(str, o))])
                    cur = x
                log.append([chrom, gi, A[0] if A else "START", B[0] if B else "END", direction, "stop",
                            cur[0], ".", ".", why])
            walks[(chrom, gi)] = res

    # Reject contigs claimed by more than one walk; truncate walks before them.
    claims = Counter(x[0] for r in walks.values() for d in ("fwd", "bwd") for x in r[d])
    conflicted = set(n for n, c in claims.items() if c > 1)
    for r in walks.values():
        for d in ("fwd", "bwd"):
            for i, x in enumerate(r[d]):
                if x[0] in conflicted:
                    r[d] = r[d][:i]
                    break

    added = 0
    with open(a.out, "w") as o:
        for chrom, placed in order.items():
            seq = []
            for gi in range(len(placed) + 1):
                r = walks[(chrom, gi)]
                seq.extend(r["fwd"])
                seq.extend(x for x in reversed(r["bwd"]) if x not in r["fwd"])
                added += len(r["fwd"]) + len([x for x in r["bwd"] if x not in r["fwd"]])
                if gi < len(placed):
                    seq.append(placed[gi])
            for n, st in seq:
                o.write("%s\t0\t0\t%s\t.\t%s\t0\t0\n" % (chrom, n, st))
    with open(a.extents, "w") as f:
        f.write("contig\tdirection\tanchor\toverlap\n")
        for (chrom, gi), r in walks.items():
            for d in ("fwd", "bwd"):
                for x in r[d]:
                    anchor, o = ovl[(chrom, gi, d, x[0])]
                    f.write("%s\t%s\t%s\t%d\n" % (x[0], d, anchor, o))
    with open(a.log, "w") as f:
        f.write("# overlap_extension_prototype.py %s %s\n" % (__version__, " ".join(sys.argv[1:])))
        f.write("chrom\tgap\tleft\tright\tdirection\tstep\tfrom\tadded\tstrand\toverlaps_or_reason\n")
        for r in log:
            f.write("\t".join(map(str, r)) + "\n")
    sys.stderr.write("added %d contigs; %d rejected as claimed by several walks\n" % (added, len(conflicted)))


if __name__ == "__main__":
    main()
