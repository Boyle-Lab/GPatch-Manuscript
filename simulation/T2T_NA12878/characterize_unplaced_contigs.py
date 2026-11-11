#!/usr/bin/env python3
"""
Characterize simulated contigs that GPatch left unplaced (unplaced.tsv), as was
done for the T2T-HG002 simulation, for a benchmark split into several contigs.

For every simulated contig the script records GPatch's status (placed, or the
unplaced reason) and:
  - primary alignment to the patching reference (BAM): chromosome, MAPQ,
    aligned fraction, supplementary alignments (MAPQ >= 10 and >= 10 kb);
  - repeat context at the primary alignment, from reference (CHM13) CenSat
    (excluding ct transition and gap entries) and SD annotations, and from the
    contig sequence itself (21-mer diversity, HSat2/3 CATTC/GAATG content),
    because the benchmark has no repeat annotation;
  - location relative to its placed benchmark neighbours: the nearest placed
    contigs before and after it in the same benchmark contig, on the reference
    chromosome named after the benchmark chromosome. The expected region spans
    their primary alignments, +/- --bracket-pad; with one placed neighbour it
    is that neighbour's alignment +/- (benchmark distance + pad);
  - overlap (>= --min-overlap bp) with its true benchmark neighbours, from the
    simulation junctions;
  - for nested contigs, whether its primary alignment lies within the primary
    alignment of the containing contig ("nested by aligned coordinates") or
    only within the containing contig's chained/extrapolated extent.
Overlap extension results come from extension.tsv (inserted contigs and walk
stop reasons).

Outputs a per-contig TSV (-o) and a summary on stdout.
"""

import re
import bisect
import statistics as st
from argparse import ArgumentParser
from collections import defaultdict, Counter

import pysam

__version__ = "0.1.1"


def merged_intervals(bed_f, keep=None):
    iv = defaultdict(list)
    for l in open(bed_f):
        if l.startswith(("track", "#")) or not l.strip():
            continue
        x = l.rstrip("\n").split("\t")
        if keep and not keep(x):
            continue
        iv[x[0]].append((int(x[1]), int(x[2])))
    out = {}
    for c, v in iv.items():
        v.sort()
        m = []
        for s, e in v:
            if m and s <= m[-1][1]:
                m[-1][1] = max(m[-1][1], e)
            else:
                m.append([s, e])
        out[c] = m
    return out


def frac(m, c, s, e):
    if e <= s:
        return 0.0
    v = m.get(c, [])
    i = max(0, bisect.bisect_left(v, [s, 0]) - 1)
    tot = 0
    while i < len(v) and v[i][0] < e:
        tot += max(0, min(e, v[i][1]) - max(s, v[i][0]))
        i += 1
    return tot / (e - s)


def composition(seq, k=21):
    s = seq.upper()
    n = len(s)
    if n < k:
        return float("nan"), float("nan")
    div = len(set(s[i:i + k] for i in range(n - k + 1))) / (n - k + 1)
    hs = (s.count("CATTC") + s.count("GAATG")) * 5 / n
    return div, hs


def main():
    p = ArgumentParser(description="Characterize simulated contigs left unplaced by GPatch.")
    p.add_argument("-s", "--sim-bed", required=True, help="Simulated contig BED (benchmark coordinates).")
    p.add_argument("-j", "--junctions", required=True, help="Simulation junctions.tsv.")
    p.add_argument("-a", "--assembly", required=True, help="Simulated contig FASTA.")
    p.add_argument("-g", "--gpatch-bed", required=True, help="GPatch contigs.bed.")
    p.add_argument("-u", "--unplaced", required=True, help="GPatch unplaced.tsv.")
    p.add_argument("-e", "--extension", default=None, help="GPatch extension.tsv.")
    p.add_argument("-b", "--bam", required=True, help="BAM used for patching.")
    p.add_argument("--censat", required=True, help="Reference CenSat BED.")
    p.add_argument("--sd", required=True, help="Reference segmental-duplication BED.")
    p.add_argument("--sep", default=".", help="Benchmark contig ID separator (default '.').")
    p.add_argument("--bracket-pad", type=int, default=2000000)
    p.add_argument("--min-overlap", type=int, default=500)
    p.add_argument("-o", "--out", required=True, help="Per-contig TSV.")
    a = p.parse_args()

    truth, by_bench = {}, defaultdict(list)
    for l in open(a.sim_bed):
        x = l.rstrip("\n").split("\t")
        truth[x[3]] = dict(bench=x[0], chrom=x[0].split(a.sep)[0], s=int(x[1]), e=int(x[2]), strand=x[5])
        by_bench[x[0]].append(x[3])
    idx = {}
    for b in by_bench:
        by_bench[b].sort(key=lambda n: (truth[n]["s"], -truth[n]["e"]))
        for i, n in enumerate(by_bench[b]):
            idx[n] = i
    placed = {l.split("\t")[3] for l in open(a.gpatch_bed)}
    unpl = {}
    with open(a.unplaced) as f:
        next(f)
        for l in f:
            x = l.rstrip("\n").split("\t")
            unpl[x[0]] = dict(reason=x[2], detail=x[3])
    status = {n: ("placed" if n in placed else unpl.get(n, {}).get("reason", "missing")) for n in truth}

    prim, supp = {}, defaultdict(list)
    for r in pysam.AlignmentFile(a.bam):
        if r.is_unmapped:
            prim[r.query_name] = None
            continue
        if r.is_secondary:
            continue
        d = dict(chrom=r.reference_name, rs=r.reference_start, re=r.reference_end, mapq=r.mapping_quality,
                 alen=r.query_alignment_length, qlen=r.infer_read_length(), strand="-" if r.is_reverse else "+")
        if r.is_supplementary:
            supp[r.query_name].append(d)
        else:
            prim[r.query_name] = d

    cen = merged_intervals(a.censat, keep=lambda x: not re.match(r"(ct|gap)", x[3], re.I))
    sd = merged_intervals(a.sd)

    ovl = defaultdict(list)   # contig -> overlaps with true neighbours (bp)
    for l in open(a.junctions):
        if l.startswith("chrom"):
            continue
        x = l.rstrip("\n").split("\t")
        d = int(x[5])
        if d <= -a.min_overlap:
            ovl[x[1]].append(-d)
            ovl[x[2]].append(-d)

    inserted, stops = set(), Counter()
    if a.extension:
        for l in open(a.extension):
            if l.startswith(("#", "chrom\t")):
                continue
            x = l.rstrip("\n").split("\t")
            if x[5] == "stop":
                stops[re.sub(r":.*", "", re.sub(r"reached .*", "reached flanking contig", x[9]))] += 1
            elif x[5] != "rejected" and x[7] != ".":
                inserted.add(x[7])

    fa = pysam.FastaFile(a.assembly)

    def neighbour(n, step):
        lst = by_bench[truth[n]["bench"]]
        i = idx[n] + step
        while 0 <= i < len(lst):
            m = lst[i]
            if status[m] == "placed" and prim.get(m) and prim[m]["chrom"] == truth[n]["chrom"]:
                return m
            i += step
        return None

    def location(n):
        pr = prim.get(n)
        if not pr:
            return "unaligned"
        if pr["chrom"] != truth[n]["chrom"]:
            return "wrong_chrom"
        L, R = neighbour(n, -1), neighbour(n, 1)
        if not L and not R:
            return "right_chrom_no_placed_neighbour"
        if L and R:
            lo = min(prim[L]["rs"], prim[R]["rs"]) - a.bracket_pad
            hi = max(prim[L]["re"], prim[R]["re"]) + a.bracket_pad
        else:
            m = L or R
            dist = abs(truth[m]["s"] - truth[n]["s"]) + (truth[n]["e"] - truth[n]["s"])
            lo, hi = prim[m]["rs"] - dist - a.bracket_pad, prim[m]["re"] + dist + a.bracket_pad
        return "right_chrom_in_bracket" if lo <= pr["rs"] and pr["re"] <= hi else "right_chrom_elsewhere"

    rows = []
    for n, t in truth.items():
        pr = prim.get(n)
        L = t["e"] - t["s"]
        cf = frac(cen, pr["chrom"], pr["rs"], pr["re"]) if pr else float("nan")
        sf = frac(sd, pr["chrom"], pr["rs"], pr["re"]) if pr else float("nan")
        # Sequence composition only for unplaced contigs (placed ones include
        # chromosome-scale contigs and are not summarized by composition).
        div, hs = composition(fa.fetch(n)) if status[n] != "placed" else (float("nan"), float("nan"))
        nest_aln = "NA"
        if status[n] == "nested":
            m = re.search(r"of (\S+)", unpl[n]["detail"])
            c = prim.get(m.group(1)) if m else None
            nest_aln = bool(c and pr and c["chrom"] == pr["chrom"] and c["rs"] <= pr["rs"] and pr["re"] <= c["re"])
        bench_any_placed = any(status[m] == "placed" for m in by_bench[t["bench"]])
        rows.append(dict(
            contig=n, bench=t["bench"], length=L, status=status[n],
            inserted_by_extension=(n in inserted and status[n] == "placed"),
            mapq=pr["mapq"] if pr else -1, ref_chrom=pr["chrom"] if pr else ".",
            ref_start=pr["rs"] if pr else -1, ref_end=pr["re"] if pr else -1,
            aligned_frac=round(pr["alen"] / pr["qlen"], 4) if pr else 0,
            n_supp=len(supp.get(n, [])),
            hq_supp=sum(1 for s in supp.get(n, []) if s["mapq"] >= 10 and s["alen"] >= 10000),
            censat_frac=round(cf, 4), sd_frac=round(sf, 4), kmer21_div=round(div, 4), hsat23_frac=round(hs, 4),
            satellite_like=(div < 0.5 or hs > 0.2), location=location(n),
            max_true_overlap=max(ovl.get(n, [0])), nested_by_aligned=nest_aln,
            bench_contig_has_placed=bench_any_placed))

    cols = list(rows[0])
    with open(a.out, "w") as o:
        o.write("# characterize_unplaced_contigs.py %s\n" % __version__)
        o.write("\t".join(cols) + "\n")
        for r in rows:
            o.write("\t".join(str(r[c]) for c in cols) + "\n")

    # ---- summary ---------------------------------------------------------
    w = print
    tot_bp = sum(r["length"] for r in rows)
    un = [r for r in rows if r["status"] != "placed"]
    w("simulated contigs: %d (%.2f Mb); placed %d (%.2f Mb, %d by overlap extension); unplaced %d (%.2f Mb, %.2f%%)"
      % (len(rows), tot_bp / 1e6, len(rows) - len(un), (tot_bp - sum(r["length"] for r in un)) / 1e6,
         sum(r["inserted_by_extension"] for r in rows), len(un), sum(r["length"] for r in un) / 1e6,
         100 * sum(r["length"] for r in un) / tot_bp))

    def lb(l):
        return "<10kb" if l < 1e4 else "10-100kb" if l < 1e5 else "100kb-1Mb" if l < 1e6 else ">=1Mb"

    def mb(q):
        return "unaligned" if q < 0 else "0" if q == 0 else "1-9" if q < 10 else "10-29" if q < 30 else "30-59" if q < 60 else "60"

    for stt in ("placed", "low_mapq", "nested", "unmapped"):
        R = [r for r in rows if r["status"] == stt]
        if not R:
            continue
        Ls = [r["length"] for r in R]
        w("\n== %s: n=%d total %.2f Mb; median length %.1f kb (min %d bp, max %.2f Mb)"
          % (stt, len(R), sum(Ls) / 1e6, st.median(Ls) / 1e3, min(Ls), max(Ls) / 1e6))
        w("   length bins: %s" % dict(Counter(lb(l) for l in Ls)))
        w("   MAPQ: %s" % dict(sorted(Counter(mb(r["mapq"]) for r in R).items())))
        al = [r for r in R if r["mapq"] >= 0]
        if al:
            cs = sum(r["censat_frac"] >= .5 for r in al)
            sdn = sum(r["sd_frac"] >= .5 for r in al)
            either = [r for r in al if r["censat_frac"] >= .5 or r["sd_frac"] >= .5]
            w("   reference repeat context at alignment (>=50%%): CenSat %d, SD %d, either %d (%.2f Mb), neither %d"
              % (cs, sdn, len(either), sum(r["length"] for r in either) / 1e6, len(al) - len(either)))
        sat = [r for r in R if r["satellite_like"]]
        if stt != "placed":
            w("   satellite-like sequence (21-mer diversity < 0.5 or HSat2/3 > 20%%): %d (%.2f Mb)"
            % (len(sat), sum(r["length"] for r in sat) / 1e6))
            w("   CenSat/SD at alignment OR satellite-like: %d (%.2f Mb)"
            % (sum(1 for r in R if r["satellite_like"] or r["censat_frac"] >= .5 or r["sd_frac"] >= .5),
               sum(r["length"] for r in R if r["satellite_like"] or r["censat_frac"] >= .5 or r["sd_frac"] >= .5) / 1e6))
        if al:
            w("   primary aligned fraction median %.3f; with >=1 supplementary (MAPQ>=10, >=10 kb): %d"
              % (st.median(r["aligned_frac"] for r in al), sum(r["hq_supp"] > 0 for r in al)))
        w("   location vs placed benchmark neighbours: %s" % dict(Counter(r["location"] for r in R)))
        w("   by reference chromosome (top 10): %s" % dict(Counter(r["ref_chrom"] for r in R).most_common(10)))
        if stt != "placed":
            ov = [r for r in R if r["max_true_overlap"] >= 5000]
            w("   overlap a true neighbour by >=%d bp: %d; by >=5 kb: %d (%.2f Mb)"
              % (a.min_overlap, sum(r["max_true_overlap"] > 0 for r in R), len(ov), sum(r["length"] for r in ov) / 1e6))
            nb = [r for r in R if not r["bench_contig_has_placed"]]
            w("   in benchmark contigs with no placed contig: %d (%.2f Mb)" % (len(nb), sum(r["length"] for r in nb) / 1e6))
        if stt == "nested":
            w("   nested by aligned coordinates: %d; only by chained/extrapolated extent: %d"
              % (sum(r["nested_by_aligned"] is True for r in R), sum(r["nested_by_aligned"] is False for r in R)))
    if stops:
        w("\nextension walk stops: %s" % dict(stops))


if __name__ == "__main__":
    main()
