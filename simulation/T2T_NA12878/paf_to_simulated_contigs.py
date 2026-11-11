#!/usr/bin/env python3
"""
Build a simulated contig-level assembly with known order and orientation from
minimap2 alignments of a model contig-level assembly to a chromosome-scale
benchmark genome.

Each primary alignment block (PAF tp:A:P) that passes filters defines one
simulated contig spanning [tstart, tend) on the benchmark chromosome. Contig
sequence is extracted from the benchmark genome and reverse-complemented when
the block aligned on the - strand. Contigs are named <chrom>_<rank>, with rank
the 0-based order of the contig along the chromosome (sorted by start, then
longest first). Overlaps and unfilled gaps between adjacent contigs are
inherited from the alignment boundaries and are not altered.

Outputs (for --out-prefix P):
  P.contigs.fa         simulated contig sequences
  P.contigs.bed        8-column BED (no header):
                         1 chrom  2 start  3 end  4 name  5 length  6 strand
                         7 mapq   8 semicolon-delimited attributes:
                           src (model contig), qlen, qs, qe (query coords),
                           AS, NM, matches, blen, de (gap-compressed divergence),
                           nblk (alignment blocks merged into this contig)
                       Column 5 holds length (not a 0-1000 score) and the file
                       is capped at 8 columns for compatibility with
                       compare_contig_ordering.py / compare_contig_positions.py.
  P.junctions.tsv      one row per adjacent contig pair: signed distance
                       (negative = overlap, 0 = abutting, positive = gap)
  P.summary.txt        summary statistics
"""

import sys
import statistics
from argparse import ArgumentParser
from collections import defaultdict, Counter

from Bio import SeqIO

__version__ = "0.2.0"

PAF_COLS = ["qname", "qlen", "qstart", "qend", "strand", "tname", "tlen",
            "tstart", "tend", "matches", "blen", "mapq"]
INT_COLS = {"qlen", "qstart", "qend", "tlen", "tstart", "tend", "matches",
            "blen", "mapq"}
FASTA_WIDTH = 80


def parse_paf(paf_f):
    with open(paf_f) as paf:
        for line in paf:
            f = line.rstrip("\n").split("\t")
            rec = {k: (int(v) if k in INT_COLS else v)
                   for k, v in zip(PAF_COLS, f[:12])}
            rec["tags"] = {}
            for t in f[12:]:
                key, typ, val = t.split(":", 2)
                rec["tags"][key] = val
            yield rec


def n50(lengths):
    tot = sum(lengths)
    acc = 0
    for l in sorted(lengths, reverse=True):
        acc += l
        if acc >= tot / 2:
            return l
    return 0


def describe(vals):
    if not vals:
        return "n=0"
    return "n=%d min=%d median=%d mean=%.1f max=%d" % (
        len(vals), min(vals), statistics.median(vals),
        statistics.mean(vals), max(vals))


def merged_coverage(intervals):
    """Number of bases covered by a list of (start, end) intervals."""
    cov, cur_s, cur_e = 0, None, None
    for s, e in sorted(intervals):
        if cur_e is None or s > cur_e:
            if cur_e is not None:
                cov += cur_e - cur_s
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        cov += cur_e - cur_s
    return cov


def merge_same_source(alns, max_dist):
    """
    Merge consecutive alignment blocks from the same model contig and strand
    that are colinear on both query and target and separated by at most
    max_dist bp (gap or overlap) on both. Merged blocks span the union of
    their parts; matches/blen/NM/AS are summed, mapq is the minimum and de is
    the blen-weighted mean. Returns (merged list, number of merge events).
    """
    groups = defaultdict(list)
    for a in alns:
        a.setdefault("nblk", 1)
        groups[(a["qname"], a["strand"])].append(a)
    out, n_merges = [], 0
    for (q, strand), g in groups.items():
        g.sort(key=lambda a: (a["tstart"], -a["tend"]))
        cur = g[0]
        for nxt in g[1:]:
            tdist = nxt["tstart"] - cur["tend"]
            if strand == "+":
                qdist = nxt["qstart"] - cur["qend"]
                colinear = nxt["qend"] > cur["qend"]
            else:
                qdist = cur["qstart"] - nxt["qend"]
                colinear = nxt["qstart"] < cur["qstart"]
            if (colinear and nxt["tend"] > cur["tend"]
                    and abs(tdist) <= max_dist and abs(qdist) <= max_dist):
                cur = combine(cur, nxt)
                n_merges += 1
            else:
                out.append(cur)
                cur = nxt
        out.append(cur)
    return out, n_merges


def combine(a, b):
    m = dict(a)
    for k in ("tstart", "qstart"):
        m[k] = min(a[k], b[k])
    for k in ("tend", "qend"):
        m[k] = max(a[k], b[k])
    m["matches"] = a["matches"] + b["matches"]
    m["blen"] = a["blen"] + b["blen"]
    m["mapq"] = min(a["mapq"], b["mapq"])
    m["nblk"] = a["nblk"] + b["nblk"]
    ta, tb, t = a["tags"], b["tags"], dict(a["tags"])
    for k in ("AS", "NM"):
        t[k] = str(int(ta[k]) + int(tb[k])) if k in ta and k in tb else "NA"
    if "de" in ta and "de" in tb:
        t["de"] = "%.4f" % ((float(ta["de"]) * a["blen"] + float(tb["de"]) * b["blen"])
                            / (a["blen"] + b["blen"]))
    else:
        t["de"] = "NA"
    m["tags"] = t
    return m


def write_fasta(out, name, seq):
    out.write(">%s\n" % name)
    for i in range(0, len(seq), FASTA_WIDTH):
        out.write(seq[i:i + FASTA_WIDTH])
        out.write("\n")


def main():
    p = ArgumentParser(description=__doc__.split("\n\n")[1])
    p.add_argument("-p", "--paf", required=True,
                   help="minimap2 PAF of model contigs (query) vs benchmark (target); run with -c.")
    p.add_argument("-r", "--benchmark", required=True,
                   help="Benchmark genome FASTA (the PAF target).")
    p.add_argument("-o", "--out-prefix", required=True)
    p.add_argument("--min-mapq", type=int, default=0,
                   help="Minimum MAPQ for an alignment block (default 0).")
    p.add_argument("--min-len", type=int, default=1000,
                   help="Minimum benchmark span of an alignment block (default 1000).")
    p.add_argument("--merge-dist", type=int, default=None, metavar="N",
                   help="Merge consecutive colinear same-strand blocks of the same model contig "
                        "separated by <= N bp (gap or overlap) on both benchmark and model "
                        "coordinates. Applied before --min-len. Default: no merging.")
    p.add_argument("--keep-contained", action="store_true",
                   help="Keep blocks whose benchmark interval lies entirely within another kept block (default: drop).")
    p.add_argument("--version", action="version", version=__version__)
    args = p.parse_args()

    # ---- Collect and filter alignment blocks --------------------------------
    filt = Counter()
    query_lens = {}
    blocks = defaultdict(list)
    chrom_lens = {}
    for a in parse_paf(args.paf):
        query_lens[a["qname"]] = a["qlen"]
        chrom_lens[a["tname"]] = a["tlen"]
        filt["total"] += 1
        if a["tags"].get("tp") != "P":
            filt["not_primary"] += 1
            continue
        if a["mapq"] < args.min_mapq:
            filt["low_mapq"] += 1
            continue
        blocks[a["tname"]].append(a)

    kept = {}
    for chrom, alns in blocks.items():
        if args.merge_dist is not None:
            alns, n = merge_same_source(alns, args.merge_dist)
            filt["merged"] += n
        filt["short"] += sum(1 for a in alns if a["tend"] - a["tstart"] < args.min_len)
        alns = [a for a in alns if a["tend"] - a["tstart"] >= args.min_len]
        alns.sort(key=lambda a: (a["tstart"], -a["tend"]))
        out, max_end = [], -1
        for a in alns:
            if not args.keep_contained and a["tend"] <= max_end:
                filt["contained"] += 1
                continue
            out.append(a)
            max_end = max(max_end, a["tend"])
        for i, a in enumerate(out):
            a["name"] = "%s_%d" % (chrom, i)
        kept[chrom] = out

    # ---- Extract sequences from benchmark ----------------------------------
    fa_f = args.out_prefix + ".contigs.fa"
    written = 0
    chrom_order = []
    with open(fa_f, "w") as fa:
        for rec in SeqIO.parse(args.benchmark, "fasta"):
            if rec.id not in kept:
                continue
            chrom_order.append(rec.id)
            if len(rec.seq) != chrom_lens[rec.id]:
                sys.exit("ERROR: length of %s in benchmark FASTA (%d) != PAF tlen (%d)"
                         % (rec.id, len(rec.seq), chrom_lens[rec.id]))
            for a in kept[rec.id]:
                s = rec.seq[a["tstart"]:a["tend"]]
                if a["strand"] == "-":
                    s = s.reverse_complement()
                write_fasta(fa, a["name"], str(s))
                written += 1
            sys.stderr.write("extracted %s (%d contigs)\n" % (rec.id, len(kept[rec.id])))
    if written != sum(len(v) for v in kept.values()):
        sys.exit("ERROR: wrote %d sequences but %d contigs defined; benchmark FASTA is missing chromosomes"
                 % (written, sum(len(v) for v in kept.values())))

    # ---- Write BED and junctions in benchmark chromosome order ------------------------
    bed_f = args.out_prefix + ".contigs.bed"
    junc_f = args.out_prefix + ".junctions.tsv"
    junc_dists = []
    end_gaps = []
    src_blocks = Counter()
    lengths = []
    strands = Counter()
    covered = 0
    with open(bed_f, "w") as bed, open(junc_f, "w") as junc:
        junc.write("chrom\tleft_contig\tright_contig\tleft_end\tright_start\tdistance\ttype\tsame_source\n")
        for chrom in chrom_order:
            alns = kept[chrom]
            for a in alns:
                L = a["tend"] - a["tstart"]
                lengths.append(L)
                strands[a["strand"]] += 1
                src_blocks[a["qname"]] += 1
                t = a["tags"]
                attrs = "src=%s;qlen=%d;qs=%d;qe=%d;AS=%s;NM=%s;matches=%d;blen=%d;de=%s;nblk=%d" % (
                    a["qname"], a["qlen"], a["qstart"], a["qend"],
                    t.get("AS", "NA"), t.get("NM", "NA"), a["matches"],
                    a["blen"], t.get("de", "NA"), a.get("nblk", 1))
                bed.write("%s\t%d\t%d\t%s\t%d\t%s\t%d\t%s\n" % (
                    chrom, a["tstart"], a["tend"], a["name"], L,
                    a["strand"], a["mapq"], attrs))
            for left, right in zip(alns, alns[1:]):
                d = right["tstart"] - left["tend"]
                typ = "overlap" if d < 0 else ("abut" if d == 0 else "gap")
                junc_dists.append(d)
                junc.write("%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\n" % (
                    chrom, left["name"], right["name"], left["tend"],
                    right["tstart"], d, typ, left["qname"] == right["qname"]))
            end_gaps.append(alns[0]["tstart"])
            end_gaps.append(chrom_lens[chrom] - alns[-1]["tend"])
            covered += merged_coverage([(a["tstart"], a["tend"]) for a in alns])

    # ---- Summary ------------------------------------------------------------
    genome_len = sum(chrom_lens.values())
    jt = Counter("overlap" if d < 0 else ("abut" if d == 0 else "gap") for d in junc_dists)
    nj = max(len(junc_dists), 1)
    split_src = sum(1 for q, n in src_blocks.items() if n > 1)
    with open(args.out_prefix + ".summary.txt", "w") as s:
        w = lambda *x: s.write("\t".join(str(v) for v in x) + "\n")
        w("# paf_to_simulated_contigs.py", __version__)
        w("# options", "min_mapq=%d min_len=%d merge_dist=%s keep_contained=%s" % (args.min_mapq, args.min_len, args.merge_dist, args.keep_contained))
        w("paf_records", filt["total"])
        w("filtered_not_primary", filt["not_primary"])
        w("filtered_low_mapq", filt["low_mapq"])
        w("blocks_merged_away", filt["merged"])
        w("filtered_short", filt["short"])
        w("filtered_contained", filt["contained"])
        w("simulated_contigs", len(lengths))
        w("strand_plus", strands["+"])
        w("strand_minus", strands["-"])
        w("contig_length", describe(lengths))
        w("contig_N50", n50(lengths))
        w("model_contigs_in_paf", len(query_lens))
        w("model_contigs_used", len(src_blocks))
        w("model_contigs_split_into_multiple", split_src)
        w("model_assembly_N50", n50(list(query_lens.values())))
        w("benchmark_chroms_with_contigs", "%d/%d" % (len(kept), len(chrom_lens)))
        w("benchmark_bases_covered", "%d/%d (%.3f%%)" % (covered, genome_len, 100 * covered / genome_len))
        w("junctions", len(junc_dists))
        w("junctions_overlap", "%d (%.2f%%)" % (jt["overlap"], 100 * jt["overlap"] / nj))
        w("junctions_abut", "%d (%.2f%%)" % (jt["abut"], 100 * jt["abut"] / nj))
        w("junctions_gap", "%d (%.2f%%)" % (jt["gap"], 100 * jt["gap"] / nj))
        w("overlap_size", describe([-d for d in junc_dists if d < 0]))
        w("gap_size", describe([d for d in junc_dists if d > 0]))
        w("chrom_end_uncovered", describe(end_gaps))
    sys.stderr.write("Done: %d simulated contigs\n" % len(lengths))


if __name__ == "__main__":
    main()
