#!/usr/bin/env python3
"""Summarize GPatch runs at several MAPQ thresholds (-m) on the T2T-HG002 simulation.

Expects GPatch_v0.6.0/mapq_sweep/m<N>/ to hold, for prefix P, the GPatch outputs
(P.contigs.bed, P.unplaced.tsv, P.extension.tsv) and the evaluation outputs
P.trim-vs-overlap.summary.txt (compare_overlap_trimming.py), P.contigs.ordering-stats.txt
(compare_contig_ordering.py) and P.contigs.mapping-summary.txt (compile_mapping_stats.Rscript).

Writes one tab-separated row per threshold to stdout.
"""

__version__ = "0.1.0"

import argparse
import ast
import os
import re
import sys
from collections import Counter

P = "T2T-HG002.HG002.1.simulated.merge100000"


def summary_vector(path):
    """Parse a print()ed named R vector (alternating name / value lines)."""
    lines = [l.split() for l in open(path) if l.strip()]
    out = {}
    for names, vals in zip(lines[0::2], lines[1::2]):
        out.update(zip(names, map(float, vals)))
    return out


def junction_error(block):
    m = re.search(r"ALL units with truth.*?extra=(\d+) missing=(\d+)", block)
    return int(m.group(1)), int(m.group(2))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-d", "--sweep-dir", default="GPatch_v0.6.0/mapq_sweep")
    ap.add_argument("-s", "--sim-bed", default=P + ".contigs.bed")
    ap.add_argument("-m", "--mapq", default="0,1,10,20,30,40")
    ap.add_argument("-b", "--baseline", type=int, default=10,
                    help="Threshold that the placed-sequence difference is relative to.")
    args = ap.parse_args()

    truth = {}
    for l in open(args.sim_bed):
        x = l.split()
        truth[x[3]] = int(x[2]) - int(x[1])
    n_total = len(truth)

    rows = []
    for m in map(int, args.mapq.split(",")):
        d = os.path.join(args.sweep_dir, "m%d" % m, P)
        placed = [l.split()[3] for l in open(d + ".contigs.bed")]
        unplaced = Counter(l.split("\t")[2] for l in open(d + ".unplaced.tsv")
                           if not l.startswith(("#", "contig\t")))
        # Contigs inserted by an extension walk (numeric step), excluding those a later
        # "rejected" row withdrew, and counted only if they ended up placed.
        added, rejected = set(), set()
        for l in open(d + ".extension.tsv"):
            if l.startswith(("#", "chrom\t")):
                continue
            x = l.split("\t")
            if x[5] == "rejected":
                rejected.add(x[7])
            elif x[5] != "stop" and x[7] != ".":
                added.add(x[7])
        ext = len((added - rejected) & set(placed))
        s = open(d + ".trim-vs-overlap.summary.txt").read()
        units = ast.literal_eval(re.search(r"units: (\{.*?\})", s).group(1))
        misor = int(re.search(r"opposite orientation to the benchmark: (\d+)", s).group(1))
        all_blk, auto_blk = s.split("-- excluding chrY --")
        ae, am = junction_error(auto_blk)
        ge, gm = junction_error(all_blk)
        o = [l.rstrip("\n").split("\t") for l in open(d + ".contigs.ordering-stats.txt")]
        avg = dict(zip(o[0], next(r for r in o if r[0] == "Average")))
        ms = summary_vector(d + ".contigs.mapping-summary.txt")
        rows.append({
            "mapq": m,
            "placed": len(placed),
            "pct_contigs_placed": len(placed) / n_total,
            "placed_by_extension": ext,
            "placed_mb": sum(truth[n] for n in placed) / 1e6,
            "unplaced_low_mapq": unplaced.get("low_mapq", 0),
            "unplaced_nested": unplaced.get("nested", 0),
            "diff_chrom_joins": units.get("diff_chrom", 0),
            "not_adjacent_joins": units.get("not_adjacent", 0),
            "misoriented": misor,
            "auto_extra_mb": ae / 1e6,
            "auto_missing_mb": am / 1e6,
            "auto_junction_error_mb": (ae + am) / 1e6,
            "all_junction_error_mb": (ge + gm) / 1e6,
            "clustering_accuracy": ms["clustering_accuracy"],
            "orientation_accuracy": ms["orientation_accuracy"],
            "norm_edit_dist": float(avg["norm_edit_dist"]),
            "pct_pairs_recovered": float(avg["pct_pairs_recovered"]),
            "pct_pairs_correct": float(avg["pct_pairs_correct"]),
            "switch_err_rate": float(avg["switch_err_rate"]),
        })
    base = next((r["placed_mb"] for r in rows if r["mapq"] == args.baseline), None)
    cols = list(rows[0]) + ["placed_mb_vs_baseline"]
    print("\t".join(cols))
    for r in rows:
        r["placed_mb_vs_baseline"] = r["placed_mb"] - base if base is not None else float("nan")
        print("\t".join(("%.4f" % r[c]) if isinstance(r[c], float) else str(r[c]) for c in cols))


if __name__ == "__main__":
    sys.exit(main())
