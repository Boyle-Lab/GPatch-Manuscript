#!/usr/bin/env python3
"""Build formatted summary tables (as presented in the manuscript discussion)
from the full tables written by loop_length_dists.py and loop_enrichment_fisher.R.

Usage: make_display_tables.py <prefix>
Reads  <prefix>.by_T2T_sharing.tsv, <prefix>.scaffold_vs_patched.tsv,
       <prefix>.fisher_long_loops.tsv
Writes <prefix>.summary.by_T2T_sharing.tsv
       <prefix>.summary.scaffold_not_in_patched.tsv
       <prefix>.summary.fisher_long_loops.tsv"""
import csv, sys

P = sys.argv[1]
read = lambda s: list(csv.DictReader(open(f"{P}.{s}.tsv"), delimiter="\t"))
ASM = {"patched": "Patched", "scaffold-only": "Scaffold-only"}
CAT = {"fully_shared": "Fully shared", "partially_shared": "Partially shared", "unique": "Unique"}


def n(x):
    return f"{int(x):,}"


def length(x, median):
    x = float(x)
    if x < 1e6:
        return f"{round(x / 1000)} kb"
    return f"{x / 1e6:.3g} Mb" if median else f"{x / 1e6:.1f} Mb"


def pct(frac, zero="0"):
    v = 100 * float(frac)
    if v == 0:
        return zero
    return f"{v:.2f}%" if v < 10 else f"{v:.1f}%"


def pct1(frac):
    return f"{100 * float(frac):.1f}%"


def write(name, header, rows):
    path = f"{P}.summary.{name}.tsv"
    with open(path, "w") as fh:
        for r in [header] + rows:
            fh.write("\t".join(r) + "\n")
    sys.stderr.write(f"Wrote {path}\n")


# Table 1: length by T2T-CHM13 sharing category
rows = []
for r in read("by_T2T_sharing"):
    if r["assembly"] in ASM and r["category"] in CAT:
        rows.append([ASM[r["assembly"]], r["resolution"], CAT[r["category"]], n(r["n"]),
                     length(r["median"], True), length(r["mean"], False),
                     pct(r["frac_ge1Mb"]), pct(r["frac_ge2Mb"]),
                     pct1(r["frac_span_contig_junction"])])
write("by_T2T_sharing",
      ["Assembly", "Res", "Category", "n", "Median", "Mean", "≥1 Mb", "≥2 Mb",
       "Crosses a contig junction"], rows)

# Table 2: scaffold-only loops not found in the patched assembly
GRP = {"absent_from_reference": "Scaffold-only, not in patched",
       "present_in_reference": "In both assemblies"}
rows = []
for r in read("scaffold_vs_patched"):
    if (r["query_assembly"] == "scaffold-only" and r["match_tolerance_bins"] == "0"
            and r["T2T_sharing_category"] == "all"):
        rows.append([r["resolution"], GRP[r["group"]], n(r["n"]),
                     length(r["median"], True), length(r["mean"], False),
                     pct(r["frac_ge1Mb"]), pct(r["frac_ge2Mb"]),
                     pct1(r["frac_span_contig_junction"])])
write("scaffold_not_in_patched",
      ["Res", "Group", "n", "Median", "Mean", "≥1 Mb", "≥2 Mb",
       "Crosses a contig junction"], rows)

# Table 3: Fisher tests vs scaffold-only loops also found in patched
rows = []
for r in read("fisher_long_loops"):
    if r["background"] != "scaffold-only loops also in patched":
        continue
    orr = "∞" if r["odds_ratio"] == "Inf" else f"{round(float(r['odds_ratio']))}"
    p = "< 1e-300" if float(r["p_value"]) == 0 else f"{float(r['p_value']):.2g}"
    rows.append([r["resolution"], r["threshold"].replace(">=", "≥").replace("Mb", " Mb"),
                 f"{n(r['test_n_ge'])}/{n(r['test_n'])} ({pct(r['test_frac'], '0%')})",
                 f"{n(r['bg_n_ge'])}/{n(r['bg_n'])} ({pct(r['bg_frac'], '0%')})",
                 f"{orr} (≥{round(float(r['OR_95CI_lower']))})", p])
write("fisher_long_loops",
      ["Res", "Threshold", "Not in patched", "Also in patched",
       "Odds ratio (95% CI lower bound)", "p"], rows)
