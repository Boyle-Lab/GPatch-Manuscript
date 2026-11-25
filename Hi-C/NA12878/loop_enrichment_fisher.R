#!/usr/bin/env Rscript
# Usage: loop_enrichment_fisher.R <prefix>
# One-sided Fisher's exact tests for enrichment of loops >=1Mb and >=2Mb among
# scaffold-only loops absent from the patched assembly (exact anchor-bin
# matching), against three backgrounds. Reads the tables written by
# loop_length_dists.py with the same <prefix>:
#   <prefix>.by_T2T_sharing.tsv, <prefix>.scaffold_vs_patched.tsv
# and writes <prefix>.fisher_long_loops.tsv.

prefix <- commandArgs(TRUE)[1]
t1 <- read.delim(paste0(prefix, ".by_T2T_sharing.tsv"), stringsAsFactors = FALSE)
t2 <- read.delim(paste0(prefix, ".scaffold_vs_patched.tsv"), stringsAsFactors = FALSE)

get2 <- function(res, group) {
  r <- t2[t2$resolution == res & t2$query_assembly == "scaffold-only" &
          t2$match_tolerance_bins == 0 & t2$group == group &
          t2$T2T_sharing_category == "all", ]
  stopifnot(nrow(r) == 1); r
}
get1 <- function(res, asm) {
  r <- t1[t1$resolution == res & t1$assembly == asm & t1$category == "all", ]
  stopifnot(nrow(r) == 1); r
}

out <- list()
for (res in c("5kb", "10kb")) {
  tst <- get2(res, "absent_from_reference")
  bgs <- list("scaffold-only loops also in patched" = get2(res, "present_in_reference"),
              "all patched loops" = get1(res, "patched"),
              "all T2T-CHM13 loops" = get1(res, "T2T-CHM13"))
  for (bg in names(bgs)) for (th in c("1Mb", "2Mb")) {
    col <- paste0("n_ge", th)
    tk <- tst[[col]]; tn <- tst$n
    bk <- bgs[[bg]][[col]]; bn <- bgs[[bg]]$n
    f <- fisher.test(matrix(c(tk, tn - tk, bk, bn - bk), 2), alternative = "greater")
    out[[length(out) + 1]] <- data.frame(
      resolution = res, threshold = paste0(">=", th),
      test_group = "scaffold-only loops absent from patched", background = bg,
      test_n_ge = tk, test_n = tn, test_frac = round(tk / tn, 4),
      bg_n_ge = bk, bg_n = bn, bg_frac = round(bk / bn, 4),
      odds_ratio = signif(unname(f$estimate), 4),
      OR_95CI_lower = signif(f$conf.int[1], 4),
      p_value = signif(f$p.value, 3))
  }
}
r <- do.call(rbind, out)
r$p_bonferroni <- signif(pmin(1, r$p_value * nrow(r)), 3)
# p-values below double precision are reported by R as exactly 0
r$p_note <- ifelse(r$p_value == 0, "p < 2.2e-308 (below machine precision)", "")
outfile <- paste0(prefix, ".fisher_long_loops.tsv")
write.table(r, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote ", outfile)
