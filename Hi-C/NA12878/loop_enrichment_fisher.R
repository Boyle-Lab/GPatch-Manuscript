#!/usr/bin/env Rscript
# Usage: loop_enrichment_fisher.R <output.tsv>
# One-sided Fisher exact tests for enrichment of loops >=1Mb and >=2Mb among
# scaffold-only loops absent from the patched assembly. Counts are hard-coded
# from the output of loop_length_dists.py (exact anchor-bin matching).
# test group: scaffold-only loops absent from patched
d <- data.frame(
  res   = c("5kb","5kb","5kb","10kb","10kb","10kb"),
  bg    = rep(c("scaffold-only loops also in patched","all patched loops","all T2T-CHM13 loops"),2),
  t_n   = c(1277,1277,1277,4264,4264,4264),
  t_1M  = c(366,366,366,2480,2480,2480),
  t_2M  = c(281,281,281,2050,2050,2050),
  b_n   = c(9591,9786,9987,19555,20056,19918),
  b_1M  = c(15,16,18,222,256,187),
  b_2M  = c(0,0,0,34,54,22))
out <- list()
for (i in seq_len(nrow(d))) for (th in c("1M","2M")) {
  tk <- d[[paste0("t_",th)]][i]; bk <- d[[paste0("b_",th)]][i]
  m <- matrix(c(tk, d$t_n[i]-tk, bk, d$b_n[i]-bk), 2)
  f <- fisher.test(m, alternative="greater")
  out[[length(out)+1]] <- data.frame(res=d$res[i], threshold=paste0(">=",th,"b"), background=d$bg[i],
    test=sprintf("%d/%d (%.1f%%)",tk,d$t_n[i],100*tk/d$t_n[i]),
    bg=sprintf("%d/%d (%.2f%%)",bk,d$b_n[i],100*bk/d$b_n[i]),
    odds_ratio=signif(f$estimate,4), OR_95CI_low=signif(f$conf.int[1],4), p=signif(f$p.value,3))
}
r <- do.call(rbind,out); rownames(r) <- NULL
r$p_bonf <- signif(pmin(1, r$p*nrow(r)),3)
print(r, right=FALSE)
write.table(r, commandArgs(TRUE)[1], sep="\t", quote=FALSE, row.names=FALSE)
