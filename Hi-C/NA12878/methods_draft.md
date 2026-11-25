### Hi-C loop sets and loop length

We compared chromatin loops called from the same NA12878 (GM12878) in situ Hi-C data (SRA runs SRR1658572 and SRR1658592) mapped to three references: T2T-CHM13 v2.0, the GPatch-patched HGSVC NA12878 assembly (NA12878.2), and a scaffold-only version of the same GPatch assembly. The scaffold-only assembly has the same contig order, orientation and coordinates as the patched assembly, but the gaps between adjacent contigs are filled with N characters instead of T2T-CHM13 sequence. Loops were called with HiCCUPS (Juicer Tools) at 5 kb and 10 kb resolution. For each loop, we defined loop length as the distance from the start of the upstream anchor bin (x1) to the end of the downstream anchor bin (y2), in the coordinates of the assembly on which the loop was called. Lengths were therefore measured without liftover, and they could be compared directly between the patched and scaffold-only assemblies because those two assemblies share a coordinate system. We report the fraction of loops at least 1 Mb and at least 2 Mb long; for reference, these fractions were 0.18% and 0% (5 kb) and 0.94% and 0.11% (10 kb) for loops called against T2T-CHM13.

To relate long loops to assembly gaps, we recorded whether each loop spans a contig junction in the GPatch assemblies. Junction positions were taken as the first and last base of every contig placed in the NA12878.2 GPatch assembly, and a loop was counted as spanning a junction if any of these positions fell within its x1–y2 interval.

### Classifying loops by sharing with T2T-CHM13

For the patched and scaffold-only assemblies, each loop was split into its two anchors, and anchors were lifted to T2T-CHM13 coordinates with UCSC liftOver using chain files relating the GPatch NA12878 assembly to T2T-CHM13. A lifted anchor was counted as shared if it overlapped (by at least 1 bp; bedtools intersect) any T2T-CHM13 loop anchor called at the same resolution. Anchors were then reassembled into loops by loop identifier. A loop was classified as fully shared if both of its anchors were shared, partially shared if exactly one anchor was shared (the other either failed to lift over or did not overlap a T2T-CHM13 anchor), and unique otherwise. Because this definition requires each anchor to overlap some T2T-CHM13 anchor but not necessarily anchors of the same T2T-CHM13 loop, we also report the subset of fully shared loops whose two anchors overlap the two anchors of a single T2T-CHM13 loop. For each assembly, resolution and category we summarized the distribution of loop lengths (median, mean, quartiles and range), the fractions of loops at least 1 Mb and at least 2 Mb long, and the fraction spanning a contig junction.

### Loops called on the scaffold-only assembly but not the patched assembly

Because the patched and scaffold-only assemblies share coordinates, loops could be matched between them without liftover. A scaffold-only loop was considered present in the patched assembly if a patched loop at the same resolution was on the same chromosome and had identical upstream and downstream anchor bins. All other scaffold-only loops were considered absent from the patched assembly. To check that the result did not depend on exact matching, we repeated the comparison allowing each anchor to differ by up to one or two bins. We also ran the reverse comparison (patched loops absent from the scaffold-only assembly). For each set we summarized loop lengths, the fractions of loops at least 1 Mb and at least 2 Mb long, and the fraction spanning a contig junction, overall and by T2T-CHM13 sharing category.

### Statistical test for enrichment of long loops

We tested whether scaffold-only loops absent from the patched assembly were enriched for long loops using one-sided Fisher's exact tests (alternative: odds ratio greater than 1), separately at 5 kb and 10 kb resolution and at thresholds of 1 Mb and 2 Mb. Each test compared the proportion of loops at or above the threshold among scaffold-only loops absent from the patched assembly with the proportion among scaffold-only loops present in the patched assembly (exact anchor-bin matching). As alternative backgrounds, we also compared against all loops called on the patched assembly and all loops called on T2T-CHM13. We report odds ratios with the lower bound of the 95% confidence interval and Bonferroni-corrected p-values across all 12 tests. P-values smaller than the smallest representable double-precision number (2.2 × 10⁻³⁰⁸) are reported as below 10⁻³⁰⁰. Tests were run in R 4.3.0 (fisher.test).

<!--
Provenance / notes for revision:
- Code (Hi-C/NA12878): loop_length_dists.py (Python 3.12, bedtools v2.29.1),
  loop_enrichment_fisher.R (R 4.3.0), make_display_tables.py, reconstitute_loops.py.
  Anchor decomposition, liftOver and intersection commands are in Hi-C/NA12878/README.md.
- Results: loop_lengths/NA12878_loops.{by_T2T_sharing,scaffold_vs_patched,fisher_long_loops}.tsv;
  formatted versions in loop_lengths/NA12878_loops.summary.*.tsv.
- Chain: biological/NA12878/liftover/alignment_results_mm2/chainnet/NA12878.T2T-CHM13.over.chain;
  liftOver run with default settings (minMatch 0.95).
- Junction positions: NA12878.2.cbreak_2.contig-boundaries.bed, derived from
  biological/NA12878/NA12878.2.cbreak_2.contigs.bed (one 1-bp record at each contig end).
- Quartiles use Python statistics.quantiles (default "exclusive" method).
- TO CHECK: exact Juicer Tools / HiCCUPS version and parameters. Matrices were built with
  juicer_tools 1.8.9 (pre, MAPQ 0); loops were called with juicer_postprocessing.sh from
  Juicer 1.6 with defaults.
- Confirmed: scaffold-only Hi-C reads were mapped to a separate N-gap reference, not the
  patched assembly: ASSEMBLY = NA12878.2.cbreak_2.scaffold-only.patched (scaffold-only/Makefile).
- TO CHECK: Hi-C data citation (Rao et al. 2014, GM12878 in situ MboI) and whether the
  read-processing description already appears elsewhere in Methods.
- Known issue: reconstitute_loops.py compares anchor_id (string) to integer 1, so every
  partially shared loop is written with its surviving anchor in the second slot. This does
  not affect category counts.
-->
