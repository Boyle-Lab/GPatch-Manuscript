### Simulation of a contig-level assembly with known contig order and orientation

To benchmark reference-guided patching against a ground truth, we simulated a contig-level assembly from the complete T2T-HG002 genome (hg002v1.1; paternal haplotype chromosomes chr1–22 and chrY, 2.948 Gb). Contig boundaries were modeled on a real contig-level assembly of the same individual: the HGSVC HG002.1 paternal haplotype assembly (hifiasm; 762 contigs, 2.964 Gb, contig N50 85.0 Mb). Breakpoints in real assemblies are not randomly distributed but concentrate in repetitive and structurally complex regions where assembly fails. We therefore derived the simulated contig boundaries from where the HG002.1 contigs align to T2T-HG002, rather than sampling breakpoints at random or from a parametric model of repeat-associated breakage. This way, breakpoint locations, contig overlaps and unfilled inter-contig gaps all follow those of a real assembly of the same genome. We then verified the repeat enrichment of the resulting breakpoints against the T2T-HG002 repeat annotations (below).

HG002.1 paternal contigs were aligned to the T2T-HG002 paternal chromosomes with minimap2 v2.28 (`-x asm5 -c --secondary=no`). Only primary alignments (`tp:A:P`) were retained, with no MAPQ threshold. Consecutive alignment blocks of the same HG002.1 contig were merged into a single block if they were on the same chromosome and strand, colinear in both query and target coordinates, and separated by no more than 100 kb in both coordinate systems (gap or overlap). This avoids introducing artificial breakpoints wherever a contiguous HG002.1 contig is interrupted in the alignment by structural variation or divergent sequence. A merged block spans the full target interval of its component blocks, and blocks that changed strand or order (e.g., inversions) were never merged. In total, 601 blocks were merged away, and 57 of the final blocks comprise two or more original alignment blocks. Blocks spanning less than 1 kb of the benchmark genome were then removed (n = 271). We also removed blocks whose benchmark interval lay entirely within that of another retained block (n = 2,188), because such blocks would have no defined position in the contig order.

Each of the remaining 332 blocks defined one simulated contig. Its sequence was extracted from the T2T-HG002 paternal chromosome at the block's target coordinates and reverse-complemented when the block aligned to the minus strand (142 of 332 contigs). Contigs were named `<chromosome>_<rank>`, where rank is the contig's 0-based position along the chromosome after sorting by start coordinate. Ground-truth locations were recorded in BED format, including strand, MAPQ, source HG002.1 contig and query coordinates, alignment score, edit distance and divergence. Overlaps and gaps between adjacent contigs were taken directly from the alignment boundaries and were not modified. Contig extraction and filtering used custom Python scripts (Python 3.12.8, Biopython 1.84).

The simulated assembly comprises 332 contigs (contig N50 76.5 Mb) covering 2.917 Gb (98.95%) of the T2T-HG002 paternal genome. Contigs were placed on all 23 chromosomes, and 26 HG002.1 contigs contributed more than one simulated contig. Of the 309 junctions between adjacent contigs, 175 (56.6%) are overlaps (median 22.6 kb, range 1.1–903 kb) and 134 (43.4%) are gaps (median 43.8 kb, range 55 bp–2.75 Mb); none abut exactly. Fifty-three junctions (17.2%) separate contigs derived from the same HG002.1 contig, at strand or order changes or at separations exceeding 100 kb.

### Repeat and segmental-duplication enrichment at simulated breakpoints

We tested whether simulated breakpoints are enriched in repetitive sequence using the HG002 Q100 annotations of the T2T-HG002 v1.1 paternal haplotype: segmental duplications (SEDEF, University of Washington; release 013025; all SDs and those with ≥98% identity), centromeric satellite annotation (CenSat v2.0; α-satellite, HSat1–3, other satellites, and rDNA; centromeric transition regions excluded), RepeatMasker (University of Connecticut; aligned blocks only, by repeat class), TRF tandem repeats, and microsatellites (mnr10). Chromosome lengths in the annotations matched the benchmark exactly, so no coordinate conversion was needed. We analyzed two units: (i) contig endpoints, excluding chromosome ends (n = 635), split into model-contig termini (n = 423; the alignment reaches within 1 kb of the end of the HG002.1 contig) and internal alignment breaks (n = 212); and (ii) junctions between adjacent contigs (n = 309; the gap or overlap interval), so that the two correlated endpoints at each junction count as a single unit. For each unit set and category, the test statistic was the fraction of units intersecting the category, with no padding and also with ±1 kb and ±10 kb windows. The null distribution came from 1,000 random placements of the units on the same chromosomes, preserving per-chromosome counts and interval lengths. We report fold enrichment over the permutation mean, two-sided empirical p-values (minimum attainable p ≈ 0.002), and Benjamini–Hochberg q-values within each unit set and window.

Simulated breakpoints were strongly concentrated in repetitive and duplicated sequence. 94.6% of endpoints fell in centromeric satellites or SDs, which together cover 12.7% of the genome. Endpoints were enriched in centromeric satellites (77.9% observed vs 15.5% expected; 5.0-fold), particularly α-satellite (12.4-fold) and rDNA (54-fold), and in segmental duplications (28.7% vs 9.7%; 3.0-fold), most strongly high-identity SDs (≥98%; 4.5-fold). They were also enriched in TRF tandem repeats (4.9-fold) and RepeatMasker satellites (6.3-fold), and depleted in interspersed repeats (LINE 0.23-fold, SINE 0.34-fold, LTR 0.20-fold, DNA transposons 0.16-fold; all q < 0.01). Only 7 of 635 endpoints (1.1%) were outside every annotated repeat or SD, compared with 37.9% expected. Model-contig termini and internal alignment breaks showed the same pattern: CenSat enrichment of 4.1- and 8.4-fold and high-identity SD enrichment of 4.9- and 3.4-fold, respectively. Junctions gave concordant results: α-satellite 11.2-fold, rDNA 16.8-fold, SDs ≥98% 3.3-fold (all q < 0.01).

### Characterization of model-assembly contigs excluded from the simulation

Of the 762 HG002.1 contigs, 251 (2.941 Gb) contributed to the simulated assembly and 511 (22.7 Mb; 0.77% of the assembly) did not. Of these, 420 (19.7 Mb; median 30.1 kb, maximum 640 kb) aligned only within the span of a larger contig and were removed by the containment filter. Of the contained contigs, 208 had a best MAPQ < 10, consistent with repeat-derived fragments, and 146 had MAPQ ≥ 50, consistent with redundant or duplicated sequence in the HG002.1 assembly. Contained contigs therefore add no benchmark sequence not already covered by retained contigs.

The remaining 91 contigs (3.07 Mb; median 28.1 kb, maximum 194 kb) produced no alignment to T2T-HG002. To determine their origin, we aligned them to T2T-CHM13 v2.0 (GCA_009914755.4, including chrX, chrY and chrM) with minimap2 using both the `asm5` and the more permissive `asm20` presets. Only one contig aligned (23.2 kb; CHM13 chr14:100,297,932–100,321,107, MAPQ 60), within the immunoglobulin heavy-chain (IGH) locus. This contig also aligned to T2T-HG002 maternal chr14 (99.8% identity) but not to paternal chr14. This is consistent with paternal IGH sequence that is absent from the paternal T2T-HG002 assembly, as expected after V(D)J recombination in a B-lymphoblastoid cell line. Because HG002 DNA derives from an Epstein–Barr virus (EBV)-transformed lymphoblastoid cell line, and several unaligned contigs carried ~3.07 kb tandem repeats with ~60% GC (characteristic of the EBV major internal repeat, IR1), we aligned the unaligned contigs to the EBV type 1 reference genome (NC_007605.1) with minimap2 (`-x asm20 -c`). Forty-six contigs (1.43 Mb) aligned to EBV over ≥93% of their length. For the remaining 44 contigs (1.62 Mb; GC 18–60%), we searched a 2 kb segment from the midpoint of each against the NCBI core nucleotide database (core_nt) with megablast through the NCBI BLAST web service (Biopython 1.84 `NCBIWWW.qblast`; E-value ≤ 1e-20; top 5 hits). All 44 had high-identity matches (98.7–100% identity over ≥1.1 kb, E ≈ 0) to non-human organisms commonly found as laboratory or reagent contaminants: *Escherichia coli* (17 contigs, 0.46 Mb), *Saccharomyces cerevisiae* (13 contigs, 0.30 Mb; nuclear, rRNA and mitochondrial), and grapevine (*Vitis vinifera*; 14 contigs, 0.86 Mb), comprising the mitochondrial genome (12 contigs), a 194 kb contig matching the chloroplast genome, and one nuclear rDNA contig. Top-5 hits were concordant for every contig at the genus level, apart from related plant mitochondrial or yeast rRNA entries. In summary, 90 of the 91 unaligned contigs are EBV or non-human contaminants, and the single human contig represents IGH sequence absent from the paternal T2T-HG002 haplotype. None of the excluded contigs therefore represents T2T-HG002 paternal sequence missing from the simulation.

### Simulation of a second contig-level assembly from a fragmented benchmark (NA12878)

To test GPatch on a more fragmented assembly, we applied the same procedure to NA12878.
* *Model assembly (NA12878.2):* the HGSVC v12 haplotype-1 assembly (PacBio CCS, Strand-seq phased; `h1-un`; 10,406 contigs, 2.983 Gb, contig N50 18.3 Mb).
* *Benchmark (NA12878.1):* the maternal chromosome-scale assembly (3.015 Gb in 54 sequences).
  * 13 chromosomes are single sequences: chr3, 4, 8, 10–13, 15–19 and X.
  * The other 10 are split into 2–14 sequences of unknown relative order and orientation.
  * All sequences of a chromosome share the same identifier, so we renamed each as `<chromosome>.<descriptor>` (e.g. `chr1.mat-0000003`) before alignment. The sequences were not changed.

Alignment, merging (≤ 100 kb), length and containment filters were as for HG002. Contigs were named `<benchmark sequence>_<rank>`, ranked within each benchmark sequence, and ground truth was recorded in benchmark-sequence coordinates. Junctions are defined only between contigs on the same benchmark sequence.
* *Filtering:* 4,395 blocks were merged away, 1,379 final blocks comprise two or more alignment blocks, 1,277 blocks shorter than 1 kb were removed, and 8,975 contained blocks were removed.
* *Simulated assembly:* 7,400 contigs (contig N50 18.3 Mb; 3,713 reverse-complemented) cover 2.896 Gb (96.1%) of the 48 benchmark sequences that received contigs. The six benchmark sequences without contigs are all < 150 kb.
* *Junctions:* of 7,352 junctions, 4,269 (58.1%) are overlaps (median 8.4 kb, range 3 bp–153 kb), 3,081 (41.9%) are gaps (median 13.4 kb, range 2 bp–1.30 Mb) and 2 abut.
* *Contigs:* most are small (median 16 kb), reflecting the model assembly; 570 model contigs contributed more than one simulated contig.

Breakpoint repeat enrichment was not tested, because NA12878.1 has no repeat annotation.

*Haplotype of the model assembly.* The HGSVC v12 haplotypes are phased with Strand-seq, without parental assignment. To determine which parent h1 represents, we aligned 31 h1 contigs of 1–5 Mb (73 Mb) to the maternal (NA12878.1) and paternal (NA12878.2 chromosome-scale) assemblies with minimap2 (`-x asm5`). We compared the length-weighted gap-compressed divergence of the MAPQ 60 alignments. 26 contigs matched one parent with divergence ≤ 0.0002 and the other at ~0.001. The parent was consistent within each chromosome but varied between chromosomes:
* maternal: chr1, 5, 8, 12, 16, 17, 18 and X;
* paternal: chr2, 4, 6, 7, 10, 11, 14, 15, 19 and 22.

The simulated contig sequences always come from the maternal benchmark. On about half of the chromosomes, however, the breakpoints come from paternal contigs aligned to the maternal benchmark, and may include alignment breaks at heterozygous structural variants.

*Model-assembly contigs excluded from the simulation.* Of the 10,406 contigs, 6,724 (2.924 Gb) contributed to the simulation and 3,682 (59.0 Mb; 1.98%) did not.

Contained contigs, 3,667 (58.8 Mb; median 14.3 kb, maximum 181 kb), account for almost all of them:
* *Containing interval:* 2,367 (35.2 Mb) lie within a single alignment block of a retained contig. The other 1,300 (23.6 Mb) lie within a span bridged by merging.
* *Satellite content:* 2,789 (46.8 Mb) are satellite-like by sequence (21-mer diversity < 0.5, or ≥ 20% HSat2/3 `CATTC`/`GAATG` motif).
* *Location:* they cluster in pericentromeric heterochromatin and acrocentric arms (e.g. chr9q12, chr1q12, chr13/14/22 short arms). 84 one-megabase windows contain 83% of their alignment blocks.
* *Mapping quality:* best MAPQ was < 10 for 1,716 and ≥ 50 for 731.
* *Parental origin:* aligning them to both parental assemblies gave mixed results. Of the 878 non-satellite contigs, 232 matched the maternal assembly better, 196 the paternal, 174 both equally and 276 neither at ≥ 99% identity.
As for HG002, contained contigs add no benchmark sequence not already in the simulation.

Fifteen contigs (0.28 Mb) produced no alignment (4) or only alignment blocks < 1 kb (11). We aligned them to T2T-CHM13, to the paternal assembly and, with `-x asm20`, to the maternal benchmark. None aligned to EBV (NC_007605.1; `-x asm20`), and all 15 were human, so no BLAST search was needed:
* Five (89 kb) are paternal-specific: they match the paternal assembly at ≥ 99.9% identity but are absent from or highly divergent in the maternal benchmark. Two of these (on chr5 and chr18) barely align to CHM13 and represent non-reference insertions.
* Nine (121 kb) are centromeric or acrocentric satellites by their CHM13 positions and sequence composition. They include chr4 42 bp AT-rich and chr16 GC-rich satellites, chr8 and chr9 α-satellite, chr14/chr21 short-arm sequence and an HSat3 array absent from both NA12878 assemblies. They align to the benchmark at 72–94% identity or not at all.
* One (66 kb) is unique sequence (CHM13 chr7:60.36 Mb, MAPQ 60) that falls in the gap between the two chr7 benchmark sequences.

The NA12878 model assembly thus contained no EBV or contaminant contigs. Its excluded contigs reflect redundant satellite fragments, other-haplotype sequence, and sequence missing from the benchmark.

### Junction handling: overlap trimming, gap patching and bookending

After contigs are placed and ordered along each reference chromosome, GPatch (v0.6.0) decides how to join each pair of adjacent contigs. If the two contigs overlap, it trims the 5′ end of the downstream contig. If there is a gap between them, it fills the gap with reference sequence ("patching"). Otherwise it joins them directly ("bookending"). Patches are written in lower case, and no N characters are introduced except in scaffold-only mode. GPatch versions up to 0.4.0 derived each contig's extent on the reference from its primary alignment, extended at each end by the full length of the soft-clipped sequence, and trimmed or patched by the overlap or distance between those extents. That extrapolation assumes unaligned contig sequence spans the same length of reference as of contig. The assumption fails when the unaligned sequence lies in a centromere, segmental duplication or structural variant, where it can be hundreds of kilobases to megabases long. We therefore replaced it with a procedure that estimates contig extents more conservatively and measures overlaps directly from the contig sequences.

*Contig extents.* The primary alignment of each placed contig is chained to that contig's supplementary alignments. A supplementary alignment is added if it has mapping quality ≥ 20, an aligned length ≥ 10 kb, the same chromosome and strand as the primary alignment, and is colinear with the growing chain in both contig and reference coordinates. Consecutive chained alignments may be separated by at most 1 Mb, or overlap by at most 10 kb, in both coordinate systems. Each contig end is then extrapolated onto the reference only by the residual unaligned sequence beyond the outermost chained alignment. The gap or overlap between adjacent contigs is estimated from these extents, as the distance from the end of the sequence placed so far to the extrapolated start of the next contig.

*Overlap measurement.* Overlaps are detected with an anchored suffix–prefix test on the contig sequences, in placement orientation. For each junction, GPatch searches a window at the 3′ end of the upstream contig and at the 5′ end of the downstream contig. The window starts at the larger of 1 Mb and twice the estimated overlap plus the residual unaligned lengths. It is doubled, up to 4 Mb, when an anchor match lies at the window's edge. Within these windows, GPatch finds exact matches of two 500 bp anchors:

- the first 500 bp of the downstream contig, located within the upstream window, which imply an overlap length ending at the upstream contig's end;
- the last 500 bp of the upstream contig, located within the downstream window, which imply an overlap length starting at the downstream contig's start.

A candidate overlap is any length supported by both anchors. It is accepted if the full overlapping sequences are ≥ 99.9% identical. Candidates for which the two anchor-implied lengths differ by ≤ 50 bp (overlaps containing indels) are verified by aligning the two overlapping segments with minimap2 (`-x asm5`), with at most 20 such candidates tested per junction. We did not use a general-purpose aligner to search for overlaps. In palindromes and tandem arrays, a short true overlap is reported as a low-scoring secondary of a longer alignment between paralogous copies, and is discarded by the aligner's secondary-alignment filter. In an earlier aligner-based implementation, this caused true overlaps of up to 101 kb to be missed. The minimum detectable overlap is the anchor length.

*Decision rules.* The extent-based estimate at a junction is classed as:

- *trusted* when the residual unaligned sequence at both contig ends is ≤ 10 kb;
- *one-sided* when only one end exceeds this;
- *untrusted* when both do.

GPatch then decides as follows:

1. If exactly one overlap is found, it is trimmed (`trim_measured`).
2. If several distinct overlap lengths are found, as occurs in tandem repeats (`trim_ambiguous`):
   - with a trusted or one-sided estimate, the candidate closest to the estimated overlap is trimmed;
   - otherwise the smallest candidate is trimmed, which never over-trims when the true overlap is among the candidates.
3. If no overlap is found and the estimated gap is non-negative, the gap is patched with the corresponding reference interval: `patch` with a trusted estimate, `patch_one_sided` with a one-sided one.
4. Otherwise the contigs are bookended.

Leading and trailing chromosome ends are patched when the residual at the terminal contig end is ≤ 10 kb. Contigs are never trimmed against a chromosome end. Every decision is recorded in a per-junction report (`junctions.tsv`), with the extents, residual lengths, estimate, candidate overlaps and the reason for the decision. The previous behavior remains available as `--legacy_junctions`.

### Benchmarking junction handling against the simulated assembly

To evaluate junction handling, we aligned the simulated contigs to T2T-CHM13 v2.0 with minimap2 v2.28 (`-x asm20 -a`). We then ran GPatch v0.6.0 (`-d -m 10`) in three configurations:

1. `--legacy_junctions`, which reproduces GPatch 0.4.0 with soft-clip extrapolation and trimming enabled (outputs byte-identical to GPatch 0.4.0);
2. the revised junction handling without overlap extension (`--no_extension`);
3. the default, which adds overlap extension (see the next section).

The first two runs differ in only seven contigs: the revised nested-contig check placed 155 contigs rather than 148, and contig order and orientation are otherwise the same. The default run placed 176 contigs. Because the runs place different contigs, the number of junctions in each class differs between them (Table X).

Each junction was compared with the truth from the simulation:

- *Adjacent in the simulation:* compared with the recorded overlap or gap.
- *Separated by unplaced contigs:* the true distance was computed from the two contigs' benchmark coordinates.
- *Chromosome ends:* compared with the benchmark sequence before the first or after the last placed contig.
- *Contigs from different benchmark chromosomes:* excluded (12 in each run).

For each junction, the signed error is the GPatch distance (patch length, or minus the trim length) minus the true distance (gap length, or minus the overlap length). Positive errors mean extra sequence (duplicated overlap or excess patch); negative errors mean missing sequence. T2T-CHM13 v2.0 chrY is the HG002 chrY assembly, so chrY junctions are trivially accurate. We therefore report autosomal results; Table X gives all three configurations.

*Overlaps.* The revised junction handling was compared with soft-clip extrapolation (configurations 1 and 2).

- *Trimmed exactly:* 24 of 33 true overlaps, against 0 of 31.
- *Detected:* every true overlap except one, at a contig GPatch placed in the opposite orientation to the benchmark, where no overlap exists in placement orientation.
- *Over-trimming,* which deletes real sequence: 1.61 Mb at 11 junctions under soft-clip extrapolation, against 11.3 kb at 4.
- *Under-trimming,* which duplicates sequence: similar for both (18.2 kb at 14 junctions against 18.8 kb at 4).

Under soft-clip extrapolation, the largest errors arose where unaligned contig ends were extrapolated across sequence whose length differs between HG002 and CHM13. At the chr10 centromere, the two contigs' unaligned ends (3.68 and 3.13 Mb) were extrapolated across 5.82 Mb of CHM13 sequence, producing a 993 kb trim at a true 2.6 kb overlap. Soft-clip extrapolation also patched 6 true overlaps (2.9–903 kb), duplicating 5.6 Mb. It trimmed 10 true gaps by 1.6–619 kb; the revised procedure trimmed 3, by 1.7–8.0 kb, through false overlaps in tandem arrays.

*Gaps.* Of the 18 true gaps patched by the revised procedure, 9 were within 1 kb of the true gap length (median absolute error 1.5 kb). Under soft-clip extrapolation, 2 of 11 were (median 28 kb). Most remaining patch-length errors reflect genuine differences in length between the HG002 gap sequence and the CHM13 sequence used to fill it, such as an 845 kb excess at chr7_7/chr7_8.

*Totals.* Across all autosomal junctions and chromosome ends:

| Measure | Soft-clip extrapolation | Revised procedure |
|---|---|---|
| Extra sequence | 29.2 Mb | 10.3 Mb |
| Missing sequence | 37.7 Mb | 33.4 Mb |
| Total absolute junction error | 66.8 Mb | 43.7 Mb |

Most of the remaining error (41.8 Mb) lies at junctions flanking contigs that GPatch did not place. There, the true distance includes the sequence of the unplaced contigs, so junction handling cannot recover it.

*With overlap extension (configuration 3):*

- *Accuracy:* 39 of 49 true overlaps were trimmed exactly, and total junction error fell to 42.6 Mb, with less missing sequence (32.3 Mb).
- *New true gaps:* extension turns some junctions that previously flanked unplaced contigs into junctions between true neighbours. Some of these gaps are patched with large errors, for example 1.10 Mb at chr14_14/chr14_15 and 316 kb at chr9_3/chr9_4. Positions carried through overlaps from the walk's starting contig accumulate the HG002–CHM13 length difference across the walk, and those errors now appear under true gaps.
- *Decisions:* of the 153 junctions between adjacent placed contigs, GPatch trimmed 82 (70 single overlaps and 12 ambiguous), patched 51 (37 with a trusted estimate and 14 one-sided) and bookended 20.
- *Sequence check:* within ±20 kb of each of the 60 exactly trimmed junctions, the patched sequence was identical to T2T-HG002.

**Table X.** Junction accuracy on the simulated assembly (autosomes), GPatch v0.6.0. Configuration 1 is `--legacy_junctions` (equivalent to GPatch 0.4.0 with trimming); 2 is `--no_extension`; 3 is the default. n is the number of junctions in each class, which differs between configurations because they place different contigs. Exact and within-1 kb are junction counts; extra, missing and total error are summed over junctions.

| Junction class | Configuration | n | Exact | ≤ 1 kb | Extra (bp) | Missing (bp) | Total error (bp) |
|---|---|---|---|---|---|---|---|
| True overlap | 1 | 31 | 0 | 10 | 5,579,963 | 1,612,778 | 7,192,741 |
| | 2 | 33 | 24 | 26 | 119,773 | 11,276 | 131,049 |
| | 3 | 49 | 39 | 41 | 129,969 | 11,276 | 141,245 |
| True gap | 1 | 21 | 0 | 2 | 1,388,600 | 1,979,288 | 3,367,888 |
| | 2 | 24 | 0 | 9 | 1,052,053 | 316,759 | 1,368,812 |
| | 3 | 32 | 0 | 9 | 2,467,034 | 708,465 | 3,175,499 |
| Chromosome end | 1 | 39 | 7 | 14 | 1,459,457 | 499,657 | 1,959,114 |
| | 2 | 39 | 11 | 26 | 297,671 | 65,110 | 362,781 |
| | 3 | 39 | 11 | 26 | 297,671 | 65,110 | 362,781 |
| Flanking unplaced contigs | 1 | 32 | 0 | 0 | 20,722,546 | 33,566,957 | 54,289,503 |
| | 2 | 33 | 1 | 2 | 8,877,639 | 32,963,786 | 41,841,425 |
| | 3 | 27 | 1 | 2 | 7,462,658 | 31,505,801 | 38,968,459 |
| All | 1 | 123 | 7 | 26 | 29,150,566 | 37,658,680 | 66,809,246 |
| | 2 | 130 | 37 | 64 | 10,347,136 | 33,356,931 | 43,704,067 |
| | 3 | 148 | 52 | 79 | 10,357,332 | 32,290,652 | 42,647,984 |

### Contig placement: unplaced contigs, nested-contig detection and overlap-based extension

*Unplaced contigs.* GPatch places only contigs whose primary alignment has mapping quality at or above a threshold (`-m`; 10 here) and that are not nested. A contig is nested when its reference extent lies entirely within that of a contig sorted before it. GPatch now records every contig left out of the pseudoassembly, with the reason (unmapped, low mapping quality, outside whitelist regions, or nested, with the containing contig), in a report (`unplaced.tsv`). Their sequences are written in original orientation to `unplaced.fasta`, so no assembled sequence is lost silently.

Of the 332 simulated contigs, the original implementation left 184 (26.3 Mb; 0.9% of the benchmark) unplaced; none were unaligned.

- *Low mapping quality: 161 contigs, 22.4 Mb.*
  - These were small (median 59 kb) and almost entirely centromeric: 152 were ≥ 50% CenSat.
  - Their primary alignments were complete (median aligned fraction 100%) but had mapping quality 0–9, indicating multiple equally good placements within satellite arrays.
  - 131 nevertheless aligned to the correct chromosome between the CHM13 positions of their placed benchmark neighbours.
- *Nested: 23 contigs, 3.85 Mb.* 17 of these were nested only because the containing contig's extent had been extrapolated by its soft-clipped ends (up to 27.3 Mb). Their own alignments did not overlap the containing contig's aligned interval. The remaining 6 were nested by their aligned coordinates, mostly through acrocentric short-arm homology.

*Nested-contig detection.* GPatch now identifies nested contigs using the chained extents used for junction handling (primary alignment plus colinear supplementary alignments, extrapolated only by the residual unaligned sequence). Previously it used extents extended by the full soft-clips; that behavior is retained with `--legacy_junctions`. This placed 155 contigs instead of 148 without changing the order of the others. Total junction error was unchanged (43.71 vs 43.69 Mb, all chromosomes).

*Mapping-quality threshold.* Lowering the threshold recovered contigs at the cost of accuracy.

| Threshold | Contigs placed | Extra sequence placed | Autosomal junction error |
|---|---|---|---|
| `-m 10` (default) | 148 | — | 43.7 Mb |
| `-m 1` | 184 | +7.3 Mb | 45.0 Mb |
| `-m 0` | 191 | +7.7 Mb | 52.3 Mb |

The extra error comes from joins between chromosomes and misoriented contigs. Within satellite arrays, the CHM13 alignment of a low-mapping-quality contig does not reliably determine its order relative to other contigs, because HG002 and CHM13 arrays differ in structure. We therefore retained `-m 10`.

*Overlap-based extension.* Contigs that cannot be placed by alignment can often be placed by sequence overlap with a placed neighbour. Of the 184 unplaced contigs, 141 (20.9 Mb) overlapped at least one benchmark neighbour by ≥ 500 bp. We tested every unplaced contig, in both orientations, against all other contigs with the anchored overlap test (above). The results were:

- every true overlap was detected;
- 55 contigs (10.2 Mb) had only correct overlap partners;
- 86 (10.6 Mb) also had false partners in tandem arrays;
- requiring consistency with each contig's CHM13 alignment (same chromosome, within 1 Mb) reduced false overlaps from 6,966 to 334.

On this basis, GPatch (v0.6.0) includes a conservative extension step that runs after alignment-based placement. It is enabled by default and disabled with `--no_extension`. Contigs left out for low mapping quality or nesting are candidates. For each gap between adjacent placed contigs, and for each chromosome end, walks extend from both flanking contigs into the gap. At each step, an unplaced contig is added if:

1. it overlaps the current contig end by ≥ 5 kb;
2. it is the only unplaced contig to do so;
3. its primary alignment is on the same chromosome, within 1 Mb of the gap, and in the orientation in which it would be placed.

A walk stops when no candidate or more than one qualifies, or when it reaches the contig on the other side of the gap. Contigs claimed by more than one walk are rejected. The reference coordinates of inserted contigs are derived from the placed contig at the start of the walk, through the measured overlaps, rather than from their own alignments. These positions are used to estimate gaps, and inherit that contig's residual unaligned length for the trust rules. Every walk step, stop reason and rejected contig is recorded in a report (`extension.tsv`). Inserted contigs are included in the placement and junction outputs and removed from the unplaced-contig report.

Two design choices were essential.

- *Position from overlaps, not alignments.* Using the inserted contigs' own low-mapping-quality alignments produced large errors in the estimated gaps at the ends of walks, for example a 3.0 Mb patch for a 79 kb gap. Total junction error then rose to 48.7 Mb.
- *Minimum overlap of 5 kb.* With a 500 bp minimum, a 585 bp false overlap in a chrY palindrome inserted a contig out of order. All overlaps used for insertion were 7.2–136 kb (median 19.4 kb).

Starting from the 155 contigs placed by GPatch, extension inserted 21 contigs (5.4 Mb), all in the correct relative order. As a result (all chromosomes):

| Measure | Before extension | After extension |
|---|---|---|
| Contigs placed | 155 | 176 |
| Junctions between benchmark neighbours | 88 | 116 |
| Exactly trimmed overlaps | 41 of 50 | 59 of 69 |
| Total junction error | 43.71 Mb | 42.66 Mb (42.65 Mb on autosomes) |

Of the 254 walk ends, 194 stopped for lack of a qualifying candidate, 4 at ambiguous candidates, and 56 on reaching the opposite flanking contig. The remaining unplaced sequence (20.7 Mb) lies almost entirely in centromeric arrays. There, walks stop at ambiguous candidates or false closures, where short false overlaps between flanking contigs make a gap appear closed. Some candidates are also excluded because their alignment is on another chromosome or in the opposite orientation.

### Unplaced contigs after GPatch 0.6.0 placement (HG002 and NA12878)

For both simulations, we aligned the simulated contigs to T2T-CHM13 v2.0 (minimap2 `-x asm20 -a`) and ran GPatch v0.6.0 with default settings (`-d -m 10`, revised junction handling and overlap extension). We then characterized every contig in the unplaced-contig report.

NA12878.1 has no repeat annotation, so for both genomes repeat context was taken from two sources:
* the T2T-CHM13 v2.0 CenSat (v2.1, excluding centromeric transition entries) and SD annotations, over each contig's primary CHM13 alignment;
* the contig sequence itself: 21-mer diversity < 0.5, or ≥ 20% HSat2/3 motif.

For each unplaced contig we also recorded:
* *location:* whether the primary alignment lies on the correct chromosome between the CHM13 positions of its nearest placed neighbours from the same benchmark sequence (± 2 Mb);
* *overlaps:* its true overlaps with benchmark neighbours;
* *nesting type:* for nested contigs, whether the contig is nested by its aligned coordinates or only by the containing contig's chained extent.

To understand the limits of overlap extension, we did two further analyses:
* *Walk stops:* for each walk that stopped without a candidate, we classified the true benchmark neighbour of its last contig.
* *Reachability:* we determined which unplaced contigs can be reached from a placed contig of the same benchmark sequence through a chain of true overlaps ≥ 5 kb.

| Measure | HG002 | NA12878 |
|---|---|---|
| Simulated contigs | 332 | 7,400 |
| Placed (of which by overlap extension) | 176 (21) | 2,704 (296) |
| Unplaced | 156 (20.7 Mb; 0.71%) | 4,696 (93.1 Mb; 3.17%) |
| Low mapping quality | 141 (17.1 Mb; median 49 kb) | 4,516 (89.6 Mb; median 15 kb) |
| Nested | 15 (3.6 Mb) | 174 (3.5 Mb) |
| Unaligned | 0 | 6 (13 kb) |
| Low-MAPQ contigs ≥ 50% CenSat/SD at CHM13 alignment | 141 of 141 | 4,515 of 4,516 |
| Low-MAPQ contigs between placed neighbours | 117 | 3,207 |
| Low-MAPQ contigs on another chromosome | 22 | 326 |
| Unplaced, reachable through overlaps ≥ 5 kb | 37 (6.4 Mb) | 376 (11.6 Mb) |
| Unplaced, cut off by gaps or overlaps < 5 kb | 119 (14.3 Mb) | 4,030 (76.8 Mb) |
| Overlap islands of unplaced contigs (≥ 100 kb span) | 87 (33; 15.9 Mb) | 2,644 (110; 19.4 Mb) |

In both genomes, low mapping quality accounts for most unplaced contigs, and these contigs lie almost exclusively in centromeric satellites or SDs. Their alignments are complete (median aligned fraction 1.00) but not unique.

NA12878 differs in scale. Its model assembly breaks centromeric and pericentromeric satellite arrays into thousands of ~15 kb contigs, most prominently in chr1q12 (1,020 unplaced contigs on chr1), chr9, chr16, chr15, chr3 and chr19.
* Thirteen small benchmark sequences on chr1, chr5, chr6 and chr20 have no placed contig at all (290 contigs, 4.8 Mb).
* Lowering the mapping-quality threshold would misplace many of these contigs: 326 low-MAPQ contigs align to another chromosome and 698 elsewhere on the correct chromosome.

*Nested contigs.*
* *HG002:* 4 of the 15 are nested by their aligned coordinates, all through acrocentric short-arm homology. The other 11 are nested only by the containing contig's chained extent:
  * six lie within the span of the container's high-quality alignments. They include chr1_3 (2.54 Mb, MAPQ 60), nested within the extent of the 19.2 Mb chr1q12 contig chr1_2, of which only 1.08 Mb is in the primary alignment.
  * five are nested only through residual extrapolation.
* *NA12878:* 134 of the 174 are nested by aligned coordinates and 40 only by chained extent.

*Limits of overlap extension.* Overlap extension can only walk inward from placed contigs. Within runs of unplaced contigs, a third of NA12878 junctions are true gaps (1,398 of 4,193) and 743 are overlaps < 5 kb.
* *NA12878:* only 376 unplaced contigs (11.6 Mb) are reachable from a placed contig through true overlaps ≥ 5 kb, although 2,892 (68.1 Mb) overlap a true neighbour by ≥ 5 kb. Of 4,240 walk ends:
  * 577 reached the flanking contig;
  * 3,363 stopped at a true gap or overlap < 5 kb, of which 2,597 were beside an already placed neighbour, where there was nothing to insert;
  * 59 stopped at ambiguous candidates;
  * 72 stopped because the alignment-consistency rules excluded a truly overlapping neighbour (46 outside the gap region, 17 on another chromosome, 9 in the opposite orientation).
* *HG002:* of 254 walk ends, 56 reached the flanking contig, 130 met a true gap or overlap < 5 kb, 41 reached a chromosome end and 9 were blocked by the alignment rules.

In both genomes, the remaining unplaced sequence forms many small overlap islands, separated from placed contigs by gaps. Placing it would require assembling these islands and positioning them as units, rather than extending further from placed contigs.

### Evaluating against a fragmented benchmark

Because the relative order and orientation of NA12878.1 sequences from the same chromosome are unknown, contig order, contig placement and sequence accuracy are evaluated within each benchmark sequence. Each benchmark sequence may lie in either orientation relative to the patched genome (`compare_contig_ordering_split.py`, `compare_contig_positions_split.py` and `compare_fasta_by_editdist_split.py`).

For each benchmark sequence:
* *Orientation:* + if the length-weighted majority of its placed contigs lie on the same strand as in the simulation, − otherwise. For − sequences, the expected strand of each contig is reversed, as is the expected order of its contigs, and the patched sequence is reverse-complemented before comparison.
* *Target:* the patched chromosome holding most of its placed contig sequence.

*Ordering.* Edit distance, adjacent-pair recovery and switch errors compare the benchmark order of the sequence's contigs with their order on the target. The comparison spans from the first to the last of these contigs. Contigs from other benchmark sequences placed between them count as errors; contigs outside the span are not compared.

*Placement.* The exact-coordinate test of the original scripts is replaced by a relative one. A contig is correctly placed if:
* it lies on the target chromosome named after its benchmark chromosome; and
* its start and end offsets (patched minus benchmark coordinate, after orientation) both equal the most frequent offset among contigs of the same benchmark sequence.

Order is checked against the previous placed contig of the same benchmark sequence.

*Sequence.* A benchmark sequence that is the only one for its chromosome is compared with the whole patched chromosome, as for HG002. Fragmented ones are compared with the patched interval spanning their placed contigs.

Benchmark sequences with no placed contigs are skipped, as unplaced chromosomes are in the original scripts. On the HG002 simulation, whose benchmark is not fragmented, the per-benchmark-sequence versions reproduce the original edit distances, orientation calls and sequence distances. Results are unchanged when a patched chromosome is reversed (coordinates, strands and sequence).

*Contig mapping accuracy.* Each simulated contig was scored as mapped (placed in the patched genome), correctly mapped (relative-offset test above), on the correct chromosome, in the correct orientation and in the correct order. Mapping statistics were summarized with `compile_mapping_stats_split.Rscript`, which computes clustering accuracy from the chromosome part of the benchmark sequence ID. For HG002 we report both the per-benchmark-sequence definitions and the original exact-coordinate definitions; clustering and orientation accuracy are identical under both.

| Measure (GPatch v0.6.0, default) | HG002 | NA12878 |
|---|---|---|
| Simulated contigs | 332 | 7,400 |
| Contigs mapped (placed) | 176 (53.0%; 2,905.8 Mb) | 2,704 (36.5%; 2,845.2 Mb) |
| Clustering accuracy (correct chromosome) | 0.943 | 0.950 |
| Orientation accuracy | 0.960 | 0.926 |
| Correct order | 0.909 (original definition: 0.636) | 0.878 |
| Correctly mapped (common offset) | 0.170 (original exact coordinates: 0.023) | 0.034 |

Fractions are of mapped contigs. Clustering and orientation accuracy are similar for the two genomes. Although NA12878 has more than 20 times as many contigs, only 5.0% of its placed contigs are on the wrong chromosome and 7.4% in the wrong orientation (HG002: 5.7% and 4.0%). Placed contigs cover similar amounts of sequence in both genomes, because the unplaced NA12878 contigs are small.

The correctly-mapped fraction is low in both genomes. That test requires every junction between a contig and the most common offset of its benchmark sequence to be exact. A single mis-sized patch or trim upstream shifts all subsequent contigs, so the junction-level evaluation above is the more informative measure of positional accuracy.

Order accuracy is lower for NA12878 than for HG002. Many NA12878 contigs are short satellite fragments placed by overlap extension or by low-uniqueness alignments within arrays, where order is hardest to resolve.

<!--
Notes for authors (delete before submission):
- Enrichment p-values assume independent units. Breakpoints cluster within the same centromeres
  and SD blocks, so p-values are likely anti-conservative; the fold enrichments are the robust
  quantity. The junction-level test partly addresses this.
- BLAST identification of the 44 contaminant contigs is based on one 2 kb segment per contig.
- The IGH interpretation is inferred, not tested.
- Junction benchmark baseline: the "soft-clip" run used GPatch 0.4.0 with trimming enabled. The
  released 0.4.0 default disabled trimming (--no_trim on by default), which bookends overlapping
  contigs untrimmed; that configuration was not benchmarked here. Decide which baseline the paper
  should report.
- Several "errors" next to the 5 contigs placed in the opposite orientation to the benchmark
  (chr2_6, chr7_5, chr9_21, chr17_9, chr20_1; centromeric/pericentromeric or 17q21.31) reflect the
  orientation, not junction handling (e.g., the paired ±0.8–1.4 kb trims and the missed chr17_8/9
  overlap). Placement accuracy should be reported separately.
- Thresholds (10 kb residual, 500 bp anchors, 99.9% identity) were chosen on this simulation, the
  same data used to evaluate them; consider validating on an independent assembly.
- GPatch versions. All changes described here are in GPatch 0.6.0 (commit fda475d): revised
  junction handling (commit b02e259), unplaced-contig report and chained-extent nested check
  (9008618; also in 0.5.0), and overlap extension (e72d0cb; first released in 0.6.0, on by
  default, --no_extension disables it). The soft-clip baseline is --legacy_junctions, which
  reproduces GPatch 0.4.0 outputs byte-for-byte when trimming is enabled (the 0.4.0 default
  disabled trimming; see the baseline note above).
- Which run produced which numbers:
  - Table X and the junction-benchmark text: GPatch 0.6.0 (fda475d) run as --legacy_junctions,
    --no_extension and default (outputs in GPatch_v0.6.0/; table source
    GPatch_v0.6.0/T2T-HG002.HG002.1.simulated.merge100000.benchmark-table.txt). The legacy run is
    byte-identical to the original GPatch 0.4.0 (+trimming) run; the other two are identical to
    the validated test runs.
  - Placement section: 155 contigs = 0.6.0 --no_extension; 176 contigs = 0.6.0 default. Both
    were verified identical to the prototypes (overlap_extension_prototype.py v0.1.0 +
    gpatch_junction_prototype.py --propagated) on these data (identical contigs.bed,
    patches.bed and junction decisions).
  - Evaluation script: compare_overlap_trimming.py v0.2.0.
- The 5 kb minimum overlap for extension was chosen after observing one spurious 585 bp overlap
  on this same simulation.
- The -m 1 / -m 0 junction errors exclude junctions between contigs from different benchmark
  chromosomes, so they understate the cost of lowering the threshold.
- NA12878: the h1 parent assignment is from a 31-contig sample; the per-chromosome pattern was
  consistent but not every chromosome was sampled (chr3, 9, 13, 20, 21 untested). A
  parent-consistent model (h1 or h2 per chromosome, matched to the maternal benchmark) would make
  the NA12878 simulation equivalent to HG002's; not done.
- NA12878 repeat context uses CHM13 annotations at each contig's CHM13 alignment (and sequence
  composition), not benchmark annotations; low-MAPQ contigs may align to the wrong copy of an
  array, but the repeat class is unlikely to change. For HG002, the earlier unplaced-contig
  numbers (184 contigs, 152 ≥ 50% CenSat) used HG002 benchmark annotations at true coordinates and
  the original GPatch placement; the 0.6.0 comparison table uses CHM13 annotations for both genomes.
- HG002 chr1_3 (2.54 Mb, MAPQ 60) is dropped as nested in chr1_2's chained extent; worth
  revisiting the nested check (e.g. require the nested contig's alignment to lie within the
  container's chained aligned blocks, not the extrapolated residual).
- NA12878 -m 1 / -m 0 runs and an island-based placement strategy were not tested.
- Mapping statistics: NA12878 from compare_contig_positions_split.py +
  compile_mapping_stats_split.Rscript (T2T_NA12878/GPatch/*.mapping-stats.txt, *.mapping-summary.txt).
  HG002 "per-benchmark-sequence" values were computed with compare_contig_positions_split.py on
  T2T-HG002/GPatch (not written to the results folder); the "original definition" values are from
  T2T-HG002/GPatch/*.mapping-stats.txt (compare_contig_positions.py). The original order test
  compares a contig's patched start with the previous reference contig's benchmark end (mixing
  coordinate frames, and counting unplaced predecessors), which is why it gives 0.636.
- NA12878 code/results: simulation/T2T_NA12878/ (rename_benchmark_contigs.sh,
  simulate_contigs_from_alignments.sh -M 100000, characterize_unplaced_contigs.py v0.1.1,
  analyze_extension_limits.py v0.1.0; excluded_contigs/; GPatch/*.unplaced-characterization.*,
  GPatch/*.extension-limits.txt). HG002 0.6.0 unplaced results: T2T-HG002/GPatch/. Split
  comparison scripts: simulation/*_split.py + split_benchmark.py v0.1.0.
- Code: simulation/T2T-HG002/simulate_contigs_from_alignments.sh (-M 100000),
  paf_to_simulated_contigs.py v0.2.0, get_hg002_annotations.sh, breakpoint_repeat_enrichment.py
  v0.1.0 (seed 20261109). Results: T2T-HG002.HG002.1.simulated.merge100000.breakpoint-enrichment.tsv
-->
