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

### Junction handling: overlap trimming, gap patching and bookending

After contigs are placed and ordered along each reference chromosome, GPatch decides how to join each pair of adjacent contigs. If the two contigs overlap, it trims the 5′ end of the downstream contig. If there is a gap between them, it fills the gap with reference sequence ("patching"). Otherwise it joins them directly ("bookending"). Patches are written in lower case, and no N characters are introduced except in scaffold-only mode. Earlier versions of GPatch derived each contig's extent on the reference from its primary alignment, extended at each end by the full length of the soft-clipped sequence, and trimmed or patched by the overlap or distance between those extents. That extrapolation assumes unaligned contig sequence spans the same length of reference as of contig. The assumption fails when the unaligned sequence lies in a centromere, segmental duplication or structural variant, where it can be hundreds of kilobases to megabases long. We therefore replaced it with a procedure that estimates contig extents more conservatively and measures overlaps directly from the contig sequences.

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

To evaluate junction handling, we aligned the simulated contigs to T2T-CHM13 v2.0 with minimap2 v2.28 (`-x asm20 -a`) and ran GPatch (`-d -m 10`) with both the previous (soft-clip) and the revised junction handling. The two runs place the same contigs in the same order and orientation, so they differ only at junctions. For this comparison, nested contigs were identified with soft-clip extents in both runs, as in the original implementation (see the next section for the revised nested-contig check). GPatch placed 148 of the 332 simulated contigs, giving 125 junctions between adjacent placed contigs plus 46 chromosome ends.

Each junction was compared with the truth from the simulation. Junctions between contigs that are adjacent in the simulation (83) were compared with the recorded overlap or gap. For pairs separated by unplaced contigs (35), the true distance was computed from the contigs' benchmark coordinates. Chromosome ends (41) were compared with the benchmark sequence preceding the first or following the last placed contig. Junctions between contigs from different benchmark chromosomes (12) were excluded. For each junction, the signed error is the GPatch distance (patch length, or minus the trim length) minus the true distance (gap length, or minus the overlap length). Positive errors represent extra sequence (duplicated overlap or excess patch) and negative errors represent missing sequence. Because T2T-CHM13 v2.0 chrY is the HG002 chrY assembly, chrY junctions are trivially accurate, so we report results for the autosomes (123 junctions with truth).

The revised procedure trimmed 30 of 31 true overlaps, 22 of them exactly, compared with none exactly under soft-clip extrapolation (Table X). Over-trimming, which deletes real sequence, fell from 1.61 Mb at 11 junctions to 11 kb at 4. Under-trimming, which duplicates sequence, was similar for the two methods (18.8 kb at 4 junctions vs 18.2 kb at 14). Under soft-clip extrapolation, the largest errors arose where unaligned contig ends were extrapolated across sequence whose length differs between HG002 and CHM13. At the chr10 centromere, the two contigs' unaligned ends (3.68 and 3.13 Mb) were extrapolated across 5.82 Mb of CHM13 sequence. This produced a 993 kb trim at a true 2.6 kb overlap. Soft-clip extrapolation also patched 6 true overlaps (2.9–903 kb), duplicating 5.6 Mb, and trimmed 10 true gaps by 1.6–619 kb. The revised procedure detected all but one of the true overlaps. The exception was a contig that GPatch placed in the opposite orientation to the benchmark, so no overlap exists in placement orientation.

Patch lengths were also more accurate (Table X). Of the 16 true gaps the revised procedure patched, 9 were within 1 kb of the true gap length (median absolute error 127 bp). Of the 11 true gaps patched under soft-clip extrapolation, 2 were (median absolute error 28 kb). Most remaining patch-length errors reflect genuine differences in length between the HG002 gap sequence and the CHM13 sequence used to fill it, such as an 845 kb excess at chr7_7/chr7_8. Across all 123 autosomal junctions and chromosome ends, extra sequence fell from 29.2 Mb to 10.3 Mb and missing sequence from 37.7 Mb to 33.3 Mb. The total absolute junction error fell from 66.8 Mb to 43.7 Mb. Most of the remaining error (41.8 Mb) lies at junctions flanking contigs that GPatch did not place. There, the true distance includes the sequence of the unplaced contigs, so it cannot be recovered by junction handling.

Of the 125 junctions, the revised procedure trimmed 54 (45 single overlaps and 9 ambiguous), patched 50 (36 trusted and 14 one-sided) and bookended 21. Within ±20 kb of each of the 40 exactly trimmed junctions, the patched sequence was identical to T2T-HG002.

**Table X.** Junction accuracy on the simulated assembly (autosomes). Exact and within-1 kb counts; extra and missing sequence summed over junctions.

| Junction class | n | Method | Exact | ≤ 1 kb | Extra (bp) | Missing (bp) |
|---|---|---|---|---|---|---|
| True overlap | 31 | Soft-clip | 0 | 10 | 5,579,963 | 1,612,778 |
| | | Revised | 22 | 24 | 119,773 | 11,276 |
| True gap | 21 | Soft-clip | 0 | 2 | 1,388,600 | 1,979,288 |
| | | Revised | 0 | 9 | 1,046,076 | 312,430 |
| Chromosome end | 39 | Soft-clip | 7 | 14 | 1,459,457 | 499,657 |
| | | Revised | 11 | 26 | 297,671 | 65,110 |
| Flanking unplaced contigs | 32 | Soft-clip | 0 | 0 | 20,722,546 | 33,566,957 |
| | | Revised | 2 | 2 | 8,873,435 | 32,955,618 |
| All | 123 | Soft-clip | 7 | 26 | 29,150,566 | 37,658,680 |
| | | Revised | 35 | 61 | 10,336,955 | 33,344,434 |

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

On this basis, GPatch includes a conservative extension step that runs after alignment-based placement. It is enabled by default and disabled with `--no_extension`. Contigs left out for low mapping quality or nesting are candidates. For each gap between adjacent placed contigs, and for each chromosome end, walks extend from both flanking contigs into the gap. At each step, an unplaced contig is added if:

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
- Junction code: GPatch commit b02e259 (default mode; soft-clip baseline = --legacy_junctions).
  Evaluation: compare_overlap_trimming.py v0.2.0; prototype gpatch_junction_prototype.py v0.3.0
  produced identical outputs to GPatch b02e259 on these data.
- Placement changes: GPatch commit 9008618 (unplaced report; nested check on chained extents) and
  commit e72d0cb (overlap extension). e72d0cb is after the v0.5.0 tag, so the version that
  includes extension still needs a release/version number. GPatch e72d0cb reproduces the
  prototype (overlap_extension_prototype.py v0.1.0 + gpatch_junction_prototype.py --propagated)
  exactly on these data (176 contigs placed; identical contigs.bed, patches.bed and junction
  decisions). The 5 kb minimum overlap was chosen after observing one spurious 585 bp overlap on
  this same simulation.
- The -m 1 / -m 0 junction errors exclude junctions between contigs from different benchmark
  chromosomes, so they understate the cost of lowering the threshold.
- The Table X benchmark (148 placed) predates the chained-extent nested check; rerunning it on the
  155-contig placement changes the totals only slightly (43.69 vs 43.71 Mb, all chromosomes) but
  the numbers in the text would need updating if the paper reports the final GPatch version.
- Code: simulation/T2T-HG002/simulate_contigs_from_alignments.sh (-M 100000),
  paf_to_simulated_contigs.py v0.2.0, get_hg002_annotations.sh, breakpoint_repeat_enrichment.py
  v0.1.0 (seed 20261109). Results: T2T-HG002.HG002.1.simulated.merge100000.breakpoint-enrichment.tsv
-->
