# Alignment-guided contig simulation (NA12878 benchmark)

This is a simulated contig-level assembly with known order and orientation, for benchmarking
reference-guided patching against T2T-CHM13. It uses the same method as `../T2T-HG002/`.

* **Model assembly (NA12878.2):** HGSVC v12, h1 contigs:
  `manuscript/data/v12_NA12878_giab_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta`
* **Benchmark genome (NA12878.1):** maternal chromosome-scale assembly, `manuscript/data/NA12878.1.chroms.fasta`.
  It has 54 records:
  * 13 chromosomes are complete (chr3, 4, 8, 10–13, 15–19, X);
  * the other 10 are split into 2–14 contigs.
* **Patching reference:** T2T-CHM13

## Benchmark contig IDs

Every benchmark record is named only by its chromosome (`>chr1 mat-0000003`, `>chr1 mat-0000006`, …).
minimap2 needs unique target names, so `rename_benchmark_contigs.sh` writes
`NA12878.1.chroms.renamed.fa`, with IDs `<chrom>.<descriptor>` (e.g. `chr1.mat-0000003`), and indexes it.
The sequences are unchanged. Simulated contigs are therefore named `<chrom>.<descriptor>_<rank>`,
e.g. `chr1.mat-0000033_12`. Ranks are 0-based within each benchmark contig, and each contig's BED
coordinates are in its benchmark contig.

## Run

```
./rename_benchmark_contigs.sh
./simulate_contigs_from_alignments.sh -t 32
./simulate_contigs_from_alignments.sh -a NA12878.1.NA12878.2.simulated.paf \
    -p NA12878.1.NA12878.2.simulated.merge100000 -M 100000
```

`simulate_contigs_from_alignments.sh` and `paf_to_simulated_contigs.py` are copies of the T2T-HG002
scripts. Only the defaults differ: model, benchmark and output prefix. The output formats are
described in `../T2T-HG002/README.md`. Junctions are only defined between contigs on the same
benchmark contig.

## Results

| `-M` | contigs | N50 | overlaps / gaps | benchmark covered |
|---|---|---|---|---|
| off | 10,606 | 6.38 Mb | 50% / 50% | 95.83% |
| 100 kb | 7,400 | 18.30 Mb | 58% / 42% | 96.08% |

* The model assembly has an N50 of 18.32 Mb and 10,406 contigs. Its median contig length is short,
  so most simulated contigs are short: the median is 16 kb with `-M 100000`.
* 48 of the 54 benchmark contigs receive simulated contigs. Six small ones (13–147 kb) receive none:
  * chr1.mat-0000047, chr1.mat-0000126;
  * chr2.mat-0000119, chr2.mat-0000125, chr2.mat-0000127;
  * chr7.mat-0000124.
* The blen-weighted divergence between model and benchmark is 0.07%.

## Model-assembly contigs excluded from the simulation (`-M 100000`)

This repeats the T2T-HG002 analysis. Results are in `excluded_contigs/`: per-contig table
`NA12878.2.excluded_contigs.tsv` and the alignments behind it.

Of the 10,406 model contigs (2,983 Mb), 6,724 (2,924 Mb) contributed to the simulation. The other
3,682 (59.0 Mb; 1.98% of the model assembly; median 14.3 kb, maximum 181 kb) did not.

* **Contained: 3,667 contigs, 58.8 Mb.** Each aligns only inside a retained contig's benchmark
  interval, so no benchmark sequence is lost.
  * Containing interval: 2,367 (35.2 Mb) lie within a single aligned block of a retained contig.
    1,300 (23.6 Mb) lie only within a span bridged by `-M 100000` merging.
  * Satellite content: 2,789 (46.8 Mb) are satellite-like (21-mer diversity < 0.5 or ≥ 20%
    HSat2/3 `CATTC` motif).
  * Location: they pile up in pericentromeric heterochromatin and acrocentric arms (e.g. chr9q12,
    chr1q12, chr7 centromere, chr13/14/22 short arms). 84 one-Mb windows hold 83% of their
    alignment blocks.
  * MAPQ: best MAPQ < 10 for 1,716; ≥ 50 for 731.
  * Haplotype: aligning them to both parental NA12878 assemblies (asm5) gives mixed calls. Of the
    878 non-satellite contigs, 232 match maternal better, 196 paternal, 174 both equally, and 276
    neither at ≥ 99%. Other-haplotype contigs contribute, but most contained contigs are repeat
    fragments.
* **No alignment: 4 contigs (62 kb). Primary blocks all < 1 kb: 11 (0.22 Mb).** All 15 are human.
  None align to EBV (NC_007605.1, asm20), so BLAST was not needed. This differs from HG002.1, which
  had 46 EBV and 44 contaminant contigs. Using CHM13, the paternal NA12878.2 assembly and
  `asm20` against NA12878.1:
  * 5 (89 kb) are paternal-specific. They are absent from or highly divergent in the maternal
    benchmark but match NA12878.2 at ≥ 99.9% identity. Two of them (chr5, chr18) are non-reference
    insertions that barely align to CHM13.
  * 9 (121 kb) are centromeric/acrocentric satellites, where CHM13 positions place them:
    * chr4 42 bp AT-rich, chr16 GC-rich, chr8 and chr9 α-satellite;
    * chr14/chr21 short arms;
    * HSat3 absent from both NA12878 assemblies.
    They align to the benchmark only at 72–94% identity, or not at all.
  * 1 (66 kb) is unique sequence (CHM13 chr7:60.36 Mb, MAPQ 60) that falls in the NA12878.1 chr7
    pericentromeric gap between chr7.mat-0000001 and chr7.mat-0000036.

  The benchmark is thus missing a little real NA12878 sequence: one 66 kb pericentromeric gap
  segment, plus satellite arrays that it does not represent.

**The h1 model is not a single parental haplotype.** HGSVC v12 h1 is phased by Strand-seq without
parental assignment. In a 31-contig sample (1–5 Mb each, 73 Mb), 16 contigs match the maternal
benchmark better and 15 the paternal assembly. Most match one parent with gap-compressed
divergence ≤ 0.0002 and the other at ~0.001, and the match is consistent within a chromosome:
* maternal: chr1, 5, 8, 12, 16, 17, 18, X;
* paternal: chr2, 4, 6, 7, 10, 11, 14, 15, 19, 22.

For about half the chromosomes, the simulated breakpoints therefore come from paternal contigs
aligned to the maternal benchmark. Alignment breaks at heterozygous SVs may then add breakpoints
that a maternal assembly would not have. The simulated sequence itself always comes from the
maternal benchmark.

## Simulated contigs left unplaced by GPatch

This repeats the T2T-HG002 analysis for GPatch 0.6.0 (`-m 10`, default junction handling and
overlap extension) on the `-M 100000` simulation. NA12878.1 has no repeat annotation, so repeat
context comes from two sources:
* CHM13 CenSat (v2.1, excluding `ct` transition entries) and SD annotations at each contig's
  primary CHM13 alignment, in `annotation/raw/`;
* the contig sequence itself: 21-mer diversity < 0.5 or ≥ 20% HSat2/3 `CATTC`.

```
P=NA12878.1.NA12878.2.simulated.merge100000
./characterize_unplaced_contigs.py -s $P.contigs.bed -j $P.junctions.tsv -a $P.contigs.fa \
    -g GPatch/$P.contigs.bed -u GPatch/$P.unplaced.tsv -e GPatch/$P.extension.tsv -b GPatch/$P.bam \
    --censat annotation/raw/chm13v2.0_censat_v2.1.bed --sd annotation/raw/chm13v2.0_SD.bed \
    -o GPatch/$P.unplaced-characterization.tsv > GPatch/$P.unplaced-characterization.summary.txt
./analyze_extension_limits.py -s $P.contigs.bed -j $P.junctions.tsv -g GPatch/$P.contigs.bed \
    -u GPatch/$P.unplaced.tsv -e GPatch/$P.extension.tsv -b GPatch/$P.bam > GPatch/$P.extension-limits.txt
```

GPatch placed 2,704 of the 7,400 simulated contigs (2,845 Mb), 296 of them by overlap extension.
It left 4,696 unplaced (93.1 Mb, 3.2% of the simulated sequence). For HG002 the figures were 184
contigs and 26.3 Mb (0.9%).

| Reason | Contigs | Mb | Median length | Repeat context |
|---|---|---|---|---|
| Low MAPQ (< 10) | 4,516 | 89.58 | 15.1 kb | 4,515 CenSat/SD at the CHM13 alignment; 4,077 satellite-like sequence |
| Nested | 174 | 3.52 | 13.2 kb | 172 CenSat/SD; 131 satellite-like |
| Unmapped | 6 | 0.01 | 1.3 kb | 3 satellite-like |

* **Low MAPQ.** Every contig but one lies in CHM13 CenSat (4,342) or SD (491) sequence.
  * MAPQ is 0 for 350 and 1–9 for 4,166.
  * Alignments are complete (median aligned fraction 1.00) but not unique.
  * Location relative to placed benchmark neighbours:
    * 3,207 align between their CHM13 positions (± 2 Mb);
    * 698 align elsewhere on the right chromosome;
    * 326 align to the wrong chromosome;
    * 285 have no placed neighbour.

  Most are on chr1 (1,020, mainly the 1q12 HSat2 array), chr9 (494), chr16 (340), chr15, chr3
  and chr19. 286 (4.8 Mb) belong to the 13 benchmark contigs with no placed contig at all,
  which are small chr1, chr5, chr6 and chr20 fragments.
* **Nested.** 134 are nested by their aligned coordinates, mostly acrocentric and pericentromeric
  homology (chr21, chr22, chr15, chr14). 40 are nested only by the containing contig's chained
  extent, i.e. its residual unaligned ends.
* **Unmapped.** 6 tiny contigs (1–2 kb).

**Why overlap extension did not recover more.** 2,892 unplaced contigs (68.1 Mb) overlap a true
neighbour by ≥ 5 kb, but extension can only walk inward from placed contigs.
* Within runs of unplaced contigs, a third of junctions (1,398 of 4,193) are true gaps and 743
  are overlaps < 5 kb, so the chains break.
* Only 376 unplaced contigs (11.6 Mb) are reachable from a placed contig through overlaps ≥ 5 kb.
  4,030 (76.8 Mb) are cut off by gaps or short overlaps.
* Unplaced contigs form 2,644 overlap islands: 1,786 singletons, median span 16 kb, 110 islands
  ≥ 100 kb (19.4 Mb).
* Of 4,240 walk ends:
  * 577 reached the flanking contig.
  * 3,363 stopped because the next true junction is a gap or short overlap (2,597 of them next to
    an already placed neighbour, i.e. a gap with nothing to insert).
  * 59 stopped at ambiguous candidates.
  * 72 stopped because GPatch's alignment rules excluded a truly overlapping neighbour: 46 outside
    the gap region, 17 on another chromosome, 9 in the opposite orientation.
  * The rest stopped at a benchmark-contig end (40), a walk later rejected as claimed by several
    walks (11), an unmapped neighbour (1), or a true neighbour already placed elsewhere despite an
    overlap ≥ 5 kb (117).

Compared with HG002, the unplaced sequence is 3.5× larger and in far more, smaller pieces. The
HGSVC v12 CCS assembly breaks centromeric and pericentromeric satellite arrays into thousands of
~15 kb contigs separated by gaps. Alignment to CHM13 cannot order these contigs (MAPQ < 10), and
overlap walks from placed contigs reach only the ends of each run.

## Benchmarking against the split benchmark

The relative order and orientation of benchmark contigs from the same chromosome are unknown. Use the
`*_split.py` comparison scripts in `..` (described in `../README.md`). They compare within one benchmark
contig at a time and let each benchmark contig be reversed relative to CHM13.

```
P=NA12878.1.NA12878.2.simulated.merge100000
S=..
$S/compare_contig_ordering_split.py  -r $P.contigs.bed -p GPatch/<out>.contigs.bed -d GPatch/<out>.ordering-details.tsv > GPatch/<out>.contigs.ordering-stats.txt
$S/compare_contig_positions_split.py -r $P.contigs.bed -p GPatch/<out>.contigs.bed -b GPatch/<out>.bam > GPatch/<out>.contigs.mapping-stats.txt
$S/compare_fasta_by_editdist_split.py -r NA12878.1.chroms.renamed.fa -p GPatch/<out>.patched.fasta \
    -s $P.contigs.bed -g GPatch/<out>.contigs.bed -d GPatch/<out>.editdist-details.tsv > GPatch/<out>.contigs.editdist-stats.txt
```

The comparison FASTA must be the renamed benchmark. For RagTag-derived BEDs, which have a 9th scaffold
column, pass `--fasta-id-col 9` to the edit-distance script if the patched FASTA is keyed by scaffold.
