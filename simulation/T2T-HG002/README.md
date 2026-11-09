# Alignment-guided contig simulation (T2T-HG002 benchmark)

Simulated contig-level assembly with known order and orientation, for benchmarking
reference-guided patching (GPatch) against T2T-CHM13.

* **Model assembly:** HG002.1 (HGSVC), paternal contigs: `manuscript/data/HG002-full-0.14.1.pat.fa`
* **Benchmark genome:** T2T-HG002 paternal chromosomes: `manuscript/data/HG002_paternal.chroms.fa`
* **Patching reference:** T2T-CHM13

## Strategy

Model contigs are aligned to the benchmark genome with `minimap2 -x asm5 -c --secondary=no`.
Every primary alignment block becomes one simulated contig spanning the block's
target interval. Its sequence comes from the benchmark genome and is reverse-complemented
for `-` strand blocks. The breakpoints are where real contigs break, or stop aligning,
against the benchmark genome. That makes them non-random and enriched in repeats/SDs
without needing repeat annotations for T2T-HG002. Overlaps and unfilled gaps between
adjacent contigs come straight from the alignment boundaries.

Filters (defaults): primary blocks only (`tp:A:P`), MAPQ >= 0, benchmark span >= 1 kb.
Blocks whose benchmark interval lies entirely inside another block are dropped (`-k` keeps them).
Contigs are named `<chrom>_<rank>`, 0-based in benchmark coordinate order.

**Merging split model contigs (`-M N`, off by default).** By default every alignment block is
its own contig, so a model contig that aligns in several pieces (around SVs or divergent
sequence) is split even though it is continuous in HG002.1. With `-M N`, consecutive blocks
from the same model contig are joined if they lie on the same chromosome and strand, are
colinear in both genomes, and are <= N bp apart (gap or overlap) in both. The merged contig
spans the whole benchmark interval, including the sequence between the pieces. Merging happens
before the length and containment filters, so small pieces can be bridged, and a different
contig lying inside the merged span is then dropped as contained. Pieces that change strand or
order (inversions, rearrangements) are never merged.

| `-M` | contigs | N50 | overlaps / gaps | same-contig junctions left |
|---|---|---|---|---|
| off | 430 | 55.6 Mb | 47% / 53% | 140 |
| 10 kb | 375 | 76.5 Mb | 52% / 48% | 87 |
| 100 kb | 332 | 76.5 Mb | 57% / 43% | 53 |
| 1 Mb | 300 | 81.5 Mb | 61% / 39% | 34 |

(HG002.1 paternal model assembly N50: 85.0 Mb.)

## Run

```
./simulate_contigs_from_alignments.sh -t 32
# merged variant, reusing the alignment:
./simulate_contigs_from_alignments.sh -a T2T-HG002.HG002.1.simulated.paf \
    -p T2T-HG002.HG002.1.simulated.merge100000 -M 100000
```

Options: `-q` model FASTA, `-r` benchmark FASTA, `-p` output prefix, `-m` min MAPQ,
`-l` min length, `-M` merge distance, `-a` existing PAF to reuse, `-k` keep contained blocks, `-f` force re-alignment (an existing PAF is
reused otherwise). Tools come from `GPatch_env` (minimap2 2.28, Biopython 1.84). Commands,
versions and git commit are recorded in `<prefix>.log`.

## Breakpoint repeat/SD enrichment

```
./get_hg002_annotations.sh        # HG002 Q100 v1.1 paternal SD, CenSat, RepeatMasker, TRF, mnr10 -> annotation/raw/
R=annotation/raw
python3 breakpoint_repeat_enrichment.py \
    -b T2T-HG002.HG002.1.simulated.merge100000.contigs.bed -p T2T-HG002.HG002.1.simulated.paf \
    --sd $R/hg002v1.1.pat.SDs.013025.bed --censat $R/hg002v1.1.cenSatv2.0.PAT.bed \
    --rmsk $R/hg002v1.1.paternal.rmsk.bed --trf $R/hg002v1.1.simpleRepeat.bed --mnr $R/v1.1.mnr10.pat.bed \
    -o T2T-HG002.HG002.1.simulated.merge100000.breakpoint-enrichment.tsv
```
The test uses contig endpoints (all, true HG002.1 termini, internal alignment breaks) and
junctions. The null is 1,000 within-chromosome random placements (seed 20261109). The RepeatMasker
bigBed (joined format) is converted with pyBigWig from `HiCExplorer_env`, and only aligned blocks
are used. With `-M 100000`, 94.6% of breakpoints fall in CenSat or SDs (12.7% of the genome):
α-satellite 12.4×, rDNA 54×, SD ≥98% 4.5×, while interspersed repeats are depleted.

## Outputs (prefix `T2T-HG002.HG002.1.simulated`)

| File | Contents |
|---|---|
| `.paf` | raw minimap2 alignments |
| `.contigs.fa` | simulated contigs |
| `.contigs.bed` | contig locations in T2T-HG002 (see below) |
| `.junctions.tsv` | each adjacent pair: signed distance (<0 overlap, 0 abut, >0 gap) |
| `.summary.txt` | filter counts, contig/N50 stats, overlap & gap rates and sizes, coverage |

`contigs.bed` has 8 columns and no header:
`chrom start end name length strand mapq attributes`. The attributes are
`src` (model contig), `qlen`, `qs`, `qe`, `AS`, `NM`, `matches`, `blen`, `de` and `nblk`
(the number of alignment blocks merged into the contig). For merged contigs, `AS`, `NM`,
`matches` and `blen` are summed over the blocks, `mapq` is the lowest block MAPQ, and `de`
is the mean weighted by `blen`.
Column 5 is the length, not a BED score. The file stops at 8 columns to stay compatible with
`../compare_contig_ordering.py`, which treats >8 columns as RagTag scaffolds, and with
`../compare_contig_positions.py`.
