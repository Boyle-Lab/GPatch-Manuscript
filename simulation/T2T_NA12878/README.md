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
