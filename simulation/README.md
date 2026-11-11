# Simulation Analysis

This folder contains scripts used for the simulated data analyses performed in NA12878 and HG002 subfolders.

### Compiling overall stats for no-indel set
```
paste NA12878/GPatch/CHM13.NA12878.1.contigs.ordering-stats.txt <(awk '{for (i=2; i<NF; i++) printf $i "\t"; print $NF}' NA12878/GPatch/CHM13.NA12878.1.contigs.editdist-stats.txt) | column -s $'\t' -t | awk -v SAMPLE="NA12878" -v METHOD="GPatch" '{if (NR == 1) {printf "Sample\tMethod\t%s\tpct_assembly_localized\n", $0} else {printf "%s\t%s\t%s\t-1\n", SAMPLE, METHOD, $0}}' > simulation_stats.combined.txt

paste NA12878/RagTag/no_indel/*.ordering-stats.txt <(awk '{for (i=2; i<NF; i++) printf $i "\t"; print $NF}' NA12878/RagTag/no_indel/*.editdist-stats.txt) | column -s $'\t' -t | awk -v SAMPLE="NA12878" -v METHOD="RagTag" '{if (NR > 1) printf "%s\t%s\t%s\t-1\n", SAMPLE, METHOD, $0}' >> simulation_stats.combined.txt

paste HG002/GPatch/CHM13.HG002.1.contigs.ordering-stats.txt <(awk '{for (i=2; i<NF; i++) printf $i "\t"; print $NF}' HG002/GPatch/CHM13.HG002.1.contigs.editdist-stats.txt) | column -s $'\t' -t | awk -v SAMPLE="HG002" -v METHOD="GPatch" '{if (NR > 1) printf "%s\t%s\t%s\t-1\n", SAMPLE, METHOD, $0}' >> simulation_stats.combined.txt

paste HG002/RagTag/no_indel/*.ordering-stats.txt <(awk '{for (i=2; i<NF; i++) printf $i "\t"; pri
nt $NF}' HG002/RagTag/no_indel/*.editdist-stats.txt) | column -s $'\t' -t | awk -v SAMPLE="HG002" -v METHOD="RagTag" '{if (NR > 1) printf "%s\t%s\t%s\t-1\n", SAMPLE, METHOD, $0}' >> simulation_stats.combined.txt
```

### Compiling overall stats for	SURVIVOR set
```
paste NA12878/GPatch/CHM13.NA12878.1.SURVIVOR.contigs.ordering-stats.txt <(awk '{for (i=2; i<NF; i++) printf $i "\t"; print $NF}' NA12878/GPatch/CHM13.NA12878.1.SURVIVOR.contigs.editdist-stats.txt) | column -s $'\t' -t | awk -v SAMPLE="NA12878" -v METHOD="GPatch" '{if (NR == 1) {printf "Sample\tMethod\t%s\tpct_assembly_localized\n", $0} else {printf "%s\t%s\t%s\t-1\n", SAMPLE, METHOD, $0}}' > simulation_stats.SURVIVOR.combined.txt

paste HG002/GPatch/CHM13.HG002.1.SURVIVOR.contigs.ordering-stats.txt <(awk '{for (i=2; i<NF; i++) printf $i "\t"; print $NF}' HG002/GPatch/CHM13.HG002.1.SURVIVOR.contigs.editdist-stats.txt) | column -s $'\t' -t | awk -v SAMPLE="HG002" -v METHOD="GPatch" '{if (NR > 1) printf "%s\t%s\t%s\t-1\n", SAMPLE, METHOD, $0}' >> simulation_stats.SURVIVOR.combined.txt
```

### Produce bar plots for figure.
```
./do_barplots.Rscript simulation_stats.combined.txt "simulation_stats.combined"
```

### Data for Tables S2 and S3
Genome lengths in Tables S2 and S3 can be found in {NA12878,HG002}/{GPatch,RagTag}/\*.stats

Counts needed to calculate patched genome contig coverage %, and % assembly localized (Table S3) are also found in {NA12878,HG002}/{GPatch,RagTag}/\*.stats
Remaining rows, with the exception of runtime and max memory, are from simulation_stats.*.combined.txt

### Comparisons against fragmented benchmarks (`*_split.py`)
`compare_contig_ordering_split.py`, `compare_contig_positions_split.py` and
`compare_fasta_by_editdist_split.py` are versions of the three comparison scripts for benchmark
assemblies whose chromosomes are split into several contigs (NA12878.1; see `T2T_NA12878/`). They
share `split_benchmark.py`. The originals are unchanged.

* The simulated-contig BED must use unique benchmark contig IDs `<chrom>.<descriptor>`
  (`--sep` changes the separator). Simulated contigs are named `<benchmark contig>_<rank>`.
* Comparisons are made within one benchmark contig at a time, because the order and orientation
  of benchmark contigs from the same chromosome are unknown.
* Each benchmark contig gets an orientation from the length-weighted majority of its placed
  contigs: `-` when they mostly lie on the opposite strand to their simulated strand. Expected
  strands and order are reversed for `-` contigs, and patched sequence is reverse-complemented
  before the sequence comparison.
* The *target* sequence of a benchmark contig is the patched sequence holding most of its placed
  contigs.
* **Ordering:** the patched vector runs from the benchmark contig's first to its last contig on
  the target. Contigs from other benchmark contigs placed in between are kept and count as
  errors. Contigs outside that span are not compared, unlike the whole-chromosome original.
* **Positions:** the original's exact-coordinate test becomes a relative one. A contig is
  correctly mapped if it is on the target, the target is named after the benchmark chromosome,
  and its start and end offsets (patched − benchmark, after orientation) equal the benchmark
  contig's most common offset. Order is checked against the previous *placed* contig of the same
  benchmark contig.
* **Sequence:** a complete chromosome (the only benchmark contig for its chromosome) is compared
  over the whole patched chromosome, as in the original. Fragmented ones are compared over the
  span of their contigs on the target (`--scope span` uses spans for all).
* `compile_mapping_stats_split.Rscript` summarizes `compare_contig_positions_split.py` output. It is
  the same as `compile_mapping_stats.Rscript`, except that clustering accuracy compares only the
  chromosome part of `ref_chrom` (before `.`) with `mapped_chrom`.
* Output columns are the same as the originals, with one row per benchmark contig. `-d` writes
  per-benchmark-contig orientation, target, strand agreement and scope.

On the T2T-HG002 simulation (unfragmented benchmark), the split versions reproduce the originals'
edit distances, orientation calls and sequence distances exactly. Pair accuracy differs on
chromosomes with foreign contigs outside the span, and position/order calls differ by
definition. Flipping a patched chromosome (coordinates, strands and sequence) leaves all
results unchanged.
