# Biological Data Analysis

This folder contains scripts used for the biological data analyses performed in NA12878 and HG002 subfolders.

biological_patching_stats.combined.txt was compiled manually from data in NA12878 and HG002 subdirectories.

### Compile counts/pct for biological_patching_stats.combined.txt (Presented in Figure 3A-C)

Pct contigs localized = placed contigs / total contigs
```
# NA12878 placed contigs
wc -l NA12878.2.cbreak_2.contigs.bed
# Total NA12878 contigs
grep ">" ../../data/v12_NA12878_giab_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta | wc -l

# HG002 placed contigs
wc -l HG002/HG002.1.contigs.bed
# Total HG002 contigs
grep ">" ../../data/HG002-full-0.14.1.pat.fa | wc -l
```

Pct assembly localized
```
# NA12878 = post-break_2 contig base count / HGSVC NA12878 total nucleotide count
grep "post-break_2 contig base count" NA12878/NA12878.2.cbreak_2.contigs.bed
grep "unpatched assembly genome size" NA12878/NA12878.2.cbreak_2.contigs.bed
# HG002 = pre-break contig base count / HPRC HG002 total nucleotide count
grep "pre-break contig base count" HG002/HG002.1.contigs.bed
grep "unpatched	assembly genome	size" HG002/HG002.1.contigs.bed
```

Pct contig coverage
```
# NA12878 = post-break_2 contig base count / GPatc NA12878 total nucleotide count
grep "post-break_2 contig base count" NA12878/NA12878.2.cbreak_2.contigs.bed
grep "post-break_2 total patched genome size" NA12878/NA12878.2.cbreak_2.contigs.bed
# HG002 = pre-break contig base count / GPatch HG002 total nucleotide count
grep "pre-break contig base count" HG002/HG002.1.contigs.bed
grep "post-break_2 total patched genome size" HG002/HG002.1.contigs.bed
```

### Bar plots for figures 3A-C
./do_barplots.Rscript biological_patching_stats.combined.txt "biological_patching_stats.combined"
