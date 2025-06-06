## Simulated Data Generation

### Simulate the data
```
grep -v ">" ../../data/v12_NA12878_giab_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta | awk '{print length($0)}' > NA12878.1.contig-lengths.txt
../simulate_contigs.py -f ../..//data/GenBank.GCA_009914755.4/GCA_009914755.4_T2T-CHM13v2.0_genomic.chroms.fa -l NA12878.1.contig-lengths.txt -r -p CHM13.pseudocontigs.NA12878.1.rc
```

### Simulated indels
```
SURVIVOR simSV CHM13.pseudocontigs.NA12878.1.rc.contigs.fa SURVIVOR.params.NA12878.1 0.01 0 CHM13.pseudocontigs.NA12878.1.rc.contigs.SURVIVOR
```

### Add back unmutated contigs...
```
../fix_SURVIVOR_fasta.py -r CHM13.pseudocontigs.NA12878.1.rc.contigs.fa -s CHM13.pseudocontigs.NA12878.1.rc.contigs.SURVIVOR.fasta > CHM13.pseudocontigs.NA12878.1.rc.contigs.SURVIVOR.with-missing.fasta
```

### Produce the contig-boundaries BED and concatenated FASTA
```
../build_SURVIVOR_contigs_bed.py -s CHM13.pseudocontigs.NA12878.1.rc.contigs.SURVIVOR.with-missing.fasta -r CHM13.pseudocontigs.NA12878.1.rc.contigs.bed -f CHM13.pseudocontigs.NA12878.1.rc.contigs.SURVIVOR.with-missing.concat.fasta > CHM13.pseudocontigs.NA12878.1.rc.contigs.SURVIVOR.with-missing.bed
```

