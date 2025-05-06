## Simulated Data Generation

### Simulate the data
```
awk 'BEGIN {LEN=0}; {if ($1~">") {if (NR > 1) print LEN; LEN=0} else {LEN+=length($0)}}; END {print LEN}' ../../data/HG002-full-0.14.1.pat.fa > HG002.pat.contig_lengths.txt
../simulate_contigs.py -f ../../data/GCA_009914755.4_T2T-CHM13v2.0_genomic.chroms.fa -l HG002.pat.contig_lengths.txt -p CHM13.pseudocontigs.HG002.1
```

### Simulated indels
```
SURVIVOR simSV CHM13.pseudocontigs.HG002.1.rc.contigs.fa SURVIVOR.params.HG002.1 0.01 0 CHM13.pseudocontigs.HG002.1.rc.contigs.SURVIVOR
```

### Add back unmutated contigs...
```
../fix_SURVIVOR_fasta.py -r CHM13.pseudocontigs.HG002.1.rc.contigs.fa -s CHM13.pseudocontigs.HG002.1.rc.contigs.SURVIVOR.fasta > CHM13.pseudocontigs.HG002.1.rc.contigs.SURVIVOR.with-missing.fasta
```

### Produce the contig-boundaries BED and concatenated FASTA
```
../build_SURVIVOR_contigs_bed.py -s CHM13.pseudocontigs.HG002.1.rc.contigs.SURVIVOR.with-missing.fasta -r CHM13.pseudocontigs.HG002.1.rc.contigs.bed -f CHM13.pseudocontigs.HG002.1.rc.contigs.SURVIVOR.with-missing.concat.fasta > CHM13.pseudocontigs.HG002.1.rc.contigs.SURVIVOR.with-missing.bed
```

