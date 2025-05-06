# Simulated Contig Genome Patching with GPatch
## T2T-CHM13 pseudocontigs based on HGSVC HG002 contig size distribution
## The following steps were used for patching simulated CHM13 contigs using GPatch:

## Dependencies
* Python >= v3.7
* GPatch (https://github.com/adadiehl/GPatch/)
* samtools (https://github.com/samtools/samtools)
* biopython (https://biopython.org/)
* pysam (https://github.com/pysam-developers/pysam)
* minimap2 (https://github.com/lh3/minimap2)

### GPatch patching of HG002 simulated contigs, no indels
```
/usr/bin/time -v ./patch_genome.sh ../CHM13.pseudocontigs.HG002.1.rc.contigs.fa ../../../data/GCA_009914755.4_T2T-CHM13v2.0_genomic.chroms.fa T2T-CHM13 CHM13.HG002.1
wc -l ../CHM13.pseudocontigs.HG002.1.rc.contigs.bed
wc -l CHM13.HG002.1.contigs.bed
../../compare_contig_positions.py -r ../CHM13.pseudocontigs.HG002.1.rc.contigs.bed -p CHM13.HG002.1.contigs.bed -b CHM13.HG002.1.bam > CHM13.HG002.1.contigs.mapping-stats.txt
../../compare_contig_ordering.py -r ../CHM13.pseudocontigs.HG002.1.rc.contigs.bed -p CHM13.HG002.1.contigs.bed > CHM13.HG002.1.contigs.ordering-stats.txt
../../compare_fasta_by_editdist.py -r ../../../data/GCA_009914755.4_T2T-CHM13v2.0_genomic.chroms.fa -p CHM13.HG002.1.patched.fasta > CHM13.HG002.1.contigs.editdist-stats.txt
../../compile_mapping_stats.Rscript CHM13.HG002.1.contigs.mapping-stats.txt > CHM13.HG002.1.contigs.mapping-summary.txt
n50 CHM13.HG002.1.patched.fasta 
```

### GPatch patching of HG002 simulated contigs, with simulated indels from SURVIVOR
```
/usr/bin/time -v ./patch_genome.SURVIVOR.sh ../CHM13.pseudocontigs.HG002.1.rc.contigs.SURVIVOR.with-missing.fasta ../../../data/GCA_009914755.4_T2T-CHM13v2.0_genomic.chroms.fa T2T-CHM13.SURVIVOR CHM13.HG002.1.SURVIVOR ../CHM13.pseudocontigs.HG002.1.rc.contigs.SURVIVOR.with-missing.concat.fasta
wc -l ../CHM13.pseudocontigs.HG002.1.rc.contigs.SURVIVOR.with-missing.bed
wc -l CHM13.HG002.1.SURVIVOR.contigs.bed
../../compare_contig_positions.py -r ../CHM13.pseudocontigs.HG002.1.rc.contigs.SURVIVOR.with-missing.bed -p CHM13.HG002.1.SURVIVOR.contigs.bed -b CHM13.HG002.1.SURVIVOR.bam > CHM13.HG002.1.SURVIVOR.contigs.mapping-stats.txt
../../compare_contig_ordering.py -r ../CHM13.pseudocontigs.HG002.1.rc.contigs.SURVIVOR.with-missing.bed -p CHM13.HG002.1.SURVIVOR.contigs.bed > CHM13.HG002.1.SURVIVOR.contigs.ordering-stats.txt
../../compare_fasta_by_editdist.py -r ../CHM13.pseudocontigs.HG002.1.rc.contigs.SURVIVOR.with-missing.concat.fasta -p CHM13.HG002.1.SURVIVOR.patched.fasta > CHM13.HG002.1.SURVIVOR.contigs.editdist-stats.txt
../../compile_mapping_stats.Rscript CHM13.HG002.1.contigs.mapping-stats.txt > CHM13.HG002.1.contigs.mapping-summary.txt
n50 ../CHM13.pseudocontigs.HG002.1.rc.contigs.SURVIVOR.with-missing.concat.fasta 
n50 CHM13.HG002.1.SURVIVOR.patched.fasta
```
