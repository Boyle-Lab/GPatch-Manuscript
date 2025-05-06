# Simulated Contig Genome Patching with RagTag Patch
## T2T-CHM13 pseudocontigs based on HGSVC NA12878 contig size distribution
## The following steps were used for patching simulated CHM13 contigs using RagTag Patch:

## Dependencies
* Python >= v3.7
* RagTag (https://github.com/malonge/RagTag/)
* samtools (https://github.com/samtools/samtools)
* biopython (https://biopython.org/)
* pysam (https://github.com/pysam-developers/pysam)
* minimap2 (https://github.com/lh3/minimap2)

### RagTag patching of NA12878 simulated contigs, no indels
```
/usr/bin/time -v ./do_RagTag_Patch.sh ../CHM13.pseudocontigs.NA12878.1.rc.contigs.fa ../../../data/GCA_009914755.4_T2T-CHM13v2.0_genomic.chroms.fa T2T-CHM13 CHM13.NA12878.1 no_indel
wc -l ../CHM13.pseudocontigs.NA12878.1.rc.contigs.bed
wc -l no_indel/ragtag.patch.contigs.bed
../../compare_contig_positions.py -r ../CHM13.pseudocontigs.NA12878.1.rc.contigs.bed -p no_indel/ragtag.patch.contigs.bed > no_indel/ragtag.patch.contigs.mapping-stats.txt
../../compare_contig_ordering.py -r ../CHM13.pseudocontigs.NA12878.1.rc.contigs.bed -p no_indel/ragtag.patch.contigs.bed > no_indel/ragtag.patch.contigs.ordering-stats.txt
../../compare_fasta_by_editdist.py -r ../../../data/GCA_009914755.4_T2T-CHM13v2.0_genomic.chroms.fa -p no_indel/ragtag.patch.scf.fasta > no_indel/ragtag.patch.contigs.editdist-stats.txt
../../compile_mapping_stats.Rscript no_indel/ragtag.patch.contigs.mapping-stats.txt > no_indel/ragtag.patch.contigs.mapping-summary.txt
n50 no_indel/ragtag.patch.scf.fasta
```

### RagTag patching of NA12878 simulated contigs, with simulated indels from SURVIVOR
```
/usr/bin/time -v ./do_RagTag_Patch.SURVIVOR.sh ../CHM13.pseudocontigs.NA12878.1.rc.contigs.SURVIVOR.with-missing.fasta ../../../data/GCA_009914755.4_T2T-CHM13v2.0_genomic.chroms.fa T2T-CHM13.SURVIVOR CHM13.NA12878.1.SURVIVOR SURVIVOR
```

### RagTag Patch patching attempt for SURVIVOR dataset with nucmer alignment
/usr/bin/time -v ./do_RagTag_Patch.SURVIVOR.nucmer.sh ../CHM13.pseudocontigs.HG002.1.rc.contigs.SURVIVOR.with-missing.fasta ../../../data/GCA_009914755.4_T2T-CHM13v2.0_genomic.chroms.fa T2T-CHM13.SURVIVOR CHM13.HG002.1.SURVIVOR SURVIVOR.nucmer
