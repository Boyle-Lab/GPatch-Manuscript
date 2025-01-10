# Simulated Contig Genome Patching with GPatch
## T2T-CHM13 pseudocontigs based on HGSVC NA12878 contig size distribution
## The following steps were used for patching simulated CHM13 contigs using GPatch:

## Dependencies
* Python >= v3.7
* GPatch (https://github.com/adadiehl/GPatch/)
* samtools (https://github.com/samtools/samtools)
* biopython (https://biopython.org/)
* pysam (https://github.com/pysam-developers/pysam)
* minimap2 (https://github.com/lh3/minimap2)

./patch_genome.sh ../CHM13.pseudocontigs.NA12878.1.contigs.fa ../../../data/GenBank.GCA_009914755.4/GCA_009914755.4_T2T-CHM13v2.0_genomic.chroms.fa T2T-CHM13 CHM13.NA12878.1 ../../../data/whitelist.T2T_CHM13.bed
