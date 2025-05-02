# HGSVC NA12878 patching and analysis with GPatch
Analysis steps for producing GPatch NA12878 from HGSVC NA12878.

## Dependencies
* Python >= v3.7
* GPatch https://github.com/adadiehl/GPatch/
* samtools (https://github.com/samtools/samtools)
* biopython (https://biopython.org/)
* pysam (https://github.com/pysam-developers/pysam)
* minimap2 (https://github.com/lh3/minimap2
* N50 (https://anaconda.org/bioconda/n50)


### Alignment, analysis, and dot-plot generation for HGSVC NA12878
/usr/bin/time -v ./patch_genome.sh ../../data/v12_NA12878_giab_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta ../../data/GCA_009914755.4_T2T-CHM13v2.0_genomic.chroms.fa T2T-CHM13 NA12878.2 "" "-d -t -m 10"

### Generate dot-plots to T2T-NA12878
./do_target_dotplots.sh NA12878.2.cbreak_2.patched.fasta ../../data/NA12878.1.chroms.complete.fasta NA12878_MAT NA12878.2.cbreak_2
./do_target_dotplots.sh NA12878.2.cbreak.patched.fasta ../../data/NA12878.1.chroms.complete.fasta NA12878_MAT NA12878.2.cbreak
./do_target_dotplots.sh NA12878.2.patched.fasta ../../data/NA12878.1.chroms.complete.fasta NA12878_MAT NA12878.2

### Calculate N50s
n50 ../../data/v12_NA12878_giab_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta
n50 ../../data/NA12878.1.chroms.complete.fasta
n50 ../../data/NA12878.1.chroms.fasta
n50 ../../data/GCA_009914755.4_T2T-CHM13v2.0_genomic.chroms.fa
n50 NA12878.2.cbreak_2.patched.fasta 
n50 NA12878.2.cbreak.patched.fasta 
n50 NA12878.2.patched.fasta

### Generate data for phenogram
awk 'BEGIN {printf "chr\tpos\tendpos\tpheno\tposcolor\n"}; {gsub(/chr/, "", $1); printf "%s\t%d\t%d\tCHM13\t2\n", $1, $2, $3}' HG002.1.patches.bed | grep -v "^Y" | grep -v "^M" > HG002.1.patches.phenogram.txt
