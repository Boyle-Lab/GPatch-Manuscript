# HPRC HG002 patching and analysis with GPatch
Analysis steps for producing GPatch HG002 from HPRC HG002.

## Dependencies
* Python >= v3.7
* GPatch https://github.com/adadiehl/GPatch/
* samtools (https://github.com/samtools/samtools)
* biopython (https://biopython.org/)
* pysam (https://github.com/pysam-developers/pysam)
* minimap2 (https://github.com/lh3/minimap2
* N50 (https://anaconda.org/bioconda/n50)


### Alignment, analysis, and dot-plot generation for HPRC HG002
/usr/bin/time -v ./patch_genome.sh ../../data/HG002-full-0.14.1.pat.fa ../../data/GCA_009914755.4_T2T-CHM13v2.0_genomic.chroms.fa T2T-CHM13 HG002.1 "" "-d -t -m 10"

### Generate dot-plots to T2T-NA12878
./do_target_dotplots.sh HG002.1.cbreak.patched.fasta ../../data/HG002_paternal.chroms.fa T2T-HG002_PAT HG002.1.cbreak
./do_target_dotplots.sh HG002.1.patched.fasta ../../data/HG002_paternal.chroms.fa T2T-HG002_PAT HG002.1

### Calculate N50s
n50 ../../data/HG002-full-0.14.1.pat.fa 
n50 ../../data/HG002_paternal.chroms.fa 
n50 HG002.1.cbreak.patched.fasta 
n50 HG002.1.patched.fasta
n50 ../../data/GCA_009914755.4_T2T-CHM13v2.0_genomic.chroms.fa

### Phenogram data table
awk 'BEGIN {printf "chr\tpos\tendpos\tpheno\tposcolor\n"}; {gsub(/chr/, "", $1); printf "%s\t%d\t%d\tCHM13\t2\n", $1, $2, $3}' HG002.1.patches.bed | grep -v "^Y" | grep -v "^M" > HG002.1.patches.phenogram.txt
