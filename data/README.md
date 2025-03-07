# Public Data Retrieval and Preparation Steps

## T2T-CHM13 Genome
Downloaded the CHM13 Genome from NCBI, in zip format, from https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/
```
unzip ncbi_dataset.zip
```

Convert fasta IDs to chromosome names
```
awk '{if ($0 ~ ">") {printf ">chr%s\n", $NF} else {print $0}}' T2T-chm13v2.0/GenBank.GCA_009914755.4/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna | sed  's/chrgenome/chrM/' > GCA_009914755.4_T2T-CHM13v2.0_genomic.chroms.fa
rm -rf  T2T-chm13v2.0/
```
Generate chrom.sizes for T2T-CHM13
```
awk '{if ($1 ~ ">") {if (chrom!="") {printf "%s\t%d\n", chrom, len} chrom=$1; len=0} else {len += length($1)}}' GCA_009914755.4_T2T-CHM13v2.0_genomic.chroms.fa | sed 's/>//' > GCA_009914755.4_T2T-CHM13v2.0_genomic.chroms.chrom.sizes
```

## T2T-CHM13 Whitelist
First, download the T2T blacklist regions from Google Drive, at https://drive.google.com/drive/folders/1sF9m8Y3eZouTZ3IEEywjs2kfHOWFBSJT
```
bedtools sort -i T2T.excluderanges.bed.gz -g GCA_009914755.4_T2T-CHM13v2.0_genomic.chroms.chrom.sizes > T2T_CHM13.excluderanges.sorted.bed
bedtools complement -i T2T_CHM13.excluderanges.sorted.bed -g GCA_009914755.4_T2T-CHM13v2.0_genomic.chroms.chrom.sizes -L > whitelist.T2T_CHM13.bed
```

## NA12878 contigs from HGSVC
```
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/release/v1.0/assemblies/20200628_HHU_assembly-results_CCS_v12/assemblies/phased/v12_NA12878_giab_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/release/v1.0/assemblies/20200628_HHU_assembly-results_CCS_v12/assemblies/phased/v12_NA12878_giab_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta
```

## HG002 contigs from HPRC
```
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG002/assemblies/hifiasm_v0.14.1_raw/HG002-full-0.14.1.pat.fa.gz
```

## T2T-HG002 Genome, v1.1
```
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.1.fasta.gz
gunzip hg002v1.1.fasta.gz
awk -v RS='>' '$1 ~ "PATERNAL" {print RS $0}' hg002v1.1.fasta | sed 's/_PATERNAL//' > HG002_paternal.chroms.fa
```

## NA12878 T2T long-read genome from https://www.biorxiv.org/content/10.1101/2024.08.05.606142v1.full.pdf
```
aws s3 cp --no-sign-request s3://platinum-pedigree-data/assemblies/NA12878/verkko/1.3.1/assembly.haplotype1.fasta NA12878.haplotype1.fasta NA12878.1.fasta
aws s3 cp --no-sign-request s3://platinum-pedigree-data/assemblies/NA12878/verkko/1.3.1/assembly.haplotype2.fasta NA12878.haplotype2.fasta NA12878.2.fasta
# Add chromosome names by aligning to reference
minimap2 -x asm5 -a -t 24 /data/genomes/T2T-CHM13v2.0/GenBank.GCA_009914755.4/GCA_009914755.4_T2T-CHM13v2.0_genomic.chroms.fa NA12878.1.fasta | samtools view -h -b > NA12878.1.CHM13.bam
minimap2 -x asm5 -a -t 24 /data/genomes/T2T-CHM13v2.0/GenBank.GCA_009914755.4/GCA_009914755.4_T2T-CHM13v2.0_genomic.chroms.fa NA12878.2.fasta | samtools view -h -b > NA12878.2.CHM13.bam
./scripts/rename_NA12878_fasta_chroms.py -f NA12878.1.fasta -b NA12878.1.CHM13.bam > NA12878.1.chroms.fasta
./scripts/rename_NA12878_fasta_chroms.py -f NA12878.2.fasta -b NA12878.2.CHM13.bam > NA12878.2.chroms.fasta
# Extract complete chromosomes into their own fasta
./scripts/extract_complete_NA12878_chroms.py -f NA12878.1.chroms.fasta > NA12878.1.chroms.complete.fasta
./scripts/extract_complete_NA12878_chroms.py -f NA12878.2.chroms.fasta > NA12878.2.chroms.complete.fasta
```
