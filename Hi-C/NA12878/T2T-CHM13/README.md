# Hi-C Data Analysis

This folder contains scripts used for the analysis of NA12878 Hi-C data using T2T-CHM13 as the reference.

### Dependencies
* Juicer Tools (https://github.com/aidenlab/JuicerTools)
* SAMtools (https://www.htslib.org/)
* R (https://www.r-project.org/)
* JCuda (http://www.jcuda.org/jcuda/JCuda.html)
* CUDA Toolkit (https://developer.nvidia.com/cuda-toolkit)
* BEDtools (https://bedtools.readthedocs.io/en/latest/)

### Gathering basic alignment stats
```
samtools flagstat -@ 24 bam/SRR1658592.bam
samtools flagstat -@ 24 bam/SRR1658572.bam
```

### Gather mapping qualities and calculate averages
```
samtools view SRR1658572.bam | awk '{print $5}' > SRR1658572.quals.txt
samtools view SRR1658592.bam | awk '{print $5}' > SRR1658592.quals.txt
```
These vectors were read into R and averaged per each replicate:
```
qs1 = scan("SRR1658572.quals.txt")
qs2 = scan("SRR1658592.quals.txt")
summary(qs1)
summary(qs2)
```

### Gather Pairs stats
```
# Final, filtered, UU, UR, and RU pairs.
zcat pairsam-filter/T2T-CHM13.dedup.pairs.gz | grep -v "#" | wc -l
# Unmapped, duplicate, and other non-UU, -UR, and -RU pairs.
zcat pairsam-filter/T2T-CHM13.unmapped.sam.pairs.gz | grep -v "#" | wc -l
# Final, filtered pairs with mate on different chromosome
zcat pairsam-filter/T2T-CHM13.dedup.pairs.gz | grep -v "#" | awk '{if ($2 != $4) print $0}' | wc -l
```
Total pairs is the sum of unmapped and filtered pairs.

### Loop prediction at 5kb and 10kb resolution with Juicer Tools
```
make postProcessJuicer
```

# Gather stats on loop numbers and loop splitting across contigs.
```
cd juicer/T2T-CHM13_loops

# Number of loops
wc -l enriched_pixels_10000.bedpe
wc -l enriched_pixels_5000.bedpe

# Average loop length
cat enriched_pixels_5000.bedpe | awk 'BEGIN {N=1}; {if (NR > 1) {printf "chr%s\t%d\t%d\tloop%d\n", $1, $2, $6, N; N++}}' > loop_intervals.5kb.bed
cat loop_intervals.5kb.bed | awk 'BEGIN {SUM=0;N=0}; {SUM+=($3-$2);N++}; END {print SUM/N}'

cat enriched_pixels_10000.bedpe | awk 'BEGIN {N=1}; {if (NR > 1) {printf "chr%s\t%d\t%d\tloop%d\n", $1, $2, $6, N; N++}}' > loop_intervals.10kb.bed
cat loop_intervals.10kb.bed | awk 'BEGIN {SUM=0;N=0}; {SUM+=($3-$2);N++}; END {print SUM/N}'
