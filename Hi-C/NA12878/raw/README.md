# Hi-C Data Analysis

This folder contains scripts used for the analysis of NA12878 Hi-C data using HGSVC NA12878 as the reference.

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
zcat pairsam-filter/NA12878_raw.dedup.pairs.gz | grep -v "#" | wc -l
# Unmapped, duplicate, and other non-UU, -UR, and -RU pairs.
zcat pairsam-filter/NA12878_raw.unmapped.sam.pairs.gz | grep -v "#" | wc -l
# Final, filtered pairs with mate on different chromosome
zcat pairsam-filter/NA12878_patched.dedup.pairs.gz | grep -v "#" | awk '{if ($2 != $4) print $0}' | wc -l
```
Total pairs is the sum of unmapped and filtered pairs.
