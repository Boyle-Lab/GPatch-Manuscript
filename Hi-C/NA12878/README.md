# Hi-C Data Analysis

This folder contains scripts used for the analysis of NA12878 Hi-C data using GPatch NA12878, T2T-CHM13, and HGSVC NA12878 (unpatched) assemblies as the reference.

### Dependencies
* conda (https://anaconda.org/anaconda/conda)
* BWA (https://github.com/lh3/bwa)
* Juicer Tools (https://github.com/aidenlab/JuicerTools)
* SAMtools (https://www.htslib.org/)
* FastQC (https://github.com/s-andrews/FastQC)
* R (https://www.r-project.org/)
* JCuda (http://www.jcuda.org/jcuda/JCuda.html)
* CUDA Toolkit (https://developer.nvidia.com/cuda-toolkit)

### Copy GPatch NA12878 to the central /data folder for convenience.
```
cp ../../../biological/NA12878/NA12878.2.cbreak_2.patched.fasta ../../../data/
```

### Generate chrom.sizes for T2T-CHM13, GPatch NA12878, and HGSVC NA12878.
```
cd ../../../data
cat NA12878.2.cbreak_2.patched.fasta | awk 'BEGIN {len=0}; {if ($1 ~ ">") {if (len > 0) {printf "%d\n", len}; printf "%s\t", $1; len=0} else {len += length($1)} }; END {printf "%d\n", len}' | sed 's/>//' > NA12878.2.cbreak_2.patched.chrom.sizes
cat v12_NA12878_giab_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta | awk 'BEGIN {len=0}; {if ($1 ~ ">") {if (len > 0) {printf "%d\n", len}; printf "%s\t", $1; len=0} else {len += length($1)} }; END {printf "%d\n", len}' | sed 's/>//' > v12_NA12878_giab_pbsq2-ccs_1000-pereg.h2-un.racon-p2.chrom.sizes
cat GCA_009914755.4_T2T-CHM13v2.0_genomic.chroms.fa | awk 'BEGIN {len=0}; {if ($1 ~ ">") {if (len > 0) {printf "%d\n", len}; printf "%s\t", $1; len=0} else {len += length($1)} }; END {printf "%d\n", len}' | sed 's/>//' > GCA_009914755.4_T2T-CHM13v2.0_genomic.chroms.chrom.sizes
```

### Generate BWA indexes for the three assemblies.
```
bwa index NA12878.2.cbreak_2.patched.fasta
bwa index v12_NA12878_giab_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta 
bwa index GCA_009914755.4_T2T-CHM13v2.0_genomic.chroms.fa 
```

### Generate restriction enzyme sites files for the three assemblies.
Juicer tools scripts should be located at the path indicated in $JUICER_TOOLS_PATH.
```
python3 $JUICER_TOOLS_PATH/generate_site_positions.py MboI NA12878.2.cbreak_2.patched NA12878.2.cbreak_2.patched.fasta

python3 $JUICER_TOOLS_PATH/generate_site_positions.py MboI v12_NA12878_giab_pbsq2-ccs_1000-pereg.h2-un.racon-p2 v12_NA12878_giab_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta

python3 $JUICER_TOOLS_PATH/generate_site_positions.py MboI GCA_009914755.4_T2T-CHM13v2.0_genomic.chroms GCA_009914755.4_T2T-CHM13v2.0_genomic.chroms.fa
```
## The following steps are performed on the compute cluster.

### Prerequisites
Make sure BWA, Juicer Tools, and SAMtools, are accessible on your system. All other dependencies will be installed in the conda environment.

### Set up the conda environment for Hi-C processing
conda create -n 4DN_HiC_env -c conda-forge -c bioconda pairtools pairix cooler numpy cython

### Initial processing with 4DN pipeline and Juicer Tools
Initial mapping, pairs extraction, filtering, Hi-C matrix construction, and normalization steps are automated in `Makefile.*`. These files contain target-assembly-specific paths and variables, most-importantly, cluster-accessible paths to the Hi-C read data and fasta target assembly sequences, chrom.sizes, and BWA indexes, and RE sites files, needed for Hi-C processing.

Jobs for each target assembly were processed on memory-optimized cluster nodes allocated with 12 processor cores and 500G memory, with jobs provisioned through Slurm:
```
sbatch <slurm_script.target.slurm.sh> <Makefile.target>
```
Results are written to a location specified within Makefile.target. These were then copied to a local server for loop prediction steps.

Additional processing steps to gather data presented in the manuscript and figures are documented in the patched (GPatch NA12878), raw (HGSVC NA12878), and T2T-CHM13 subdirectories.

### Loop prediction at 5kb and 10kb resolution with Juicer Tools
Loop prediction steps were performed on a local, GPU-enabled server, with Makefiles with the '.loops' suffix used to automate the process. Steps are documented in each target's subdirectory.

