# Hi-C Data Analysis

This folder contains scripts used for the analysis of NA12878 Hi-C data using GPatch NA12878, T2T-CHM13, and HGSVC NA12878 (unpatched) assemblies as the reference.

### Dependencies


### Initial processing with 4DN pipeline and Juicer Tools
Initial mapping, pairs extraction, filtering, Hi-C matrix construction, and normalization steps are automated in the `Makefile.*` files. These were processed on memory-optimized cluster nodes allocated with 12 processor cores and 500G memory, with jobs provisioned through Slurm:
```
sbatch <slurm_script.target.slurm.sh> <Makefile.target>
```

### Loop prediction at 5kb and 10kb resolution with Juicer Tools
Loop prediction steps were performed on a local, GPU-enabled server, with Makefiles with the '.loops' suffix used to automate the process. Steps are documented in each target's subdirectory.

