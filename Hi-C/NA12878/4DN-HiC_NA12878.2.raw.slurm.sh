#!/bin/bash

#SBATCH --partition=largemem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=1503g
#SBATCH --time=07-00:00:00
#SBATCH --output=/home/%u/GPatch_Hi-C/%x-%j.log
#SBATCH --account=apboyle3
#SBATCH --mail-user=adadiehl@umich.edu
#SBATCH --mail-type=END

# SLURM script for Hi-C analysis of Rao and Huntley data for NA12878,
# mapped to the genome defined in the given Makefile..

MAKEFILE_PATH=$1

module load Bioinformatics
module load bwa
module load fastqc
module load samtools

conda init bash
source ~/.bashrc
conda activate 4DN_HiC_env

mkdir -pv /tmp_data/$SLURM_JOB_ID
cd /tmp_data/$SLURM_JOB_ID
cp $MAKEFILE_PATH Makefile

/usr/bin/time -v make

# Copy results to the output directory
