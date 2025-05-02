#!/bin/bash

# Automates PAF alignment of patched genome against target assembly,
# using only chromosomes that are complete in the target assembly.

PATCHED_FASTA=$1
TARGET_FASTA=$2
TARGET_NAME=$3
PREFIX=$4

SCRIPTS_PATH=/data/projects/adadiehl/genome_patching/GPatch/scripts

# Filter chromosomes complete in NA12878 from patched output
#../../scripts/extract_NA12878-complete_patched_chroms.py -p $PATCHED_FASTA -t $TARGET_FASTA > $PREFIX.patched.target-complete-chroms.fasta

# Alignment to reference and dot-plots.
echo "Mapping the patched genome against the target..."
#minimap2 -x asm20 -t 24 $TARGET_FASTA $PREFIX.patched.target-complete-chroms.fasta > $PREFIX.patched.target-complete-chroms.$TARGET_NAME.paf
echo "Generating dot plots"
$SCRIPTS_PATH/create_dotplots.Rscript $PREFIX.patched.target-complete-chroms.$TARGET_NAME.paf $PREFIX.target-complete-chroms $TARGET_NAME
