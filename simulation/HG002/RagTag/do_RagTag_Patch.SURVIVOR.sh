#!/bin/bash

# Run RagTag Patch using the given reference and query genomes and
# convert results to follow GPatch conventions.

ASSEMBLY_FASTA=$1
REFERENCE_FASTA=$2
REFERENCE_NAME=$3
PREFIX=$4
OUT_DIR=$5
CONCAT_FASTA=$6

SCRIPTS_PATH=/data/projects/adadiehl/genome_patching/GPatch/scripts

# Run RagTag Patch
/usr/bin/time -v ragtag.py patch --aligner minimap2 --mm2-params "-x asm20 -c -t 24" -o $OUT_DIR -f 100 -s 10000 $ASSEMBLY_FASTA $REFERENCE_FASTA

# Convert AGP output to contigs.bed
../../RagTag_agp_to_contigs_bed.py -r $OUT_DIR/ragtag.patch.rename.agp -c $OUT_DIR/ragtag.patch.ctg.agp -p $OUT_DIR/ragtag.patch.agp > $OUT_DIR/ragtag.patch.contigs.bed

# Filter out only scaffolds from FASTA output
../../process_ragtag_fasta.py -f $OUT_DIR/ragtag.patch.fasta -r $OUT_DIR/ragtag.patch.rename.agp -p $OUT_DIR/ragtag.patch.agp -c $OUT_DIR/ragtag.patch.ctg.agp > $OUT_DIR/ragtag.patch.scf.fasta

# Pull out some quick assembly stats
printf "pre-break contig base count: %d\n" $(awk 'BEGIN {len=0} {len += ($3-$2)} END {print len}' $OUT_DIR/ragtag.patch.contigs.bed) > $OUT_DIR/ragtag.patch.stats
#printf "pre-break patch base count: %d\n" $(awk 'BEGIN {len=0} {len += ($3-$2)} END {print len}' $OUT_DIR/$PREFIX.patches.bed) >> $OUT_DIR/ragtag.patch.stats
printf "Initial unpatched assembly genome size (bp): %d\n" $(grep -v ">" $ASSEMBLY_FASTA | awk 'BEGIN{N=0}{N+=length($0)}END{print N}') >> $OUT_DIR/ragtag.patch.stats
printf "pre-break total patched genome size (bp): %d\n" $(grep -v ">" $OUT_DIR/ragtag.patch.scf.fasta | awk 'BEGIN{N=0}{N+=length($0)}END{print N}') >> $OUT_DIR/ragtag.patch.stats

# Alignment to reference and dot-plots.
minimap2 -x asm20 -t 24 $CONCAT_FASTA $OUT_DIR/ragtag.patch.scf.fasta > $OUT_DIR/ragtag.patch.$REFERENCE_NAME.paf
$SCRIPTS_PATH/create_dotplots.Rscript $OUT_DIR/ragtag.patch.$REFERENCE_NAME.paf $PREFIX $REFERENCE_NAME
mv *.pdf $OUT_DIR
