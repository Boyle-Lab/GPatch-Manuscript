#!/bin/bash

# Postprocess chains produced by the nf-LO pipeline to convert to liftover chains.
# Based on steps described here:
# https://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=1390827233_st3mhvM8SnGwC9IK4FZ9ysMqCN54&db=hg38&c=chrX&g=chm13LiftOver

INPUT_CHAIN=$1
QUERY_NAME=$2
TARGET_NAME=$3
OUT_DIR=$4
QUERY_FASTA=$5
TARGET_FASTA=$6

QUERY_2BIT=$(dirname $QUERY_FASTA)/$(basename $QUERY_FASTA ."${QUERY_FASTA##*.}").2bit
TARGET_2BIT=$(dirname $TARGET_FASTA)/$(basename $TARGET_FASTA ."${TARGET_FASTA##*.}").2bit

if ! [ -f QUERY_2BIT ];then
    # Need to generate twobit for the query genome
    faToTwoBit $QUERY_FASTA $QUERY_2BIT
fi

if ! [ -f QUERY_2BIT ];then
    # Need to generate twobit for the target genome
    faToTwoBit $TARGET_FASTA $TARGET_2BIT
fi

SPLIT_CHAIN=$QUERY_NAME.$TARGET_NAME.split.chain
SPLIT_PAF=$QUERY_NAME.$TARGET_NAME.split.paf
TRIMMED_PAF=$QUERY_NAME.$TARGET_NAME.trimmed.paf
TRIMMED_CHAIN=$QUERY_NAME.$TARGET_NAME.chain
INV_CHAIN=$TARGET_NAME.$QUERY_NAME.chain
FINAL_CHAIN=$QUERY_NAME.$TARGET_NAME.over.chain
FINAL_INV_CHAIN=$TARGET_NAME.$QUERY_NAME.over.chain

module load chaintools
module load rustybam
echo "Splitting chains"
split.py -c $INPUT_CHAIN -o $OUT_DIR/$SPLIT_CHAIN
echo "Converting to paf"
# For some reason, to_paf seems to reverse query and target designations.
to_paf.py -c $OUT_DIR/$SPLIT_CHAIN -q $TARGET_FASTA -t $QUERY_FASTA -o $OUT_DIR/$SPLIT_PAF
echo "Trimming paf aligments"
cat $OUT_DIR/$SPLIT_PAF | rb break-paf --max-size 10000  | rb trim-paf -r | rb invert | rb trim-paf -r | rb invert > $OUT_DIR/$TRIMMED_PAF

echo "Creating chains from paf trimmed paf"
module load pafTools
paf2chain -i $OUT_DIR/$TRIMMED_PAF > $OUT_DIR/$TRIMMED_CHAIN
invert.py -c $OUT_DIR/$TRIMMED_CHAIN -o $OUT_DIR/$INV_CHAIN
echo "Adding scores to $FINAL_CHAIN"
chainMergeSort $OUT_DIR/$TRIMMED_CHAIN | chainScore stdin $QUERY_2BIT $TARGET_2BIT $OUT_DIR/$FINAL_CHAIN
echo "Adding scores to $FINAL_INV_CHAIN"
chainMergeSort $OUT_DIR/$INV_CHAIN | chainScore stdin $TARGET_2BIT $QUERY_2BIT $OUT_DIR/$FINAL_INV_CHAIN
