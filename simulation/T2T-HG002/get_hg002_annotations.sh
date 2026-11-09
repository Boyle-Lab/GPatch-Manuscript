#!/usr/bin/env bash
#
# Download HG002 Q100 (hg002v1.1) paternal repeat / segmental duplication
# annotations into annotation/raw/ and convert the RepeatMasker bigBed to BED.
# Chromosome names are left as chrN_PATERNAL; breakpoint_repeat_enrichment.py
# strips the suffix to match HG002_paternal.chroms.fa.

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

S3=https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/annotation
# pyBigWig (bigBed reader) is available in HiCExplorer_env
BBPY=/data/software/Anaconda3/2022.10/envs/HiCExplorer_env/bin/python3
RAW=annotation/raw
mkdir -p $RAW annotation

fetch() {
    if [[ ! -s $2 ]]; then
        curl -fsSL --retry 3 -o "$2.tmp" "$1"
        mv "$2.tmp" "$2"
    fi
}

fetch $S3/segdups/hg002v1.1.pat.SDs.013025.bed                        $RAW/hg002v1.1.pat.SDs.013025.bed
fetch $S3/repeatmasker/hg002v1.1.paternal.bb                           $RAW/hg002v1.1.paternal.rmsk.bb
fetch $S3/repeatmasker/README.txt                                      $RAW/repeatmasker.README.txt
fetch $S3/centromere/hg002v1.1_v2.0/hg002v1.1.cenSatv2.0.PAT.bed       $RAW/hg002v1.1.cenSatv2.0.PAT.bed
fetch $S3/seqfeatures/hg002v1.1.simpleRepeat.newref.sort.bed           $RAW/hg002v1.1.simpleRepeat.bed
fetch $S3/seqfeatures/v1.1.mnr10.pat.bed                               $RAW/v1.1.mnr10.pat.bed
fetch $S3/seqfeatures/README.txt                                       $RAW/seqfeatures.README.txt

if [[ ! -s $RAW/hg002v1.1.paternal.rmsk.bed ]]; then
    $BBPY - $RAW/hg002v1.1.paternal.rmsk.bb > $RAW/hg002v1.1.paternal.rmsk.bed.tmp <<'PY'
import sys, pyBigWig
bb = pyBigWig.open(sys.argv[1])
for chrom, L in bb.chroms().items():
    for s, e, rest in bb.entries(chrom, 0, L) or []:
        sys.stdout.write("%s\t%d\t%d\t%s\n" % (chrom, s, e, rest))
PY
    mv $RAW/hg002v1.1.paternal.rmsk.bed.tmp $RAW/hg002v1.1.paternal.rmsk.bed
fi

md5sum $RAW/* > annotation/raw.md5
date -Is > annotation/download_date.txt
