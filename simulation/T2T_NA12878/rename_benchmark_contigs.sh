#!/usr/bin/env bash
#
# NA12878.1.chroms.fasta labels every benchmark contig with its chromosome
# only (">chr1 mat-0000003", ">chr1 mat-0000006", ...), so IDs are not unique
# and cannot be used as minimap2 targets. Write a copy with unique IDs
# <chrom>.<descriptor> (e.g. chr1.mat-0000003) and index it. The chromosome
# is the part of the ID before the first ".".
#
# Usage: ./rename_benchmark_contigs.sh [in.fasta] [out.fa]
set -euo pipefail
IN=${1:-/data/projects/adadiehl/genome_patching/manuscript/data/NA12878.1.chroms.fasta}
OUT=${2:-NA12878.1.chroms.renamed.fa}
awk '/^>/{c=substr($1,2); id=c"."$2; if(id in seen){print "duplicate ID "id > "/dev/stderr"; exit 1}
     seen[id]=1; print ">"id" "c" "$2; next} {print}' "$IN" > "$OUT"
/data/software/Anaconda3/2022.10/envs/GPatch_env/bin/python3 -c "import pysam, sys; pysam.faidx(sys.argv[1])" "$OUT"
