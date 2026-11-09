#!/usr/bin/env bash
#
# Simulate a contig-level assembly with known order/orientation by aligning a
# real contig-level "model" assembly (HG002.1, HGSVC) to a chromosome-scale
# "benchmark" assembly (T2T-HG002 paternal) and using the boundaries of the
# aligned blocks as simulated contig boundaries. Contig sequence is extracted
# from the benchmark genome (reverse-complemented for - strand alignments), so
# breakpoints, overlaps and unfilled gaps between adjacent contigs follow the
# empirical distribution of the model assembly (enriched in repeats where real
# assemblers fail) without needing repeat annotations for the benchmark.
#
# Usage:
#   ./simulate_contigs_from_alignments.sh [-q model.fa] [-r benchmark.fa]
#       [-p prefix] [-t threads] [-m min_mapq] [-l min_len] [-M merge_dist]
#       [-a existing.paf] [-k] [-f]
#
#   -a  use this existing minimap2 PAF instead of <prefix>.paf (e.g. to try
#       different conversion options without re-aligning)
#   -M  merge consecutive colinear same-strand alignment blocks of the same
#       model contig separated by <= merge_dist bp (default: no merging)
#   -k  keep alignments whose benchmark interval is contained within another
#       (default: drop them)
#   -f  force re-running minimap2 even if the PAF already exists
#
# Outputs:
#   <prefix>.paf            raw minimap2 alignments
#   <prefix>.contigs.fa     simulated contigs (named <chrom>_<rank>)
#   <prefix>.contigs.bed    simulated contig locations in the benchmark genome
#   <prefix>.summary.txt    breakpoint/overlap/gap summary
#   <prefix>.log            commands, software versions, input checksums

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA=/data/projects/adadiehl/genome_patching/manuscript/data
ENV_BIN=/data/software/Anaconda3/2022.10/envs/GPatch_env/bin

QUERY="$DATA/HG002-full-0.14.1.pat.fa"
BENCH="$DATA/HG002_paternal.chroms.fa"
PREFIX="T2T-HG002.HG002.1.simulated"
THREADS=32
MIN_MAPQ=0
MIN_LEN=1000
KEEP_CONTAINED=""
MERGE=""
FORCE=0
PAF=""

while getopts "q:r:p:t:m:l:M:a:kfh" opt; do
    case $opt in
        q) QUERY=$OPTARG ;;
        r) BENCH=$OPTARG ;;
        p) PREFIX=$OPTARG ;;
        t) THREADS=$OPTARG ;;
        m) MIN_MAPQ=$OPTARG ;;
        l) MIN_LEN=$OPTARG ;;
        M) MERGE="--merge-dist $OPTARG" ;;
        a) PAF=$OPTARG ;;
        k) KEEP_CONTAINED="--keep-contained" ;;
        f) FORCE=1 ;;
        h|*) sed -n '2,37p' "$0"; exit 1 ;;
    esac
done

MINIMAP2="$ENV_BIN/minimap2"
PYTHON="$ENV_BIN/python3"
LOG="$PREFIX.log"

{
    echo "# $(date -Is)  $(hostname)"
    echo "# command: $0 $*"
    echo "# minimap2 $($MINIMAP2 --version)"
    echo "# python $($PYTHON --version 2>&1 | cut -d' ' -f2); biopython $($PYTHON -c 'import Bio; print(Bio.__version__)')"
    echo "# query (model assembly):   $QUERY"
    echo "# benchmark genome:         $BENCH"
    echo "# git commit: $(git -C "$SCRIPT_DIR" rev-parse --short HEAD 2>/dev/null || echo NA)"
} > "$LOG"

# 1. Align model contigs to the benchmark genome. --secondary=no keeps one
#    placement per query segment; -c computes base-level alignments so block
#    boundaries and AS/NM/de tags are exact.
if [[ -n $PAF ]]; then
    [[ -s $PAF ]] || { echo "PAF $PAF not found" >&2; exit 1; }
    echo "# using supplied PAF $PAF (md5 $(md5sum < "$PAF" | cut -c1-32))" >> "$LOG"
elif PAF="$PREFIX.paf"; [[ $FORCE -eq 1 || ! -s $PAF ]]; then
    CMD="$MINIMAP2 -x asm5 -c --secondary=no -t $THREADS $BENCH $QUERY"
    echo "$CMD > $PAF" >> "$LOG"
    $CMD > "$PAF.tmp" 2> "$PREFIX.minimap2.err"
    mv "$PAF.tmp" "$PAF"
else
    echo "# reusing existing $PAF (use -f to re-align)" >> "$LOG"
fi

# 2. Convert alignment blocks into simulated contigs + BED.
CMD="$PYTHON $SCRIPT_DIR/paf_to_simulated_contigs.py -p $PAF -r $BENCH -o $PREFIX --min-mapq $MIN_MAPQ --min-len $MIN_LEN $MERGE $KEEP_CONTAINED"
echo "$CMD" >> "$LOG"
$CMD 2>> "$LOG"

echo "# finished $(date -Is)" >> "$LOG"
