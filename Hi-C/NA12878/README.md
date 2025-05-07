# Hi-C Data Loop Comparison between GPatch NA12878 and T2T-CHM13

This document describes steps taken to directly compare loop loci between Hi-C data mapped to GPatch NA12878 or T2T-CH13 assemblies.

### Dependencies
* BEDtools (https://bedtools.readthedocs.io/en/latest/)
* Kent Tools (https://github.com/ucscGenomeBrowser/kent)

### First decompose loops into individual anchors
```
cd patched/juicer/NA12878_patched_loops

cat enriched_pixels_5000.bedpe | awk 'BEGIN {N=1}; {if (NR > 1) {printf "chr%s\t%d\t%d\tloop%d_1\nchr%s\t%d\t%d\tloop%d_2\n", $1, $2, $3, N, $4, $5, $6, N; N++}}' > loop_anchors.5kb.bed

cat enriched_pixels_10000.bedpe | awk 'BEGIN {N=1}; {if (NR > 1) {printf "chr%s\t%d\t%d\tloop%d_1\nchr%s\t%d\t%d\tloop%d_2\n", $1, $2, $3, N, $4, $5, $6, N; N++}}' > loop_anchors.10kb.bed
```

### Use previously-prepared liftover chains to convert GPatch NA12878 to T2T-CHM13 frame
```
liftOver loop_anchors.5kb.bed ../../../../../biological/NA12878/liftover/alignment_results_mm2/chainnet/NA12878.T2T-CHM13.over.chain loop_anchors.5kb.T2T-CHM13.bed loop_anchors.5kb.unmapped.bed
liftOver loop_anchors.10kb.bed ../../../../../biological/NA12878/liftover/alignment_results_mm2/chainnet/NA12878.T2T-CHM13.over.chain loop_anchors.10kb.T2T-CHM13.bed loop_anchors.10kb.unmapped.bed
```

### Intersect lifted anchors with T2T-CHM13 to isolate shared anchors.
```
bedtools intersect -a loop_anchors.5kb.T2T-CHM13.bed -b ../../../../T2T-CHM13/juicer/T2T-CHM13_loops/loop_anchors.5kb.bed -u > loop_anchors.5kb.T2T-CHM13.shared.bed
bedtools intersect -a loop_anchors.10kb.T2T-CHM13.bed -b ../../../../T2T-CHM13/juicer/T2T-CHM13_loops/loop_anchors.10kb.bed -u > loop_anchors.10kb.T2T-CHM13.shared.bed
```

### Reconstruct loop records from the sets of lifted anchors
```
../../../../reconstitute_loops.py -b loop_anchors.5kb.T2T-CHM13.shared.bed > loops.5kb.T2T-CHM13.shared.bed
../../../../reconstitute_loops.py -b loop_anchors.5kb.T2T-CHM13.bed > loops.5kb.T2T-CHM13.bed

../../../../reconstitute_loops.py -b loop_anchors.10kb.T2T-CHM13.shared.bed > loops.10kb.T2T-CHM13.shared.bed
../../../../reconstitute_loops.py -b loop_anchors.10kb.T2T-CHM13.bed > loops.10kb.T2T-CHM13.bed
```

### Gather counts
```
# Fully-conserved loops (both anchors are lifted and shared with T2T-CHM13)
grep -v '\.' loops.5kb.T2T-CHM13.shared.bed | wc -l
grep -v '\.' loops.10kb.T2T-CHM13.shared.bed | wc -l

# Partially-conserved (both anchors are lifted but only one shared with T2T-CHM13)
grep '\.' loops.5kb.T2T-CHM13.shared.bed | wc -l
grep '\.' loops.10kb.T2T-CHM13.shared.bed | wc -l
```

The NA12878-only loop count was taken as the difference between total loops and the sum of fully-conserved and partially-conserved loops.
