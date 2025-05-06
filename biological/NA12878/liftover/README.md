# GPatch NA12878 to T2T-CHM13 Alignment Chains
Steps for producing aligment chains between GPatch NA12878 and T2T-CHM13.

## Dependencies
* minimap2 (https://github.com/lh3/minimap2
* Kent Tools (https://github.com/ucscGenomeBrowser/kent)
* nextflow (https://www.nextflow.io/)
* chaintools (https://github.com/milkschen/chaintools)
* rustybam (https://github.com/mrvollger/rustybam)
* chain2paf (https://github.com/AndreaGuarracino/chain2paf)
* paf2Chain (https://github.com/AndreaGuarracino/paf2chain)

Final alignment chains to convert GPatch NA12878 to T2T-CHM13 coordinates and vice-versa are in the alignment_results_mm2/chainnet/ directory:
```
ls -alth alignment_results_mm2/chainnet/
total 75M
-rw-rw-rw-  1 adadiehl boyle-lab 5.8M Apr 23 09:59 T2T-CHM13.NA12878.over.chain
drwxrwsrwx+ 2 adadiehl boyle-lab 4.0K Apr 23 09:59 .
-rw-rw-rw-  1 adadiehl boyle-lab 5.8M Apr 23 09:59 NA12878.T2T-CHM13.over.chain
```

### nf-LO nextflow pipeline for initial chain production
```
nextflow run evotools/nf-LO --target ../../../data/GCA_009914755.4_T2T-CHM13v2.0_genomic.chroms.fa --source ../NA12878.2.cbreak_2.patched.fasta --aligner minimap2 --custom "-x asm20 -t 18" --distance near -profile conda --outdir ./alignment_results_mm2 --resume --max_memory 64.GB --max_cpus 24 --mm2_full_alignment
```

### Chain post-processing into final liftover chains
```
bash do_chain_postprocessing.sh alignment_results_mm2/chainnet/liftover.chain NA12878 T2T-CHM13 alignment_results_mm2/chainnet ../NA12878.2.cbreak_2.patched.fasta ../../../data/GCA_009914755.4_T2T-CHM13v2.0_genomic.chroms.fa
```
