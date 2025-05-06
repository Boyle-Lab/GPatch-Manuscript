# Simulation Analysis

This folder contains scripts used for the simulated data analyses performed in NA12878 and HG002 subfolders.

### Compiling overall stats for no-indel set
```
paste NA12878/GPatch/CHM13.NA12878.1.contigs.ordering-stats.txt <(awk '{for (i=2; i<NF; i++) printf $i "\t"; print $NF}' NA12878/GPatch/CHM13.NA12878.1.contigs.editdist-stats.txt) | column -s $'\t' -t | awk -v SAMPLE="NA12878" -v METHOD="GPatch" '{if (NR == 1) {printf "Sample\tMethod\t%s\tpct_assembly_localized\n", $0} else {printf "%s\t%s\t%s\t-1\n", SAMPLE, METHOD, $0}}' > simulation_stats.combined.txt

paste NA12878/RagTag/no_indel/*.ordering-stats.txt <(awk '{for (i=2; i<NF; i++) printf $i "\t"; print $NF}' NA12878/RagTag/no_indel/*.editdist-stats.txt) | column -s $'\t' -t | awk -v SAMPLE="NA12878" -v METHOD="RagTag" '{if (NR > 1) printf "%s\t%s\t%s\t-1\n", SAMPLE, METHOD, $0}' >> simulation_stats.combined.txt

paste HG002/GPatch/CHM13.HG002.1.contigs.ordering-stats.txt <(awk '{for (i=2; i<NF; i++) printf $i "\t"; print $NF}' HG002/GPatch/CHM13.HG002.1.contigs.editdist-stats.txt) | column -s $'\t' -t | awk -v SAMPLE="HG002" -v METHOD="GPatch" '{if (NR > 1) printf "%s\t%s\t%s\t-1\n", SAMPLE, METHOD, $0}' >> simulation_stats.combined.txt

paste HG002/RagTag/no_indel/*.ordering-stats.txt <(awk '{for (i=2; i<NF; i++) printf $i "\t"; pri
nt $NF}' HG002/RagTag/no_indel/*.editdist-stats.txt) | column -s $'\t' -t | awk -v SAMPLE="HG002" -v METHOD="RagTag" '{if (NR > 1) printf "%s\t%s\t%s\t-1\n", SAMPLE, METHOD, $0}' >> simulation_stats.combined.txt
```

### Compiling overall stats for	SURVIVOR set
```
paste NA12878/GPatch/CHM13.NA12878.1.SURVIVOR.contigs.ordering-stats.txt <(awk '{for (i=2; i<NF; i++) printf $i "\t"; print $NF}' NA12878/GPatch/CHM13.NA12878.1.SURVIVOR.contigs.editdist-stats.txt) | column -s $'\t' -t | awk -v SAMPLE="NA12878" -v METHOD="GPatch" '{if (NR == 1) {printf "Sample\tMethod\t%s\tpct_assembly_localized\n", $0} else {printf "%s\t%s\t%s\t-1\n", SAMPLE, METHOD, $0}}' > simulation_stats.SURVIVOR.combined.txt

paste HG002/GPatch/CHM13.HG002.1.SURVIVOR.contigs.ordering-stats.txt <(awk '{for (i=2; i<NF; i++) printf $i "\t"; print $NF}' HG002/GPatch/CHM13.HG002.1.SURVIVOR.contigs.editdist-stats.txt) | column -s $'\t' -t | awk -v SAMPLE="HG002" -v METHOD="GPatch" '{if (NR > 1) printf "%s\t%s\t%s\t-1\n", SAMPLE, METHOD, $0}' >> simulation_stats.SURVIVOR.combined.txt
```

### Produce bar plots for figure.
```
./do_barplots.Rscript simulation_stats.combined.txt "simulation_stats.combined"
```

### Data for Tables S2 and S3
Genome lengths in Tables S2 and S3 can be found in {NA12878,HG002}/{GPatch,RagTag}/*.stats
Counts required to calculate patched genome contig coverage %, and % assembly localized (Table S3) are also found in {NA12878,HG002}/{GPatch,RagTag}/*.stats
Remaining rows, with the exception of runtime and max memory, are from simulation_stats.*.combined.txt
