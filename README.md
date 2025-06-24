# GPatch-Manuscript
Data retrieval and Analysis steps for the GPatch manuscript:

Fast and Accurate Draft Genome Patching with GPatch
Adam Diehl, Alan Boyle
bioRxiv 2025.05.22.655567; doi: https://doi.org/10.1101/2025.05.22.655567

# Overview
Recent advancements in sequencing technologies have yielded numerous long-read draft genomes, promising to enhance our understanding of genomic variation. However, draft genomes are typically highly fragmented, posing significant challenges for functional genomics. We introduce GPatch, a tool that constructs chromosome-scale pseudoassemblies from fragmented drafts using alignments to a reference genome. GPatch produces complete, accurate, gap-free assemblies preserving over 95% of nucleotides from draft genomes. We show that GPatch assemblies can be used as references for Hi-C data analysis, whereas draft assemblies cannot. Until complete genome assembly becomes routine, GPatch presents a necessary tool for maximizing the utility of draft genomes.

# Organization
Data are organized according to the results sections in which they are presented:
* Simulated data analysis results are in the "simulation" subdirectory
* Patching of NA21878 and HG002 draft genomes are presented in the "biological" subdirectory
* Analysis of NA12878 Hi-C data with GPatch NA12878 and T2T-CHM13 assemblies are presented in the "Hi-C" subdirectory.
* Source data and accessions are housed in the "data" subdirectory

Each subdirectory contains a README with sufficient detail to reproduce analyses as presented in the manuscript.

# Correspondence
Please contact Adam Diehl (adadiehl@umich.edu) with any questions or issues related to this repository.
