## Flowchart of data analysis

1) Run commands listed in 1-Align.txt to align and count reads; output is a series of files named #SAMPLE_mm10_refSeq_rmskReadsPerGene.out.tab.
2) Run mergeGeneCounts.R to create the count matrix named "GeneCounts_all_samples.txt".
3) Run DESeq2 to normalize read counts (output file "rna_mael_vs_het_DESeq2.txt") , 

and generate the graph shown in Figure 6D of De Luca C, Gupta A, and Bortvin A. (2023) Ribonucleoprotein condensation driven by retrotransposon LINE-1 sustains RNA integrity and translation in mouse spermatocytes. bioRxiv https://doi.org/10.1101/2023.01.09.523313.
