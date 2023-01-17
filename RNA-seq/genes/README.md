## Flowchart of data analysis

1) Run commands listed in 1-Align.txt to align and count reads; output is a series of files named #SAMPLE_mm10_refSeq_rmskReadsPerGene.out.tab.
2) Run 2-mergeGeneCounts.R to create the count matrix named "GeneCounts_all_samples.txt".
3) Run 3-DESeq2.R to normalize read counts. The comparison between IP and BO samples is reported as an example (output file "IP_comb_vs_BO_comb_DESeq2.txt").

DESeq2 normalization was also run for the comparisons IP vs INPUT and IP vs TOTAL (output files "IP_comb_vs_INPUT_comb_DESeq2.txt" and "IP_comb_vs_TOTAL_comb_DESeq2.txt", respectively). All DESeq2 output files were submitted to NCBI's Gene Expression Omnibus. Genes with log2FoldChange > 1 and padj < 0.05 in every comparison were considered significant. The intersection group (n=1347) was considered "strongly enriched" in IP samples and it is shown in red in Figure 3B, 3C and 3D of De Luca C, Gupta A, and Bortvin A. (2023) Ribonucleoprotein condensation driven by retrotransposon LINE-1 sustains RNA integrity and translation in mouse spermatocytes. bioRxiv https://doi.org/10.1101/2023.01.09.523313.
