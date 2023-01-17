## Flowchart of data analysis

1) Run commands listed in AlignCount.txt to align and count reads; output is a series of files named "#SAMPLE_exon.count.txt" or "#SAMPLE_CDS.count.txt", containing the read counts for RNA-seq samples and Ribo-seq samples, respectively.
2) Run MergeCounts.R to create the count matrix named "counts_rna_ribo_annotated.txt".
3) Run DESeq2_deltaTE.R to calculate differences in translation efficiency of endogenous mRNAs between *Mael*-mutant and control testes (output file "TE_mael_vs_het_DESeq2.txt"), normalize read counts for RNA-seq (output file "rna_mael_vs_het_DESeq2.txt") and Ribo-seq datasets (output file "ribo_mael_vs_het_DESeq2.txt"), and generate the graph shown in Figure 6D of De Luca C, Gupta A, and Bortvin A. (2023) Ribonucleoprotein condensation driven by retrotransposon LINE-1 sustains RNA integrity and translation in mouse spermatocytes. bioRxiv https://doi.org/10.1101/2023.01.09.523313.
