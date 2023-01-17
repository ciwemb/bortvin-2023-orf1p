## Flowchart of ribosome profiling data analysis

1) Run commands listed in AlignCount.txt to align and count reads; output is a series of files named "#SAMPLE_exon.count.txt" or "#SAMPLE_CDS.count.txt" containing the read counts for RNA-seq samples and Ribo-seq samples, respectively
2) Run MergeCounts.R to create the count matrix named "counts_rna_ribo_annotated.txt"
3) Run DESeq2_deltaTE.R to calculate differences in translation efficiency between *Mael*-mutant and control testes (output file "TE_mael_vs_het_DESeq2.txt"), normalize read counts for RNA-seq datasets (output file "rna_mael_vs_het_DESeq2.txt") and Ribo-seq datasets (output file "ribo_mael_vs_het_DESeq2.txt"), and generate the graph shown in Figure 6D of 
