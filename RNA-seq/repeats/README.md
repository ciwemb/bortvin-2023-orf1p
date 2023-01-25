## Outline of data analysis

1) Run commands listed in 1-Align.txt to align and re-assign ambiguously mapped reads.
2) Run 2-mergeTelescopeCounts.R to merge Telescope counts and create a master table named "GeneCounts_telescope_allSamples_annotated.txt".
3) Run 3-collapseCount.pl to collapse counts based on common repName, repClass or repFamily.
4) Run 4-upperQuantileNormalize.R to normalize counts using EBSeq package in R, and generate the table named "repName_collapsed_counts_upper_quantile_norm.txt" that was submitted to the NCBI's Gene Expression Omnibus repository.
5) Run 5-makeHeatmap.R to generate heatmaps shown in Figures 3A and S3A of De Luca C, Gupta A, and Bortvin A. (2023) Ribonucleoprotein condensation driven by retrotransposon LINE-1 sustains RNA integrity and translation in mouse spermatocytes. bioRxiv https://doi.org/10.1101/2023.01.09.523313.
