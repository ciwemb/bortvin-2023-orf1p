## Flowchart of data analysis

1) Run commands listed in 1-Align.txt to align (STAR) and re-assign ambiguously mapped reads (Telescope).
2) Run 2-mergeTelescopeCounts.R to merge Telescope counts and create a master table named "GeneCounts_telescope_allSamples_annotated.txt".
3) Run 3-collapseCount.pl to collapse genomic repeats read counts based on same repName, repClass or repFamily.
4) Run 4-upperQuantileNormalize.R to normalize counts using EBSeq package in R.
5) Run 5-makeHeatmap.R to generate heatmaps shown in Figure 3A (for a selected subset of repeats) and in Figure S3A (for all repeats with an average expression level higher than ancient LINE-1 family Lx) of De Luca C, Gupta A, and Bortvin A. (2023) Ribonucleoprotein condensation driven by retrotransposon LINE-1 sustains RNA integrity and translation in mouse spermatocytes. bioRxiv https://doi.org/10.1101/2023.01.09.523313.
