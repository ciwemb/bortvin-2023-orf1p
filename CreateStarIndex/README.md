This folder contains the commands to generate the genome indexes files. 

Files needed are: the reference genome sequence ("mm10.fa") and the annotation file ("mm10_UCSC_refSeq_rmsk.gtf"). 
"mm10_UCSC_refSeq_rmsk.gtf" file is generated by mergeGTF.pl, where "mm10_UCSC_noAlt.gtf" is the gene annotation file obtained using iGenomes (http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz) minus the genes that are present on alternate/random chromosomes (`grep -v _random iGenomes/Genes/genes.gtf | grep -v chrUn > mm10_UCSC_noAlt.gtf`), “mm10_UCSC_rmsk.gtf” is the RepeatMasker track, and "chrlist.txt" is a text file containing the chromosome names (chr1, chr2 etc.) with one entry per row.

The generated indexed genome is utilized in the STAR mapping steps of both RNA-seq and Ribosome profiling analyses reported in this repository.
