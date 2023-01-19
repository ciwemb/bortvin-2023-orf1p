#!/bin/sh
#$ -N starindex
#$ -o starindex.txt
#$ -S /bin/sh
#$ -cwd
#$ -q zappa,largememory
#$ -l mem_free=70G,h_vmem=80G
#$ -M agupta52@jhu.edu

module load sharedapps star/2.6.0c ucsctools

#Download RepeatMasker track in bed format using UCSC table-browser
bedToGenePred mm10_UCSC_rmsk.bed mm10_UCSC_rmsk.GenePred
genePredToGtf file mm10_UCSC_rmsk.GenePred mm10_UCSC_rmsk.gtf
./mergeGTF.pl  
STAR --runMode genomeGenerate --genomeDir /genomes/mouse/mm10/STAR/ --sjdbGTFfile mm10_UCSC_refSeq_rmsk.gtf --genomeFastaFiles mm10.fa
