###################################################
### Clear Workspace
rm(list=ls())

###################################################
### Set working directory
setwd("/mnt/sequence/cdeluca/telescope/collapse_count_and_normalize/")
getwd()

###########################
### initialize command line parsing
args <- commandArgs(TRUE)


###########################
### set values to null
inputFile <- "repName_collapsed_count_all_samples_telescope.txt" # normalized on May 29th, 2020


###########################
### go through the command line arguments and look for matches to the
#flags that we've defined. When those are found, set the appropriate value accordingly.
i <- 1
while (i <  length(args))
{
  opt <- args[i]
  val <- args[i+1]
  if (opt == "-i") { inputFile <- val }
  i <- i+2
}


###########################
### make sure that the value is set, and if not, complain and quit
if (inputFile == "")
{
   message ("Usage: Rscript upperQunatileNormalize.R -i countData")
   q()
}


###################################################
### Load libraries
library(EBSeq)


###################################################
### read the data
countData <- read.table(inputFile, header=T, sep="\t", stringsAsFactors=F)
myColumnName <- colnames(countData)[1]
row.names(countData) <- countData[,1]
countData <- countData[,-1]


###################################################
### Quantile normalize the count
sizes <- QuantileNorm(countData, 0.75) # to obtain size factors for quantile normalization
countData.norm <- t(as.matrix(countData)) # traspose the rows and columns of the matrix countData
countData.norm <- countData.norm/sizes # divide counts by size factors
countData.norm <- round(countData.norm,2) # round to 2 decimal places
countData.norm <- as.data.frame(t(countData.norm)) # traspose again as data frame

countData.norm[,ncol(countData.norm)+1] <- row.names(countData.norm)  # add one column with row names (repName)
colnames(countData.norm)[ncol(countData.norm)] <- myColumnName

countData.norm <- countData.norm[,c(ncol(countData.norm),1:(ncol(countData.norm)-1))]


###################################################
### write the quantile normalized values to a file
outFile <- sub("count", "count_upper_quantile_normalized", inputFile)
write.table(countData.norm, outFile, col.names=T, row.names=F, quote=F, sep="\t")


###########################
### sessionInfo and clean and quit
sessionInfo() # R version 3.6.1, EBSeq_1.26.0
date()
rm(list=ls())

###########################
####quit
q("yes")
