# clear the workspace
rm(list=ls())

# read the names of the result files
fileList <- list.files(pattern ="ReadsPerGene\\.out\\.tab$")

# Retrieve sample names
sampleNames <-gsub("_mm10.*$","",fileList)

# read the files
myList <- list()
for (i in 1:length(fileList))
{
    message (paste("Reading", fileList[i], sep=": ")) 
    tab <- read.table(fileList[i], header=F, sep="\t", stringsAsFactors=F)
    tab <- tab[,c(1,4)]
    colnames(tab) <- c("geneID", paste("geneCount", sampleNames[i], sep="_"))
    myList[[i]] <- tab
}

# merge the files
mergedTab <- merge(myList[[1]],myList[[2]])
for (i in 3:length(fileList))
{
    mergedTab <- merge(mergedTab,myList[[i]])
}
mergedTab$repName <- sub("_\\d+$", "", mergedTab$geneID)

# add annotation
rmskAnn <- read.table("mm10.UCSC.repeatsAndVariations.gtf", header=F, sep="\t", stringsAsFactors=F)
rmskAnn <- rmskAnn[, c(11,12,13)]
colnames(rmskAnn) <- c("repName", "repClass", "repFamily")
rmskAnn <- rmskAnn[!duplicated(rmskAnn),]
mergedTab_ann <- merge(rmskAnn, mergedTab, all.y=T)
mergedTab_ann[is.na(mergedTab_ann$repClass), "repName"] <- NA
mergedTab_ann <- mergedTab_ann[,c(4,1:3,5:19)]

# write the merged values to a file
write.table(mergedTab_ann, "GeneCounts_all_samples.txt", col.names=T, row.names=F, quote=F, sep="\t")

# sessionInfo and clean and quit
sessionInfo()
date()
rm(list=ls())
q("yes")
