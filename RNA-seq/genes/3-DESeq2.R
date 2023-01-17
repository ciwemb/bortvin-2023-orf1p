# Load packages

.libPaths("/mnt/sequence/R/3.6.1/DESeq2/")
.libPaths()

library(DESeq2)
library( "dplyr" )
library( "tibble" )

# Load data
allCount <- read.table("GeneCounts_all_samples.txt", header=T, sep="\t", stringsAsFactors=F)

allCount_sub <- allCount[!allCount$repFamily%in%c("Low_complexity", "rRNA", "Satellite", "Simple_repeat", "tRNA"),]
allCount_sub$mean <- apply(allCount_sub[,-c(1:4)], 1, mean)
allCount_sub <- allCount_sub[allCount_sub$mean!=0, ]
allCount_sub$id <- 1:nrow(allCount_sub)
colnames(allCount_sub) <- sub("geneCount_", "", colnames(allCount_sub))
rownames(allCount_sub) <- allCount_sub$id

# Run DESeq2
cnts <- allCount_sub[,c("BO_1", "BO_2", "BO_3", "IP_1", "IP_2", "IP_3")]
cond <- factor(c(rep("BO_comb",3),rep("IP_comb",3)), levels=c("BO_comb", "IP_comb"))
dds <- DESeqDataSetFromMatrix(cnts, DataFrame(cond), ~ cond)

vsd <- vst( dds )
plotPCA( vsd, intgroup="cond" )

dds <- DESeq( dds )

sizeFactors( dds )

res <- results(dds)

# calculate mean of normalized counts
df_ncounts <- data.frame( counts( dds, normalized=TRUE ) )
df_ncounts$BO_mean <- apply( df_ncounts[,1:3], 1, mean )
df_ncounts$IP_mean <- apply( df_ncounts[,4:6], 1, mean )

# prepare the output
res2 <- as.data.frame(res@listData)
res2 <- cbind( res2, df_ncounts[,c("BO_mean","IP_mean")] ) # add mean of normalized counts columns
res2$id <- res@rownames
res2 <- res2[,c(9,1:8)]
for (i in 2:ncol(res2)) {
  res2[,i] <- round(res2[,i],3)
}

# add annotation
res2 <- merge(allCount_sub[,c("geneID", "repName", "repClass", "repFamily", "id")], res2)
res2 <- res2[,-1]

# sort output by padj & absolute value of log2FoldChange
res2.sorted <- res2[order(abs(res2$log2FoldChange), decreasing=T),]
res2.sorted <- res2.sorted[order(abs(res2.sorted$padj)),]
colnames(res2.sorted)[6] <- paste("log2FoldChange_", levels(cond)[2], "_vs_", levels(cond)[1], sep="")

# sanity check
filter( res2.sorted, repName == "L1Md_T" )

# write the result
write.table(res2.sorted,file=paste(levels(cond)[2], "_vs_", levels(cond)[1], "_DESeq2.txt", sep=""), sep="\t", quote=F, row.names=F)

# Document software
sessionInfo()
