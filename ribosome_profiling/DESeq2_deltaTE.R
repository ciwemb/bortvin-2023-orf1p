# Identify translationally regulated genes

# Load packages

.libPaths("/mnt/sequence/R/3.6.1/DESeq2/")
.libPaths()

library(DESeq2)
library( "dplyr" )
library( "tibble" )

# Set wd and load data
setwd("/mnt/sequence/cdeluca/bortvin-2023-orf1p/")
allCount <- read.table("counts_rna_ribo_annotated.txt", header=T, sep="\t", stringsAsFactors=F)

allCount_sub <- allCount[is.na(allCount$repName),] # subset genes
allCount_sub$mean <- apply(allCount_sub[,-c(1:4)], 1, mean)
allCount_sub <- allCount_sub[allCount_sub$mean!=0, ]
allCount_sub$id <- 1:nrow(allCount_sub)
colnames(allCount_sub) <- sub("geneCount_", "", colnames(allCount_sub))
rownames(allCount_sub) <- allCount_sub$id

counts <- allCount_sub[, c(11:16, 5:10)]
sample_info <- data.frame(sampleID = c("ribo_het_1", "ribo_het_2", "ribo_het_3", "ribo_mael_1", "ribo_mael_2", "ribo_mael_3", "mrna_het_1", "mrna_het_2", "mrna_het_3", "mrna_mael_1", "mrna_mael_2", "mrna_mael_3"), 
                          Condition = factor(c(rep("1",3), rep("2",3), rep("1",3), rep("2",3)), levels = c("1", "2")), 
                          SeqType = factor(rep(c("RIBO", "RNA"), each = 6), levels = c("RNA", "RIBO")))

# Run DESeq2
ddsMat = DESeqDataSetFromMatrix(countData = counts, colData = sample_info, design = ~ Condition + SeqType + Condition:SeqType)
ddsMat = DESeq(ddsMat)

# Obtain fold changes for TE
resultsNames(ddsMat)
res = results(ddsMat, name = "Condition2.SeqTypeRIBO")
summary(res)

# Add annotations
res <- as.data.frame(res)
res$id <- rownames(res)
res2 <- merge(allCount_sub[,c("geneID", "id")], res, by = "id")
res2 <- res2[,-1]

# Sort output by padj
res2.sorted <- res2[order(res2$padj),]

# Write the result
for (i in 2:ncol(res2.sorted)) {
  res2.sorted[,i] <- round(res2.sorted[,i],3)
}
write.table(res2.sorted,"TE_mael_vs_het_DESeq2.txt", sep="\t", quote=F, row.names=F)

######
# Visualize the global transcriptional and translational regulation

# Run DESeq2 for mRNA and RPF counts
rna_counts <- allCount_sub[, c(5:10)]
ribo_counts <- allCount_sub[, c(11:16)]

ddsMat_rna = DESeqDataSetFromMatrix(countData = rna_counts, 
                                    colData = sample_info[which(sample_info$SeqType == "RNA"),], 
                                    design = ~ Condition)

ddsMat_rna = DESeq(ddsMat_rna)
res_rna = results(ddsMat_rna, name = "Condition_2_vs_1")
summary(res_rna)

# Add annotations, sort and write the results
res_rna <- as.data.frame(res_rna)
res_rna$id <- rownames(res_rna)
res_rna_2 <- merge(allCount_sub[,c("geneID", "id")], res_rna, by = "id")
res_rna_2 <- res_rna_2[,-1]
res_rna_2.sorted <- res_rna_2[order(res_rna_2$padj),] # sneak a peak: top downregulated is Mael, OK!
for (i in 2:ncol(res_rna_2.sorted)) {
  res_rna_2.sorted[,i] <- round(res_rna_2.sorted[,i],3)
}
write.table(res_rna_2.sorted, "rna_mael_vs_het_DESeq2.txt", sep="\t", quote=F, row.names=F)

# Repeat for ribo datasets
ddsMat_ribo = DESeqDataSetFromMatrix(countData = ribo_counts, 
                                     colData = sample_info[which(sample_info$SeqType == "RIBO"),], 
                                     design = ~ Condition)

ddsMat_ribo = DESeq(ddsMat_ribo)
res_ribo = results(ddsMat_ribo, name = "Condition_2_vs_1")
summary(res_ribo)

# Add annotations, sort and write the results
res_ribo <- as.data.frame(res_ribo)
res_ribo$id <- rownames(res_ribo)
res_ribo_2 <- merge(allCount_sub[,c("geneID", "id")], res_ribo, by = "id")
res_ribo_2 <- res_ribo_2[,-1]
res_ribo_2.sorted <- res_ribo_2[order(res_ribo_2$padj),]
for (i in 2:ncol(res_ribo_2.sorted)) {
  res_ribo_2.sorted[,i] <- round(res_ribo_2.sorted[,i],3)
}
write.table(res_ribo_2.sorted, "ribo_mael_vs_het_DESeq2.txt", sep="\t", quote=F, row.names=F)

# Plot the mRNA counts vs. the RPFs (with ORF1p-bound mRNAs in evidence)

bound <- read.table("082821-intersection_significant_n1347.txt", header=T, sep="\t", stringsAsFactors=F) # ORF1p-bound, enriched n=1347

par(mai=c(1,2,1,2))
plot(x=res_rna_2[,3], y=res_ribo_2[,3], 
     xlab="RNA-seq log2 fold change",
     ylab="Ribo-seq log2 fold change",
     pch=16,
     col=rgb(128/255, 128/255, 128/255, 0.2),
     xlim = c(-12, 12), ylim = c(-12, 12),
     cex=0.6)
points(x=res_rna_2[res_rna_2$geneID%in%bound$geneID,3], y=res_ribo_2[res_ribo_2$geneID%in%bound$geneID,3], pch=16, col=rgb(1,0,0,0.8), cex=0.6)
abline(a=0, b=1, col="black")
abline(v=c(1,-1), h=c(1,-1), col=c("black", "black", "black", "black"), lty = c(2,2,2,2), lwd = c(1,1,1,1))
legend("topleft", "orf1p-bound", fill = rgb(1,0,0,0.8), bty = 'n', cex = 0.8, y.intersp=2)
legend(x=-12.5, y=11, "n=1313", bty = 'n', cex = 0.8)
legend(x=6.5, y=-10.5, c("dTE n.s."), bty = 'n', cex = 0.8, y.intersp=2)

# Document software
sessionInfo()
