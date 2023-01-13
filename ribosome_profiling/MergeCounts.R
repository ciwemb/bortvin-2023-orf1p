# Clean space
rm(list=ls())

# Import data
mrna_het_1 <- read.table("mrna_het_1_exon.count.txt", skip = 1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
mrna_het_2 <- read.table("mrna_het_2_exon.count.txt", skip = 1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
mrna_het_3 <- read.table("mrna_het_3_exon.count.txt", skip = 1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
mrna_mael_1 <- read.table("mrna_mael_1_exon.count.txt", skip = 1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
mrna_mael_2 <- read.table("mrna_mael_2_exon.count.txt", skip = 1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
mrna_mael_3 <- read.table("mrna_mael_3_exon.count.txt", skip = 1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

ribo_het_1 <- read.table("ribo_het_1_CDS.count.txt", skip = 1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ribo_het_2 <- read.table("ribo_het_2_CDS.count.txt", skip = 1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ribo_het_3 <- read.table("ribo_het_3_CDS.count.txt", skip = 1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ribo_mael_1 <- read.table("ribo_mael_1_CDS.count.txt", skip = 1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ribo_mael_2 <- read.table("ribo_mael_2_CDS.count.txt", skip = 1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ribo_mael_3 <- read.table("ribo_mael_3_CDS.count.txt", skip = 1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Get read counts and rename count columns
mrna_het_1.count <- mrna_het_1[,c(1,7)]
colnames(mrna_het_1.count) <- c("geneID", "geneCount_mrna_het_1")
mrna_het_2.count <- mrna_het_2[,c(1,7)]
colnames(mrna_het_2.count) <- c("geneID", "geneCount_mrna_het_2")
mrna_het_3.count <- mrna_het_3[,c(1,7)]
colnames(mrna_het_3.count) <- c("geneID", "geneCount_mrna_het_3")

mrna_mael_1.count <- mrna_mael_1[,c(1,7)]
colnames(mrna_mael_1.count) <- c("geneID", "geneCount_mrna_mael_1")
mrna_mael_2.count <- mrna_mael_2[,c(1,7)]
colnames(mrna_mael_2.count) <- c("geneID", "geneCount_mrna_mael_2")
mrna_mael_3.count <- mrna_mael_3[,c(1,7)]
colnames(mrna_mael_3.count) <- c("geneID", "geneCount_mrna_mael_3")

ribo_het_1.count <- ribo_het_1[,c(1,7)]
colnames(ribo_het_1.count) <- c("geneID", "geneCount_ribo_het_1")
ribo_het_2.count <- ribo_het_2[,c(1,7)]
colnames(ribo_het_2.count) <- c("geneID", "geneCount_ribo_het_2")
ribo_het_3.count <- ribo_het_3[,c(1,7)]
colnames(ribo_het_3.count) <- c("geneID", "geneCount_ribo_het_3")

ribo_mael_1.count <- ribo_mael_1[,c(1,7)]
colnames(ribo_mael_1.count) <-  c("geneID", "geneCount_ribo_mael_1")
ribo_mael_2.count <- ribo_mael_2[,c(1,7)]
colnames(ribo_mael_2.count) <- c("geneID", "geneCount_ribo_mael_2")
ribo_mael_3.count <- ribo_mael_3[,c(1,7)]
colnames(ribo_mael_3.count) <- c("geneID", "geneCount_ribo_mael_3")

# Merge mrna_het datasets
merged.mrna_het <- Reduce(function(x,y) merge(x = x, y = y, by = "geneID", all = TRUE), 
                          list(mrna_het_1.count, mrna_het_2.count, mrna_het_3.count))

# Merge mrna_mael datasets
merged.mrna_mael <- Reduce(function(x,y) merge(x = x, y = y, by = "geneID", all = TRUE), 
                           list(mrna_mael_1.count, mrna_mael_2.count, mrna_mael_3.count))

# Merge ribo_het datasets
merged.ribo_het <- Reduce(function(x,y) merge(x = x, y = y, by = "geneID", all = TRUE), 
                          list(ribo_het_1.count, ribo_het_2.count, ribo_het_3.count))

# Merge ribo_mael datasets
merged.ribo_mael <- Reduce(function(x,y) merge(x = x, y = y, by = "geneID", all = TRUE), 
                           list(ribo_mael_1.count, ribo_mael_2.count, ribo_mael_3.count))

# Merge all datasets
mergedAll <- Reduce(function(x,y) merge(x = x, y = y, by = "geneID", all = TRUE), 
                    list(merged.mrna_het, merged.mrna_mael, merged.ribo_het, merged.ribo_mael))

any(is.na(mergedAll)) # returns TRUE, mostly ncRNAs & Rik genes that are present in mrna datasets but are missing in ribo datasets

mergedAll <- mergedAll[!is.na(mergedAll$geneCount_ribo_het_1),] # exclude rows containing NAs in geneCount_ribo_het_1 column
any(is.na(mergedAll)) # now it returns FALSE

######
# Add annotation to counts (from a .gtf file)
mergedAll$repName <- sub("_\\d+$", "", mergedAll$geneID) # add repName column to counts data frame

# obtain annotations from .gtf file
rmskAnn <- read.table("mm10.UCSC.repeatsAndVariations.gtf", header=F, sep="\t", stringsAsFactors=F)
rmskAnn <- rmskAnn[, c(11,12,13)]
colnames(rmskAnn) <- c("repName", "repClass", "repFamily")
rmskAnn_dedupl <- rmskAnn[!duplicated(rmskAnn$repName),]

# merge annotations to counts
mergedCountsAnn <- merge(rmskAnn_dedupl, mergedAll, all.x=F, all.y=T, by="repName")

# substitute NA to gene names listed in "repName"
mergedCountsAnn[is.na(mergedCountsAnn$repClass), "repName"] <- NA
tail(mergedCountsAnn)

# reorder columns
mergedCountsAnn <- mergedCountsAnn[,c(4,1:3,5:16)]
head(mergedCountsAnn)

# Write data to a file
write.table(mergedCountsAnn, "counts_rna_ribo_annotated.txt", col.names=T, row.names=F, quote=F, sep="\t")

######
sessionInfo()
