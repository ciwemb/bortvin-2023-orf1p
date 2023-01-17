# Import and prepare data

BO_1 <- read.table("/mnt/sequence/cdeluca/telescope/BO_suc/telescope-telescope_report.tsv", header = TRUE, sep = "\t", skip = 1, stringsAsFactors = FALSE)
BO_2 <- read.table("/mnt/sequence/cdeluca/telescope/BO_suc_2/telescope-telescope_report.tsv", header = TRUE, sep = "\t", skip = 1, stringsAsFactors = FALSE)
BO_3 <- read.table("/mnt/sequence/cdeluca/telescope/BO_suc_3/telescope-telescope_report.tsv", header = TRUE, sep = "\t", skip = 1, stringsAsFactors = FALSE)
IP_1 <- read.table("/mnt/sequence/cdeluca/telescope/IP_suc/telescope-telescope_report.tsv", header = TRUE, sep = "\t", skip = 1, stringsAsFactors = FALSE)
IP_2 <- read.table("/mnt/sequence/cdeluca/telescope/IP_suc_2/telescope-telescope_report.tsv", header = TRUE, sep = "\t", skip = 1, stringsAsFactors = FALSE)
IP_3 <- read.table("/mnt/sequence/cdeluca/telescope/IP_suc_3/telescope-telescope_report.tsv", header = TRUE, sep = "\t", skip = 1, stringsAsFactors = FALSE)
INPUT_suc_1 <- read.table("/mnt/sequence/cdeluca/telescope/INPUT_suc/telescope-telescope_report.tsv", header = TRUE, sep = "\t", skip = 1, stringsAsFactors = FALSE)
INPUT_suc_2 <- read.table("/mnt/sequence/cdeluca/telescope/INPUT_suc_2/telescope-telescope_report.tsv", header = TRUE, sep = "\t", skip = 1, stringsAsFactors = FALSE)
INPUT_suc_3 <- read.table("/mnt/sequence/cdeluca/telescope/INPUT_suc_3/telescope-telescope_report.tsv", header = TRUE, sep = "\t", skip = 1, stringsAsFactors = FALSE)
INPUT_58_1 <- read.table("/mnt/sequence/cdeluca/telescope/INPUT_58/telescope-telescope_report.tsv", header = TRUE, sep = "\t", skip = 1, stringsAsFactors = FALSE)
INPUT_58_2 <- read.table("/mnt/sequence/cdeluca/telescope/INPUT_58_2/telescope-telescope_report.tsv", header = TRUE, sep = "\t", skip = 1, stringsAsFactors = FALSE)
INPUT_58_3 <- read.table("/mnt/sequence/cdeluca/telescope/INPUT_58_3/telescope-telescope_report.tsv", header = TRUE, sep = "\t", skip = 1, stringsAsFactors = FALSE)



# get final count and rename count columns
BO_1_count <- BO_1[,c(1,3)]
colnames(BO_1_count)[2] <- "BO_suc_count"
BO_2_count <- BO_2[,c(1,3)]
colnames(BO_2_count)[2] <- "BO_suc_2_count"
BO_3_count <- BO_3[,c(1,3)]
colnames(BO_3_count)[2] <- "BO_suc_3_count"
IP_1_count <- IP_1[,c(1,3)]
colnames(IP_1_count)[2] <- "IP_suc_count"
IP_2_count <- IP_2[,c(1,3)]
colnames(IP_2_count)[2] <- "IP_suc_2_count"
IP_3_count <- IP_3[,c(1,3)]
colnames(IP_3_count)[2] <- "IP_suc_3_count"
INPUT_suc_1_count <- INPUT_suc_1[,c(1,3)]
colnames(INPUT_suc_1_count)[2] <- "INPUT_suc_count"
INPUT_suc_2_count <- INPUT_suc_2[,c(1,3)]
colnames(INPUT_suc_2_count)[2] <- "INPUT_suc_2_count"
INPUT_suc_3_count <- INPUT_suc_3[,c(1,3)]
colnames(INPUT_suc_3_count)[2] <- "INPUT_suc_3_count"
INPUT_58_1_count <- INPUT_58_1[,c(1,3)]
colnames(INPUT_58_1_count)[2] <- "INPUT_58_count"
INPUT_58_2_count <- INPUT_58_2[,c(1,3)]
colnames(INPUT_58_2_count)[2] <- "INPUT_58_2_count"
INPUT_58_3_count <- INPUT_58_3[,c(1,3)]
colnames(INPUT_58_3_count)[2] <- "INPUT_58_3_count"



# merge BO datasets
mergedBO <- Reduce(function(x,y) merge(x = x, y = y, by = "transcript", all = TRUE), 
                   list(BO_1_count, BO_2_count, BO_3_count))
head(mergedBO)



# merge IP datasets
mergedIP <- Reduce(function(x,y) merge(x = x, y = y, by = "transcript", all = TRUE), 
                   list(IP_1_count, IP_2_count, IP_3_count))
head(mergedIP)



# merge INPUT suc datasets
mergedINPUT_suc <- Reduce(function(x,y) merge(x = x, y = y, by = "transcript", all = TRUE), 
                          list(INPUT_suc_1_count, INPUT_suc_2_count, INPUT_suc_3_count))
head(mergedINPUT_suc)



# merge INPUT 58 datasets
mergedINPUT_58 <- Reduce(function(x,y) merge(x = x, y = y, by = "transcript", all = TRUE), 
                         list(INPUT_58_1_count, INPUT_58_2_count, INPUT_58_3_count))
head(mergedINPUT_58)




# merge all datasets
mergedAll <- Reduce(function(x,y) merge(x = x, y = y, by = "transcript", all = TRUE), 
                    list(mergedBO, mergedIP, mergedINPUT_suc, mergedINPUT_58))
head(mergedAll)
tail(mergedAll)



# exclude first row ("__no_feature") and substitute NA with 0 (proceed if is.na returns FALSE)
mergedAll = mergedAll[-1,]
rownames(mergedAll) <- 1:nrow(mergedAll)
mergedAll[is.na(mergedAll)] <- 0
head(mergedAll)
tail(mergedAll)
any(is.na(mergedAll))


# Add annotation to counts (from a .gtf file)

# add repName column to counts data frame
mergedAll$repName <- sub("_\\d+$", "", mergedAll$transcript)
head(mergedAll)
tail(mergedAll)


# obtain annotations from .gtf file
rmskAnn <- read.table("/mnt/sequence/cdeluca/telescope/mm10.UCSC.repeatsAndVariations.gtf", header=F, sep="\t", stringsAsFactors=F)
rmskAnn <- rmskAnn[, c(11,12,13)]
colnames(rmskAnn) <- c("repName", "repClass", "repFamily")
rmskAnn_dedupl <- rmskAnn[!duplicated(rmskAnn$repName),]
head(rmskAnn_dedupl)


# merge annotations to counts
mergedCountsAnn <- merge(rmskAnn_dedupl, mergedAll, all.x=F, all.y=T, by="repName")
head(mergedCountsAnn)
tail(mergedCountsAnn)



# substitute NA to gene names listed in "repName"
mergedCountsAnn[is.na(mergedCountsAnn$repClass), "repName"] <- NA
tail(mergedCountsAnn)



# reorder columns
mergedCountsAnn <- mergedCountsAnn[,c(4,1:3,5:16)]
head(mergedCountsAnn)


# Write data to a file


write.table(mergedCountsAnn, "GeneCounts_telescope_allSamples_annotated.txt", col.names=T, row.names=F, quote=F, sep="\t")



sessionInfo()
