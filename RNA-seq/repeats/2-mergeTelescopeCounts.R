# Import and prepare data

BO_1 <- read.table("BO_1/telescope-telescope_report.tsv", header = TRUE, sep = "\t", skip = 1, stringsAsFactors = FALSE)
BO_2 <- read.table("BO_2/telescope-telescope_report.tsv", header = TRUE, sep = "\t", skip = 1, stringsAsFactors = FALSE)
BO_3 <- read.table("BO_3/telescope-telescope_report.tsv", header = TRUE, sep = "\t", skip = 1, stringsAsFactors = FALSE)
IP_1 <- read.table("IP_1/telescope-telescope_report.tsv", header = TRUE, sep = "\t", skip = 1, stringsAsFactors = FALSE)
IP_2 <- read.table("IP_2/telescope-telescope_report.tsv", header = TRUE, sep = "\t", skip = 1, stringsAsFactors = FALSE)
IP_3 <- read.table("IP_3/telescope-telescope_report.tsv", header = TRUE, sep = "\t", skip = 1, stringsAsFactors = FALSE)
TOTAL_1 <- read.table("TOTAL_1/telescope-telescope_report.tsv", header = TRUE, sep = "\t", skip = 1, stringsAsFactors = FALSE)
TOTAL_2 <- read.table("TOTAL_2/telescope-telescope_report.tsv", header = TRUE, sep = "\t", skip = 1, stringsAsFactors = FALSE)
TOTAL_3 <- read.table("TOTAL_3/telescope-telescope_report.tsv", header = TRUE, sep = "\t", skip = 1, stringsAsFactors = FALSE)
INPUT_1 <- read.table("INPUT_1/telescope-telescope_report.tsv", header = TRUE, sep = "\t", skip = 1, stringsAsFactors = FALSE)
INPUT_2 <- read.table("INPUT_2/telescope-telescope_report.tsv", header = TRUE, sep = "\t", skip = 1, stringsAsFactors = FALSE)
INPUT_3 <- read.table("INPUT_3/telescope-telescope_report.tsv", header = TRUE, sep = "\t", skip = 1, stringsAsFactors = FALSE)

# get final counts and rename count columns
BO_1_count <- BO_1[,c(1,3)]
colnames(BO_1_count)[2] <- "BO_1_count"
BO_2_count <- BO_2[,c(1,3)]
colnames(BO_2_count)[2] <- "BO_2_count"
BO_3_count <- BO_3[,c(1,3)]
colnames(BO_3_count)[2] <- "BO_3_count"
IP_1_count <- IP_1[,c(1,3)]
colnames(IP_1_count)[2] <- "IP_1_count"
IP_2_count <- IP_2[,c(1,3)]
colnames(IP_2_count)[2] <- "IP_2_count"
IP_3_count <- IP_3[,c(1,3)]
colnames(IP_3_count)[2] <- "IP_3_count"
TOTAL_1_count <- TOTAL_1[,c(1,3)]
colnames(TOTAL_1_count)[2] <- "TOTAL_1_count"
TOTAL_2_count <- TOTAL_2[,c(1,3)]
colnames(TOTAL_2_count)[2] <- "TOTAL_2_count"
TOTAL_3_count <- TOTAL_3[,c(1,3)]
colnames(TOTAL_3_count)[2] <- "TOTAL_3_count"
INPUT_1_count <- INPUT_1[,c(1,3)]
colnames(INPUT_1_count)[2] <- "INPUT_1_count"
INPUT_2_count <- INPUT_2[,c(1,3)]
colnames(INPUT_2_count)[2] <- "INPUT_2_count"
INPUT_3_count <- INPUT_3[,c(1,3)]
colnames(INPUT_3_count)[2] <- "INPUT_3_count"

# merge BO datasets
mergedBO <- Reduce(function(x,y) merge(x = x, y = y, by = "transcript", all = TRUE), 
                   list(BO_1_count, BO_2_count, BO_3_count))

# merge IP datasets
mergedIP <- Reduce(function(x,y) merge(x = x, y = y, by = "transcript", all = TRUE), 
                   list(IP_1_count, IP_2_count, IP_3_count))

# merge TOTAL datasets
mergedTOTAL <- Reduce(function(x,y) merge(x = x, y = y, by = "transcript", all = TRUE), 
                          list(TOTAL_1_count, TOTAL_2_count, TOTAL_3_count))

# merge INPUT datasets
mergedINPUT <- Reduce(function(x,y) merge(x = x, y = y, by = "transcript", all = TRUE), 
                         list(INPUT_1_count, INPUT_2_count, INPUT_3_count))

# merge all datasets
mergedAll <- Reduce(function(x,y) merge(x = x, y = y, by = "transcript", all = TRUE), 
                    list(mergedBO, mergedIP, mergedTOTAL, mergedINPUT))

# exclude first row ("__no_feature") and substitute NA with 0 (proceed if is.na returns FALSE)
mergedAll = mergedAll[-1,]
rownames(mergedAll) <- 1:nrow(mergedAll)
mergedAll[is.na(mergedAll)] <- 0
any(is.na(mergedAll))

# Add annotation to counts (from a .gtf file)
mergedAll$repName <- sub("_\\d+$", "", mergedAll$transcript) # add repName column to counts data frame

# obtain annotations from .gtf file
rmskAnn <- read.table("mm10.UCSC.repeatsAndVariations.gtf", header=F, sep="\t", stringsAsFactors=F)
rmskAnn <- rmskAnn[, c(11,12,13)]
colnames(rmskAnn) <- c("repName", "repClass", "repFamily")
rmskAnn_dedupl <- rmskAnn[!duplicated(rmskAnn$repName),]

# merge annotations to counts
mergedCountsAnn <- merge(rmskAnn_dedupl, mergedAll, all.x=F, all.y=T, by="repName")

# substitute NA to gene names listed in "repName"
mergedCountsAnn[is.na(mergedCountsAnn$repClass), "repName"] <- NA

# reorder columns
mergedCountsAnn <- mergedCountsAnn[,c(4,1:3,5:16)]

# Write data to a file
write.table(mergedCountsAnn, "GeneCounts_telescope_allSamples_annotated.txt", col.names=T, row.names=F, quote=F, sep="\t")

sessionInfo()
