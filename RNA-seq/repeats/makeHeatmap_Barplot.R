###################################################
# Author: Chiara De Luca
# Date: May 29, 2020

###################################################
### Clear Workspace
rm(list=ls())

###################################################
### load libraries
.libPaths( "/mnt/sequence/R/4.0.4/gplots" )
library(gplots)
library(RColorBrewer) # RColorBrewer_1.1-2
library(dplyr) # dplyr_0.8.3

###################################################
### set wd
setwd("/mnt/sequence/cdeluca/telescope/collapse_count_and_normalize/")
getwd()

###################################################
### read the data and calculate mean of collapsed counts accross triplicate
countData <- read.table("repName_collapsed_count_upper_quantile_normalized_all_samples_telescope.txt", header=T, sep="\t", stringsAsFactors=F)
countData$mean_BO_suc <- apply( countData[,2:4], 1, mean )
countData$mean_IP_suc <- apply( countData[,5:7], 1, mean )
countData$mean_INPUT_suc <- apply( countData[,8:10], 1, mean )
countData$mean_INPUT_58 <- apply( countData[,11:13], 1, mean )
for (i in 14:ncol(countData)) {
        countData[,i] <- round(countData[,i],2)
        }

###################################################
### set threshhold (based on known repressed repeats i.e. Lx)
filter( countData, repName == "Lx" ) # mean_INPUT_suc = 7085.47
expr.countData <- filter(countData, mean_INPUT_suc>7000)

###################################################
### add annotations to exclude "Low_complexity", "rRNA", "Satellite", "Simple_repeat", "tRNA" repeats
rmskAnn <- read.table("/mnt/sequence/cdeluca/telescope/mm10.UCSC.repeatsAndVariations.gtf", header=F, sep="\t", stringsAsFactors=F)
rmskAnn <- rmskAnn[, c(11,12,13)]
colnames(rmskAnn) <- c("repName", "repClass", "repFamily")
rmskAnn_dedupl <- rmskAnn[!duplicated(rmskAnn$repName),]
head(rmskAnn_dedupl)

expr.ann.countData <- merge(rmskAnn_dedupl, expr.countData, all.x=F, all.y=T, by="repName")
head(expr.ann.countData)
expr.ann.countData <- expr.ann.countData[!expr.ann.countData$repFamily%in%c("Low_complexity", "rRNA", "Satellite", "Simple_repeat", "tRNA"),]
expr.ann.countData <- expr.ann.countData[!expr.ann.countData$repName%in%"ID4_",]

rownames(expr.ann.countData) <- 1:nrow(expr.ann.countData)

# Convert data to a matrix and add repName annotation
countData_matrix_mean <- data.matrix(log2(expr.ann.countData[,-c(1:15)]))
rownames(countData_matrix_mean) <- expr.ann.countData[,1]

####################################################
####### optional: to make a matrix of a subset of repeats only 
sub.expr.ann.CountData <- expr.ann.countData[expr.ann.countData$repName%in%c("B2_Mm1a", "B2_Mm1t", "MMERVK10C-int", "L1Md_T", "L1Md_A", "L1Md_F2"),]

rownames(sub.expr.ann.CountData) <- sub.expr.ann.CountData[,1]

countData_matrix_sub <- data.matrix(log2(sub.expr.ann.CountData[,-c(1:3,16:19)]))

countData_matrix_sub_nolog <- data.matrix(sub.expr.ann.CountData[,-c(1:3,16:19)])

#####################################################
### Customizing and plotting the heat maps
### Subset of repeats:
my_palette <- colorRampPalette(c("white", "gold", "red"))(n = 299)

#### color code triplicate samples on the side (optional)
myColumns <- data.frame(RNASample=colnames(countData_matrix_sub))
myColumns$RNASample2 <- sub("_.$", "", colnames(countData_matrix_sub))
myColSideColors <- rainbow(length(unique(sub("_.$", "", colnames(countData_matrix_sub)))))
myColSideColors_mat <- data.frame(RNASample2=unique(sub("_.$", "", colnames(countData_matrix_sub))), color=myColSideColors, stringsAsFactors=FALSE)
myColSideColors_mat <- merge(myColumns,myColSideColors_mat)

#### create heatmap (subset)
pdf("/mnt/sequence/cdeluca/telescope/collapse_count_and_normalize/repName_telescope_heatmap_subset.pdf", width = 15, height = 15)
heatmap.2(countData_matrix_sub,
          main = "log2 normalized collapsed counts, all samples",  # heat map title
          notecol="black",  # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",  # turns off trace lines inside the heat map
          margins =c(16,12),   # margins around plot
          col=my_palette,   # use on color palette defined earlier
          ColSideColors = myColSideColors_mat$color,
          dendrogram= "row",
          Colv = c(1:3,10:12,4:6,7:9)) 

dev.off()

########################################################
### All repeats:
### creates a customized color palette from red to green
my_palette <- colorRampPalette(c("blue", "white", "red", "black"))(n = 299)

########### optional (to color code triplicate samples)
myColumns <- data.frame(RNASample=colnames(countData_matrix_sub))
myColumns$RNASample2 <- sub("_.$", "", colnames(countData_matrix_sub))
myColSideColors <- rainbow(length(unique(sub("_.$", "", colnames(countData_matrix_sub)))))
myColSideColors_mat <- data.frame(RNASample2=unique(sub("_.$", "", colnames(countData_matrix_sub))), color=myColSideColors, stringsAsFactors=FALSE)
myColSideColors_mat <- merge(myColumns,myColSideColors_mat)

# create heatmap (mean counts)
pdf("/mnt/sequence/cdeluca/telescope/collapse_count_and_normalize/repName_telescope_heatmap.pdf", width = 10, height = 20)
heatmap.2(countData_matrix_mean,
          main = "log2 normalized collapsed counts",  # heat map title
          notecol="black",  # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",  # turns off trace lines inside the heat map
          margins =c(16,12),   # margins around plot
          col=my_palette,   # use on color palette defined earlier
          )

dev.off()

# optional conditions for heatmap: 
ColSideColors = myColSideColors_mat$color,
Colv = c(1:3,10:12,4:6,9,7,8)) 

########################
### white - gold - red - black
my_palette <- colorRampPalette(c("white", "gold", "red", "black"))(n = 299)

pdf("/mnt/sequence/cdeluca/telescope/collapse_count_and_normalize/repName_telescope_heatmap_wgrb.pdf", width = 10, height = 20)
heatmap.2(countData_matrix_mean,
          main = "log2 normalized collapsed counts",  # heat map title
          notecol="black",  # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",  # turns off trace lines inside the heat map
          margins =c(16,12),   # margins around plot
          col=my_palette,   # use on color palette defined earlier
)

dev.off()

######################################################
### Customizing and plotting the bar plot
myPalette2 <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]

pdf("/mnt/sequence/cdeluca/chiara_de_luca/plots/041220_RepName_normalized_collapsed_counts_sub_barplot.pdf",width = 55, height = 22.5)
barplot(countData_matrix_sub_nolog, main = "normalized collapsed counts", col = myPalette2[1:100], cex.axis = 2, cex.main = 2, cex.names = 2)
legend("topright", rownames(countData_matrix_sub_nolog), fill = myPalette2[1:100])
dev.off()


###########################
### sessionInfo and clean and quit
sessionInfo() # R version 3.6.1 
date()
rm(list=ls())

