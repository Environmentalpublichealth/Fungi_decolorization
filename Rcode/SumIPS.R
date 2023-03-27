setwd("~/Desktop/Jiali/TAMU/Susie/Fungi/")
library(reshape2)

IPS <- read.csv("CAFE/cafeNov/Nov18/allOG_ips.tsv", header = F, sep = "\t")
# filter pfam and superfamily by keeping the smallest Evalue with ips IDs
IPS_interpro <- IPS[-which(IPS$V12 == "-"),c(1,9,12,13)]
IPS_interpro$V9 <- as.numeric(IPS_interpro$V9)
IPS_sort <- IPS_interpro[order(IPS_interpro$V9, decreasing = F),]
IPS_filter <- IPS_sort[!duplicated(IPS_sort[, c(1)]),]
IPS_filter$FamilyID <- paste0("OG",sprintf("%07d", (IPS_filter$V1-1)))
# sig OGs
OG_counts <- read.csv("CAFE/cafeNov/Nov18/Orthogroups.GeneCount.tsv", header = T, sep = "\t")
colnames(OG_counts) <- gsub("\\d+_.*|RSB|PC|_TatD_", "", colnames(OG_counts))
OG_counts <- merge(OG_counts, IPS_filter[,c(5,3,4)], by.x = "Orthogroup", by.y="FamilyID", all.x = T)
OGcount_filter <- na.omit(OG_counts)
FamilyIDs <- unique(OGcount_filter[,c(17,18)])

Familydf <- data.frame(matrix(ncol = 16, nrow = 3043))
x <- c("FamilyID","Description", colnames(OGcount_filter)[2:15])
colnames(Familydf) <- x
n=1
for (ID in FamilyIDs$V12) {
  oneFamily <- OGcount_filter[OGcount_filter$V12 ==ID,-c(1)]
  totalCounts <- colSums(oneFamily[,1:14])
  Familydf[n,] <- c(FamilyIDs[n,], totalCounts) 
  n = n+1
}
write.csv(Familydf, "CAFE/cafeNov/Nov18/FamilySum14species.csv", row.names = F)

# CAzymes & redox
WDGenes <- read.csv("CAFE/cafeNov/Nov18/Wood decay gene families.csv")
blast <- read.csv("CAFE/cafeNov/Nov18/OG_CAZymes.blast", header = F, sep = "\t")
blast$OG <- paste0("OG",sprintf("%07d", (blast$V1-1)))
IPS$OG <- paste0("OG",sprintf("%07d", (IPS$V1-1)))
selectOG <- paste0("OG",sprintf("%07d", as.integer(WDGenes[1, -c(1,2)]-1)))
WDGdf <- data.frame(matrix(ncol = 15, nrow = 33))
x <- c("FamilyID", colnames(OG_counts)[2:15])
colnames(WDGdf) <- x
for (i in c(1:33)) {
  selectOG <- paste0("OG",sprintf("%07d", as.integer(WDGenes[i, -c(1,2)]-1)))
  oneFamily <- OG_counts[OG_counts$Orthogroup %in% selectOG,-c(1)]
  totalCounts <- colSums(oneFamily[,1:14])
  WDGdf[i,] <- c(WDGenes[i,1], totalCounts) 
}
write.csv(WDGdf, "wood decay gene counts.csv", row.names = F)
