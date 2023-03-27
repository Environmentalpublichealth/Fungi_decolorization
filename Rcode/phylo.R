setwd("~/Desktop/Jiali/TAMU/Susie/Fungi")

library(reshape2)

fungal_list <- read.csv("list of USDA fungi 7-21-2022.xlsx - Sheet1.csv")
BLAST <- read.csv("ITS/BLAST.result.tsv", header = F, sep = "\t")
Plate_list <- read.csv("plate IDs - Sheet1.csv", header = T)
Plate_list$ITS1 <- paste0("0",Plate_list$row,Plate_list$column,"_",Plate_list$plate)
Plate_list$ITS4 <- paste0("0",Plate_list$row+6,Plate_list$column,"_",Plate_list$plate)

dyeFungi <- merge(fungal_list[,c(2,4,6)], Plate_list, by.x = "Collection.Number", by.y = "Accession", all.x =T)
dyeFungi$No.[duplicated(dyeFungi$Collection.Number)]

length(intersect(fungal_list$Collection.Number, Plate_list$Accession))
length(table(BLAST$V2))

#write.csv(dyeFungi, "Fungal list matching.csv")

# modified the IDs in excel and add genebank IDs for some species
# read the modified list 
newlist <- read.csv("fungal summary.xlsx - isolate list.csv", header = T)
# correct species names
newlist[] <- lapply(newlist, function(x) gsub("Hypsizygous","Hypsizygus", x))
newlist[] <- lapply(newlist, function(x) gsub("Cerrrena","Cerrena", x))
newlist[] <- lapply(newlist, function(x) gsub("B. subcoronatum","Botryobasidium subcoronatum", x))
newlist[] <- lapply(newlist, function(x) gsub("C. subverracisporus|C. subverracisporus 1","Crepidotus subverrucisporus", x))
newlist[] <- lapply(newlist, function(x) gsub("P. chrysosporium","Phanerochaete chrysosporium", x))
newlist[] <- lapply(newlist, function(x) gsub("Schizophyllume commune", "Schizophyllum commune", x))
newlist[] <- lapply(newlist, function(x) gsub("Serpula himantiodes", "Serpula himantioides", x))
newlist[] <- lapply(newlist, function(x) gsub("Sistostrema brinkmanii|Sistostrema brinkmannii ", "Sistostrema brinkmannii", x))
newlist[] <- lapply(newlist, function(x) gsub("Phanerochaete sordida II|Phanerochaete sordida I", "Phanerochaete sordida", x))
# remove NAs, keep one isolate for one species and remove variety names
length(table(newlist$Blast.ID.ITS.)) # 74 
# filter_list <- newlist[!is.na(newlist$ITS1),]
# length(table(filter_list$Blast.ID.ITS.)) # 73 missing Trechispora nivea
# newlist$Blast.ID.ITS.[which(!(newlist$Blast.ID.ITS. %in% filter_list$Blast.ID.ITS.))]
uniqList <- newlist[!duplicated(newlist[,3]),]
# output file and process with biopython to extract sequences in fasta
write.csv(uniqList[,c(2,3,4,8)], "unique_species_isolates.csv", row.names = F)

# add Genebank IDs and order for each species in excel. Seach from NCBI taxa search
uniqList_new <- read.csv("ITS/unique_species_isolates.csv", header = T)
table(uniqList_new$order)

# Run clustal omega and Raxml to build a phylogentic tree, visualize here
library(ggtree)
library(treeio)

alltree <- read.tree("ITS/RAxML_bestTree.GTRGAMMA_a-12345-100-hybrid-avx")
alltree$tip.label <- gsub("_", " ", alltree$tip.label)
annot <- uniqList_new[,c(3,5)]

# edit tip label, add isolate numbers in
species_df <- data.frame(table(newlist$Blast.ID.ITS.))
species_df$Var1 <- gsub('appalachiense.*',"appalachiense", species_df$Var1)
annot <- merge(annot, species_df, by.x = "Blast.ID.ITS.", by.y = "Var1")
names(annot) <- c("label", "order", "Isolates")
newtree <- full_join(alltree, annot, by="label")
newtree@phylo$tip.label <- paste0(newtree@phylo$tip.label, "(", newtree@data$Isolates,")")[1:74]
ggtree(newtree, branch.length = "none", layout = "fan") + 
  geom_tiplab(size =4, fontface=3)+theme_tree()+
  geom_tree(aes(color=order), size =1.5)+
  xlim(0,40)
ggsave("image/RaxmlTree_color_corrected.pdf", height = 11, width = 12)


