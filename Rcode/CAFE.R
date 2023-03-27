setwd("~/Desktop/Jiali/TAMU/Susie/Fungi/")
library(reshape2)

IPS <- read.csv("CAFE/cafeNov/Nov18/allOG_ips.tsv", header = F, sep = "\t")
# filter pfam and superfamily by keeping the smallest Evalue
IPS_interpro <- IPS[-which(IPS$V12 == "-"),c(1,9,12,13)]
IPS_interpro$V9 <- as.numeric(IPS_interpro$V9)
IPS_sort <- IPS_interpro[order(IPS_interpro$V9, decreasing = F),]
IPS_filter <- IPS_sort[!duplicated(IPS_sort[, c(1)]),]
FamilyIDs <- unique(IPS_filter[,c(3,4)])
# sig OGs
OGs <- read.csv("CAFE/cafeNov/Nov18/Significant_families.txt", header = T, sep = "\t")
OGs <- merge(OGs, FamilyIDs, by.x = "X.FamilyID", by.y="V12", all.x = T)

# family change data
change <- read.csv("CAFE/cafeNov/Nov18/IPR_results/Gamma_change.tab", header = T, sep = "\t")
change_large <- read.csv("CAFE/cafeNov/Nov18/IPR_large_results/Gamma_change.tab", header = T, sep = "\t")
change_all <- rbind(change_large, change)
names(change_all)
change_sig <- change_all[which(change_all$FamilyID %in% OGs$X.FamilyID), c(1,2, 3, 4,6,7, 9, 11,13,15, 18, 19, 21, 23)] #c(1,2,3,4,5,6,7,9,13,14,16,17,19,25)
GeneCounts <- read.csv("CAFE/cafeNov/Nov18/IPR.GeneCount.txt", header = T, sep = "\t")
colnames(GeneCounts) <- gsub("\\d+_.*|RSB|PC|_TatD_", "", colnames(GeneCounts))
OGs <- merge(OGs, change_sig, by.x = "X.FamilyID", by.y = "FamilyID")
colnames(OGs) <- gsub("\\d+_.*|RSB|PC|_TatD_", "", colnames(OGs))
write.csv(OGs, "CAFE/cafeNov/Nov18/significant family changes.csv")

# heatmap of gene family expansion and contraction
library(pheatmap)
library(RColorBrewer)
library(viridis)
heatmap_df <- OGs[-5,c(4:17)]
heatmap_df <- GeneCounts[which(GeneCounts$FamilyID %in% OGs$X.FamilyID),]
filter_heatmap <- na.omit(heatmap_df)
#filter_heatmap <- filter_heatmap[!duplicated(filter_heatmap[, c(2)]),]
heatmap_matrix <- as.matrix(filter_heatmap[,3:15])
rownames(heatmap_matrix) <- filter_heatmap$Description
colnames(heatmap_matrix) <- gsub("\\d+_.*", "", colnames(heatmap_matrix))
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(heatmap_matrix, n = 9)
colors <- brewer.pal(n = 9, name = "YlOrRd")
colabel <- c("Trave","Bjead","Cerun","Irplac","Pleos","Fomfom","Pycsa","Trapub","Phlgi","Lenedo", "Echtin", "Schco", "Wolco")
heatmap_matrix <- heatmap_matrix[,colabel]
colnames(heatmap_matrix) <- c("T. versicolor","B. adusta","C. unicolor","I. lacteus","P. ostreatus","F. fomentarius","T. sanguinea", "T. pubescens", "P. gigantea", "L. edodes", "E. tinctorium", "S. commune", "W. cocos")

pdf("image/genefamily_heatmap_Nov18_geneCount.pdf", width = 9, height = 7.5)
#pheatmap(heatmap_matrix, cluster_cols = F, breaks = mat_breaks, color = colors, angle_col = 45)
pheatmap(heatmap_matrix, cluster_cols = F, scale = "row", angle_col = 45)
dev.off()

# process the tree file
library(ggtree)
library(ggplot2)
library(treeio)
# treetext = "(Exigl:247.95[&&NHX:C=+],(((Wolco:88.1488[&&NHX:C=+424/-1119],Pospl:88.1488[&&NHX:C=+297/-1250]):27.7124[&&NHX:C=+65/-104],
# (((Pycsa:53.445[&&NHX:C=+385/-857],Artele:53.445[&&NHX:C=+221/-981]):6.31521[&&NHX:C=+10/-10],
# (Trabet:48.2051[&&NHX:C=+278/-905],Trave:48.2051[&&NHX:C=+427/-928]):11.5551[&&NHX:C=+69/-10]):16.8037[&&NHX:C=+74/-32],
# Fomfom:76.564[&&NHX:C=+543/-793]):39.2973[&&NHX:C=+137/-77]):13.6536[&&NHX:C=+6/-56],
# (Cerun:117.111[&&NHX:C=+400/-876],(Irplac:98.3321[&&NHX:C=+775/-919],(Bjead:81.8982[&&NHX:C=+480/-838],(Phchr:66.7765[&&NHX:C=+359/-827],
# Phlgi:66.7765[&&NHX:C=+255/-966]):15.1217[&&NHX:C=+31/-32]):16.4339[&&NHX:C=+76/-34]):18.779[&&NHX:C=+77/-38]):12.4038[&&NHX:C=+12/-23]):118.435[&&NHX:C=+30/-753]);"

# treetext = "(Pleos:149.013[&&NHX:C=+],((Cerun:82.4256[&&NHX:C=+433/-631],(Irplac:68.7797[&&NHX:C=+782/-730],(Bjead:56.0356[&&NHX:C=+423/-649],
# Phchr:56.0356[&&NHX:C=+346/-648]):12.7441[&&NHX:C=+60/-24]):13.6459[&&NHX:C=+64/-32]):8.91893[&&NHX:C=+16/-12],
# (Pospl:81.6551[&&NHX:C=+366/-1072],(Fomfom:54.2419[&&NHX:C=+550/-577],(Pycsa:42.6456[&&NHX:C=+347/-647],(Artele:37.4043[&&NHX:C=+224/-729],
# (Trave:30.1544[&&NHX:C=+393/-697],Trabet:30.1544[&&NHX:C=+264/-637]):7.24994[&&NHX:C=+54/-5]):5.24131[&&NHX:C=+10/-11]):11.5963[&&NHX:C=+39/-20]):27.4132[&&NHX:C=+144/-52]):9.68946[&&NHX:C=+6/-23]):57.6686[&&NHX:C=+16/-225]);"
treetext <- "(Pleos:149.013[&&NHX:C=+],((Cerun:80.5881[&&NHX:C=+416/-619],(Irplac:66.0893[&&NHX:C=+724/-681],Bjead:66.0893[&&NHX:C=+425/-627]):14.4988[&&NHX:C=+41/-19]):10.4747[&&NHX:C=+9/-9],
(((Trave:18.3734[&&NHX:C=+331/-676],Trapub:18.3734[&&NHX:C=+283/-782]):26.2386[&&NHX:C=+146/-39],Pycsa:44.6119[&&NHX:C=+315/-633]):13.6142[&&NHX:C=+61/-27],
Fomfom:58.2261[&&NHX:C=+521/-560]):32.8367[&&NHX:C=+112/-74]):57.9502[&&NHX:C=+12/-191]);"
tree <- read.nhx(textConnection(treetext))
tree <- read.tree("CAFE/cafeNov/Nov18/SpeciesTree_rooted.txt.ultrametric.all.tre")
head(as.data.frame(as_tibble(tree)))
# orthogroup data
OG_Stats <- read.table("CAFE/cafeNov/Nov18/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv", header = T, sep = "\t")
names(OG_Stats) <- gsub("\\d+_.*|RSB|PC|_TatD_","",names(OG_Stats)) # clean species names
OG_number <- t(OG_Stats[c(8,6),])
colnames(OG_number) <- OG_number[1,] 
OG_number <- data.frame(OG_number[-1,])
OG_number$Number.of.species.specific.orthogroups <- as.numeric(OG_number$Number.of.species.specific.orthogroups)
OG_number$Number.of.orthogroups.containing.species <- as.numeric(OG_number$Number.of.orthogroups.containing.species)
OG_number$Inter.speices <- OG_number$Number.of.orthogroups.containing.species - OG_number$Number.of.species.specific.orthogroups                   
names(OG_number) <- c("Intra-species orthogroups","Total", "Inter-species orthogroups")
OG_number$ID <- rownames(OG_number)
d1 <- reshape2::melt(OG_number[,c(1,3,4)], id = "ID")
level <- c("Trapub", "Trave", "Pycsa", "Fomfom","Bjead","Irplac","Cerun", "Pleos")
level <- gsub("\\d+_.*|RSB|PC|_TatD_","", tree$tip.label)
d1$ID <- factor(d1$ID, levels =level)
p1 <- ggplot(d1, aes(ID, value, label=value, fill = variable)) + geom_bar(stat = "identity") + 
  geom_text(size = 3, position = position_stack(vjust = 0.5))+
  coord_flip() + theme_tree2()+
  scale_fill_manual(values = c("#FFEDA0", "#BBC574"))+
  theme(legend.position='none')+xlab(NULL) + ylab("Number of orthogroups")
  
p1

# change tip label - plot tree
species <- c("Exidia glandulosa", "Schizophyllum commune", "Lentinula edodes","Pleurotus ostreatus","Echinodontium tinctorium", "Cerrena unicolor", "Irpex lacteus","Bjerkandera adusta",
             "Phlebiopsis gigantea","Wolfipora cocos","Fomes fomentarius","Trametes versicolor","Trametes pubescens",
             "Trametes sanguinea")
tree@phylo$tip.label <- species
tree$tip.label <- species
g = ggtree(tree) + geom_tiplab(size =4, fontface=3) + 
  # geom_text(aes(x=branch, label=C), size=2.5, 
  #           vjust=1.3, hjust = 0.6, color="firebrick") + 
  xlim(0,400)+theme_tree2()+
  vexpand(.01, direction = -1)+xlab("Mya")

g
cowplot::plot_grid(g, p1, ncol=2, rel_widths = c(1.2, 1)) 
ggsave("image/genefamilyPhyloNov18.pdf", height = 4, width = 10)
