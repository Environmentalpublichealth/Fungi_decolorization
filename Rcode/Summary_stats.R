setwd("~/Desktop/Jiali/TAMU/Susie/Fungi/")

library(readxl)
library(ggplot2)
library(ggpmisc)
library(hrbrthemes)

fileName <- "PDA result of USDA fungus 09262022.xlsx"
excel_sheets(fileName)
Screening<- read_excel(fileName, sheet = 1)
# convert contaminant to NA
Screening[Screening == 'contaminated'] <- NA
Screening[Screening == 'contaminate'] <- NA
Screening[Screening == 'N/A'] <- NA

# reformat the whole table
library(reshape2)
tbl_melt <- melt(Screening, id.vars = c("No.","batch","Collection Number","Blast ID(ITS)","Scientific Name", "Rot"))
# add factors
days <- vector() # culture days
for (i in 2:14) {
  day <- c(rep(i,320*9))  
  days <- c(days,day)
}
tbl_melt$Day <- days
treatment <- c(rep(c(rep("PDA",320*3),rep("PDA+Dye",320*3),rep("Dye",320*3)), 13))
tbl_melt$media <- treatment # media types
replicates <- rep(c(rep(1,320),rep(2,320),rep(3,320)),3*13)
tbl_melt$replicate <- replicates # experimental replicate ID

# analyze dye decolor
# summarize decolor rate
library(Rmisc)
dyeTable <- tbl_melt[which(tbl_melt$media == "Dye"),]
dyeTable$value <- as.numeric(dyeTable$value)
dyeColor <- summarySE(dyeTable, measurevar="value", groupvars=c("No.","Day"), na.rm = TRUE)
# remove days 11-14
dyeColor <- dyeColor[which(dyeColor$Day < 11),]
# plot heatmap
ggplot(dyeColor, aes(x= Day, y=No., fill= value)) + 
  geom_tile()+
  scale_fill_gradient(low="red", high="white")+
  ylab("")
# The plot couldn't fit in so many species, let's plot those with decolor ability
goodSpecies <- dyeColor$No.[which(dyeColor$Day==10 & dyeColor$value > 1.3)]
goodData <- dyeColor[dyeColor$No. %in% goodSpecies, ]

badData <- dyeColor[!(dyeColor$No. %in% goodSpecies), ]
badData_reshape <- reshape(badData[,c(1,2)], idvar = "No.",direction = "wide", timevar = "Day")
badData_reshape_species <- merge(badData_reshape, Screening[,c(1,2,4,5)], by = "No.", all.x = T)
library(pheatmap)
# convert data into matrix
goodData_reshape <- reshape(goodData[,c(1,2,4)], idvar = "No.",direction = "wide", timevar = "Day")
goodData_reshape_species <- merge(goodData_reshape, Screening[,c(1,2,4,5)], by = "No.", all.x = T)
rownames(goodData_reshape) <- sub(".","",goodData_reshape$No.)
goodData_heatmap <- as.matrix(goodData_reshape[,c(2:10)])
colnames(goodData_heatmap) <- c(2:10)
pdf("dyeHeatmap_good0125.pdf", width = 4, height = 7)
pheatmap(goodData_heatmap, cluster_cols = F, color=colorRampPalette(c("#B53330", "white"))(50))
dev.off()

# Summarize top decomposer species
dyeColor_select <- dyeColor_reshape_species[dyeColor_reshape_species$`Blast ID(ITS)` == "Trametes sanguinea",]
table(dyeColor_select$value.10)

# plot all dye decolor
dyeColor_reshape <- reshape(dyeColor[,c(1,2,4)], idvar = "No.",direction = "wide", timevar = "Day")
rownames(dyeColor_reshape) <- dyeColor_reshape$No.
dyeColor_heatmap <- as.matrix(dyeColor_reshape[,c(2:10)])
colnames(dyeColor_heatmap) <- c(2:10)
# remove rows with na
dyeColor_heatmap <- na.omit(dyeColor_heatmap)
pdf("dyeHeatmap_all_label.pdf", width = 4, height = 30)
pheatmap(dyeColor_heatmap, cluster_cols = F, show_rownames = T, color=colorRampPalette(c("#B53330", "white"))(50))
dev.off()

dyeColor_reshape_species$value.10 <- ceiling(dyeColor_reshape_species$value.10)
dyeColor_reshape_species <- merge(dyeColor_reshape, Screening[,c(1,2,4,5)], by = "No.", all.x = T)
write.csv(dyeColor_reshape_species, "image/decolor level all isolates.csv", row.names = F)

# plot growth data
growth <- tbl_melt[which(tbl_melt$media == "PDA+Dye"),]
growth$value <- as.numeric(growth$value)
growth_sum <- summarySE(growth, measurevar="value", groupvars=c("No.","Day"), na.rm = TRUE)
# remove day 11-14
growth_sum <- growth_sum[which(growth_sum$Day < 11),]
# reshape into heatmap matrix
growth_reshape <- reshape(growth_sum[,c(1,2,4)], idvar="No.",direction = "wide", timevar = "Day")
growth_good <- growth_reshape[growth_reshape$No. %in% goodSpecies,]
rownames(growth_good) <- sub(".","",growth_good$No.)
growth_heatmap <- as.matrix(growth_good[,c(2:10)])
heatmap_row <- as.dendrogram(hclust(dist(goodData_heatmap)))
colnames(growth_heatmap) <- c(2:10)
growth_heatmap <- growth_heatmap[order.dendrogram(heatmap_row),]
pdf("growthHeatmap_good0125.pdf", width = 3.2, height = 7)
pheatmap(growth_heatmap, cluster_cols = F, cluster_rows = F, show_rownames = T, color=colorRampPalette(c("white", "#d37243"))(50))
dev.off()

# calculate correlation between growth and dye decolor
corr_tbl <- data.frame(matrix(ncol=3, nrow = 0))
colnames(corr_tbl) <- c("ID","cor", "pvalue")
for (ID in goodSpecies) {
  one <- tbl_melt[which(tbl_melt$No. == ID & tbl_melt$Day < 11),]
  grow <- one[which(one$media == "PDA+Dye"),]
  grow$value <- as.numeric(grow$value)
  dye <- one[which(one$media == "Dye"),]
  dye$value <- as.numeric(dye$value)
  res <- cor.test(grow$value, dye$value, 
                method = "pearson")
  cc <- res$estimate
  pvalue <- res$p.value
  corr_tbl[nrow(corr_tbl)+1, ] <- c(ID, cc, pvalue)
}

chisquare_data <- merge(growth_reshape[,c(1,10)], dyeColor_reshape[,c(1,10)], by = "No.")
table(ceiling(chisquare_data$value.10.x), ceiling(chisquare_data$value.10.y))
chisq.test(ceiling(chisquare_data$value.10.x), ceiling(chisquare_data$value.10.y), correct = F)

# plot corr
corr_tbl <- corr_tbl[order(corr_tbl$cor, decreasing = T),]
corr_tbl$ID <- factor(corr_tbl$ID, levels = corr_tbl$ID)
corr_tbl$cor <- as.numeric(corr_tbl$cor)
corr_tbl$pvalue <- as.numeric(corr_tbl$pvalue)
corr_tbl$ID <- sub(".", "", corr_tbl$ID)
ggplot(data = corr_tbl, aes(x = ID, y = cor, fill = pvalue))+
  geom_bar(stat="identity")+
  scale_fill_gradient2(low = "white",mid = "brown", high = "gray", midpoint = .01)+
  theme_bw(base_size = 12)+
  theme(axis.text.x = element_text(colour="black", angle = 45, vjust = 0.8, hjust=0.8))+
  labs(x= "", y = "correlation", fill = "p-value")
ggsave("corr_plot_new.pdf", width = 11, height = 3)  

# species and types
goodData_reshape_species <- merge(goodData_reshape, Screening[,c(1,2,4,5)], by = "No.", all.x = T)
# dye barplot
library(dplyr)
data_sum <- goodData_reshape_species %>% group_by(Rot) %>%
  summarise(mean=mean(value.10),
            median=median(value.10),
            se=(sd(value.10))/sqrt(length(value.10)))
ggplot()+geom_bar(data=data_sum, aes(y = mean, x = Rot), fill="brown",
                  stat = "identity", width = 0.4)+
  geom_errorbar(data=data_sum, aes(y = mean, x = Rot,
                                   ymin = mean - se,
                                   ymax = mean + se),
                stat = "identity", width = 0.1, size =1)+
  geom_point(data = goodData_reshape_species, aes(y = value.10, x = Rot),
             size = 1.4, alpha = 0.5,
             position = position_jitter(width = 0.2, height = 0.1))+
  theme_classic(base_size = 14)+
  labs(x = "Rot type", y = "Dye decoloration")
ggsave("dye barplot.pdf", height = 3.5, width = 4)
# Welch's test on brown and white
library(rstatix)
stat.test <- goodData_reshape_species %>%
  t_test(value.10 ~ Rot) %>%
  add_significance()
stat.test

data_sum <- goodData_reshape_species %>% group_by(`Blast ID(ITS)`) %>%
  summarise(mean=mean(value.10),
            median=median(value.10),
            se=(sd(value.10))/sqrt(length(value.10)))
data_sum$`Blast ID(ITS)` <- factor(data_sum$`Blast ID(ITS)`, levels = data_sum$`Blast ID(ITS)`[order(data_sum$mean, decreasing = T)])
ggplot()+geom_bar(data=data_sum, aes(y = mean, x = `Blast ID(ITS)`), fill="brown",
                  stat = "identity", width = 0.7)+
  geom_errorbar(data=data_sum, aes(y = mean, x = `Blast ID(ITS)`,
                                   ymin = mean - se,
                                   ymax = mean + se),
                stat = "identity", width = 0.4, size =0.6)+
  geom_point(data = goodData_reshape_species, aes(y = value.10, x = `Blast ID(ITS)`),
             size = 1.4, alpha = 0.5,
             position = position_jitter(width = 0.2, height = 0.1))+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(colour="black", angle = 60, vjust = 1, hjust=0.95))+
  labs(x = "Species", y = "Dye decoloration")
ggsave("dye barplot species.pdf", height = 4.2, width = 13)


## Difference index
filtered_dyeColor <- na.omit(dyeColor_reshape_species)
filtered_dyeColor$`Blast ID(ITS)` <- gsub("Hypsizygous ulmarius", "Hypsizygus ulmarius", filtered_dyeColor$`Blast ID(ITS)`)
mean(filtered_dyeColor$value.10)
sd(filtered_dyeColor$value.10)
species <- unique(filtered_dyeColor$`Blast ID(ITS)`)
# create empty dataframe
df <- data.frame(matrix(ncol = 4, nrow = 74))
x <- c("Taxa", "SMD", "SD", "p-value")
colnames(df) <- x
n=1
for (ID in species) {
sample <- filtered_dyeColor[filtered_dyeColor$`Blast ID(ITS)` == ID,]
mean(sample$value.10)
SD <- sd(sample$value.10)
SMD <- (mean(sample$value.10)-mean(filtered_dyeColor$value.10))/sqrt((sd(filtered_dyeColor$value.10)+sd(sample$value.10))/2) 
# calculate p-value
test <- wilcox.test(filtered_dyeColor$value.10, sample$value.10)
pvalue <- test$p.value
df[n,] <- c(ID, SMD, SD, pvalue)
n=n+1
}
write.csv(df, "dyeVariations.csv")
# volcano plot
library(ggrepel)
df$`p-value` <- as.numeric(df$`p-value`)
df$logP <- -log10(df$`p-value`)
df$SD <- as.numeric(df$SD)
df$SMD <- as.numeric(df$SMD)
df$Taxa <- factor(df$Taxa, levels = df$Taxa[order(df$`p-value`, decreasing = F)])
ggplot(data = subset(df, `p-value`<0.05), aes(x = Taxa, y = SMD, fill = `p-value`))+
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=SMD-SD, ymax=SMD+SD), width=.2,
                position=position_dodge(.9))+ 
  #geom_text_repel(
  #   data = subset(df, `p-value`<0.05),
  #   aes(label = Taxa),
  #   size = 2,
  #   box.padding = unit(0.3, "lines"),
  #   point.padding = unit(0.3, "lines")
  # )+
  theme_classic()+
  theme(axis.text.x = element_text(colour="black", angle = 60, vjust = 1, hjust=0.95))+
  labs(y = paste("SMD", "SD", sep = "\u00B1"), fill="p-value")
ggsave("image/species variation.pdf", height = 4.5, width = 6)
df$SMD <- as.numeric(df$SMD)
df$SD <- as.numeric(df$SD)
pc <- prcomp(na.omit(df)[,2:3],
             center = TRUE,
             scale. = TRUE)
summary(pc)
library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)
g <- ggbiplot(pc,
              obs.scale = 1,
              var.scale = 1,
              groups = filtered_dyeColor$`Blast ID(ITS)`,
              ellipse = TRUE,
              circle = TRUE,
              ellipse.prob = 0.68)
plot(pc$x[,1], pc$x[,2])

h <- hclust(dist(dyeColor_heatmap))
h$labels[h$order]
group1 <- filtered_dyeColor[h$order[1:241],]
group2 <- filtered_dyeColor[h$order[242:304],]

# output group2
group2$No. <- sub(".","",group2$No.)
write.csv(group2, "hcluster_goodDye.csv", quote = F, row.names = F)
# scoring
filter_growth$score <- filter_growth$value.10.x + filter_growth$value.10.y
df_scores <- data.frame(table(ceiling(filter_growth$score), filter_growth$`Blast ID(ITS)`))
df_scores_positive <- df_scores[-(which(df_scores$Freq == 0)),]
ggplot(data = df_scores_positive, aes(x = Freq, y = Var1))+
  geom_point(position = position_jitter(width = 0.2, height = 0.1))

filter_growth$index <- c(rep(0, 304))  
filter_growth$index[which(filter_growth$value.10.y < 3 & filter_growth$value.10.x == 0)] = 1
filter_growth$index[which(filter_growth$value.10.y >= 3 & filter_growth$value.10.x == 0)] = 2
filter_growth$index[which(filter_growth$value.10.y < 3 & filter_growth$value.10.y>0 &  filter_growth$value.10.x < 3)] = 3
filter_growth$index[which(filter_growth$value.10.y >= 3 & filter_growth$value.10.x < 3 & filter_growth$value.10.x > 0)] = 4
filter_growth$index[which(filter_growth$value.10.y < 3 & filter_growth$value.10.y > 0  & filter_growth$value.10.x >= 3)] = 5
filter_growth$index[which(filter_growth$value.10.y >= 3 & filter_growth$value.10.x >= 3)] = 6
df_scores <- data.frame(table(ceiling(filter_growth$index), filter_growth$`Blast ID(ITS)`))
df_scores$ptg <- df_scores$Freq/rep(table(filter_growth$`Blast ID(ITS)`), each = 6) * 100
df_scores_positive <- df_scores[-(which(df_scores$Freq == 0)),]
ggplot(data = df_scores_positive, aes(x = Var1, y = ptg))+
  geom_point(position = position_jitter(width = 0.2, height = 0.1))+
  geom_text_repel(
    data = subset(df_scores_positive, Var1 == 6),
    aes(label = Var2),
    size = 2,
    box.padding = unit(0.3, "lines"),
    point.padding = unit(0.3, "lines"),
    max.overlaps = 40
  )+
  theme_classic()+
  labs(x = "index", y = "Percentage of isolates")
ggsave("image/phenotype_scores.pdf", height = 4, width = 5)
