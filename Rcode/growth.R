# set working directory
setwd("~/Desktop/Jiali/TAMU/Susie/Fungi/")
library(readxl)
library(ggplot2)
library(ggpmisc)
fileName <- "PDA result of USDA fungus08182022.xlsx"
excel_sheets(fileName)
Screening<- read_excel(fileName, sheet = 1, skip = 3)
# convert contaminant to NA
Screening[Screening == 'contaminated'] <- NA
Screening[Screening == 'contaminate'] <- NA
Screening[Screening == 'N/A'] <- NA
# 6-14 day2
# 15-23 day3
# 24-32 day4
# 33-41 day5
# 42-50 day6
# 51-59 day7
# 60-68 day8
# 69-77 day9
# 78-86 day10
# 87-95 day11
# 96-104 day12
# 105-113 day13
# 114-122 day14

# reformat the whole table
library(reshape2)
tbl_melt <- melt(Screening, id.vars = c("No.","batch","Collection Number","Blast ID(ITS)","Scientific Name"))
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

#ID = Screening$No.[10]

# plot and calculate
ScreenDye <- function(ID){
data <- tbl_melt[which(tbl_melt$No.==ID),]
data <- na.omit(data)
day1 <- data[1:3,1:6]
day1$value = c(0,0,0)
day1$Day=c(1,1,1)
day1$media=c("PDA","PDA+Dye","Dye")
day1$replicate=c(1,1,1)
data <- rbind(day1,data)
# select day reach dye decolor to 5
dye_tbl <- data[data$media =="Dye",]
cutoff <- dye_tbl$Day[dye_tbl$value==5][1]
if (is.na(cutoff) == FALSE) {
data_filter <- data[which(data$Day <= cutoff),]

# plot
plot <- ggplot(data=data_filter, aes(x=Day, y=value, color=media, group=media))+
  geom_smooth(method = "loess") +
  geom_point() +
  theme_classic()+
  labs(x="Culture days", y ="Measurement", title = ID)
print(plot)
return(output_tbl)
}}

# Loop over all strains
pdf('Growth_plot 0818.pdf')
for (strain in Screening$No.) {
ScreenDye(strain) 
}
dev.off()
# filter again and output into a table
output_tbl <- data.frame(matrix(ncol=3, nrow = 0))
colnames(output_tbl) <- c("No.","Dye.4.Day", "Grow.Max.Day")
for (strain in Screening$No.) {
  data <- tbl_melt[which(tbl_melt$No.==strain),]
  data <- na.omit(data)
  # select day reach dye decolor to 5
  dye_tbl <- data[data$media =="Dye",]
  cutoff <- dye_tbl$Day[dye_tbl$value==4][1]
  growth_tbl <- data[data$media =="PDA+Dye",]
  maxgrowth <- growth_tbl$Day[growth_tbl$value==5][1]
  if (is.na(cutoff) == FALSE) {
    output_tbl[nrow(output_tbl)+1, ] <- c(strain, cutoff, maxgrowth)
}
}
# add strain annotation
output_tbl <- merge(output_tbl, Screening[,1:5], by = "No.", all.x = T)
write.csv(output_tbl, "Screening output 0818.csv", row.names = F)
