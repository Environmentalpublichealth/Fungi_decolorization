setwd("~/Desktop/Jiali/TAMU/Susie/Fungi/")
# install
install.packages("raster")
install.packages("spData")
install.packages("spDataLarge")
install.packages("tmap")

list_fungi <- read.csv("list of USDA fungi 7-21-2022.xlsx - Sheet1.csv", header = T)
# correct some typos
list_fungi[] <- lapply(list_fungi, function(x) gsub("Hypsizygous","Hypsizygus", x))
list_fungi[] <- lapply(list_fungi, function(x) gsub("Cerrrena","Cerrena", x))
list_fungi[] <- lapply(list_fungi, function(x) gsub("B. subcoronatum","Botryobasidium subcoronatum", x))
list_fungi[] <- lapply(list_fungi, function(x) gsub("C. subverracisporus 1","Crepidotus subverrucisporus", x))
list_fungi[] <- lapply(list_fungi, function(x) gsub("P. chrysosporium","Phanerochaete chrysosporium", x))

list_fungi$Blast.ID.ITS. <- sapply(list_fungi$Blast.ID.ITS., function(y) paste(unlist(strsplit(y, " "))[1:2], collapse = " "))

species <- data.frame(table(list_fungi$Blast.ID.ITS.))
# Match state latitude and longitude
list_fungi$State <- sapply(list_fungi$Location, function(y) unlist(strsplit(y, ";|,"))[1])
states_df <- data.frame("State" = list_fungi$State)

# make maps
library(sf)
library(raster)
library(dplyr)
library(spData)
library(tmap)
data(us_states)
state_isolates <- data.frame(table(list_fungi$State))
us_states <- merge(us_states, state_isolates, by.x = "NAME", by.y = "Var1", all.x = T)
names(us_states)[7] <- "Isolates"
us_states_map = tm_shape(us_states, projection = 2163) + tm_polygons() + 
  tm_layout(frame = FALSE, legend.position = c(0.85,0.15)) + tm_dots(col = "black", size = "Isolates")

data("hawaii")
data("alaska")
hawaii <- merge(hawaii, state_isolates, by.x = "NAME", by.y = "Var1", all.x = T)
hawaii_map = tm_shape(hawaii) + tm_polygons() + tm_dots(col = "black", size = 0.014, legend.size.show=F)+
  tm_layout(title = "Hawaii", frame = FALSE, bg.color = NA, 
            title.position = c("LEFT", "BOTTOM"))
alaska <- merge(alaska, state_isolates, by.x = "NAME", by.y = "Var1", all.x = T)
alaska_map = tm_shape(alaska) + tm_polygons() + tm_dots(col = "black", size = 0.23, legend.size.show=F) +
  tm_layout(title = "Alaska", title.size = 0.9, frame = FALSE, bg.color = NA)

pdf("collection map.pdf", height = 4, width = 6.3)
us_states_map
print(hawaii_map, vp = grid::viewport(0.35, 0.1, width = 0.2, height = 0.1))
print(alaska_map, vp = grid::viewport(0.15, 0.15, width = 0.3, height = 0.3))
dev.off()

table(list_fungi$Country)

table(state_isolates)
# high decolor
list_fungi_good <- list_fungi[list_fungi$No. %in% goodSpecies,]
