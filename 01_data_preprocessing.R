library(dplyr)
library(readr)

# Must set column name as "sp", "Longitude", "Latitude"

setwd("D:/Envdata/")
SpeciesR <- read.csv("dataset.csv", stringsAsFactors = FALSE) # Set species data(whole) file 
colnames(SpeciesR)
str(SpeciesR)

sp <- SpeciesR[!is.na(SpeciesR$Longitude) & !is.na(SpeciesR$Latitude), ]

check <- sp %>% group_by(sp) %>% summarize(n = n()) %>% dplyr::filter(n < 20) # Remove species which is sample number < 20

check

check <- as.data.frame(check)

sp <- sp %>%  anti_join(check)
Splist <- unique(sp$sp)

setwd("D:/Envdata/Species")


for(i in Splist){
  readr :: write_csv(filter(SpeciesR, sp == i), file = paste0(i, ".csv"))
}

