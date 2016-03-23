setwd("D:/active/orthoptera/data")

# Libraries --------------------------------------------------------------------
library(raster)
library(rgdal)
library(sp)



# Convert csv 2 shp ------------------------------------------------------------
orthoptera <- read.table("original/lvl0300_biodiversity_data.csv", 
                         header = TRUE, sep = ";", dec = ",")

# Replace number of observations and NAs to 1/0
orthoptera[, 14:178][!is.na(orthoptera[, 14:178])] <- "yes"
orthoptera[, 14:178][is.na(orthoptera[, 14:178])] <- "no"
for(i in seq(14, 178)){
  orthoptera[, i] <- as.factor(orthoptera[, i])
}

# Compile dataset containing complete cases only
orthoptera <- 
  orthoptera[, -(which(colnames(orthoptera) == "greyval_band_11") : 
                   which(colnames(orthoptera) == "greyval_band_16"))]
any(is.na(orthoptera[, -7]))

coordinates(orthoptera) <- ~lon+lat
projection(orthoptera) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs ")

writeOGR(orthoptera, "original/lvl0300_biodiversity_data.shp", 
         "lvl0300_biodiversity", driver="ESRI Shapefile")
