# Set path ---------------------------------------------------------------------
if(Sys.info()["sysname"] == "Windows"){
  path_base <- "F:/analysis/orthoptera/"
} else {
  path_base <- "/media/tnauss/myWork/analysis/orthoptera/"
}

path_data <- paste0(path_base, "data/")
path_landsat <- paste0(path_data, "landsat/")
path_obsv <- paste0(path_data, "orthoptera/")
path_results <- paste0(path_data, "rdata/")
path_temp <- paste0(path_data, "temp/")


# Set libraries ----------------------------------------------------------------
library(rgeos)
library(ggplot2)
library(mapview)
library(raster)
library(readxl)
library(rgdal)
library(satellite)
library(satelliteTools)
library(sp)
