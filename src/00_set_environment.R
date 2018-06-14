# Set path ---------------------------------------------------------------------
if(Sys.info()["sysname"] == "Windows"){
  path_base <- "E:/tnauss/orthoptera/"
  path_base <- "E:/analysis/orthoptera/"
} else {
  path_base <- "/media/tnauss/myWork/analysis/orthoptera/"
}

path_data <- paste0(path_base, "data/")
path_landsat <- paste0(path_data, "landsat/")
path_obsv <- paste0(path_data, "orthoptera/")
path_results <- paste0(path_data, "rdata/")
path_temp <- paste0(path_data, "temp/")


path_modis_local_arc <- paste0(path_data, "modis/modis_arc/")
path_modis_out_dir <- paste0(path_data, "modis/processed/")
path_modis_crop <-paste0(path_data, "modis/croped/")



# Set libraries ----------------------------------------------------------------
# devtools::install_github("MatMatt/MODIS", ref = "develop")
# setRepositories() 
# install.packages(c(' mapdata', 'ptw '),dependencies=TRUE)

library(doParallel)
library(glcm)
library(gpm)
library(ggplot2)
library(mapview)
library(MODIS)
library(pROC)
library(raster)
library(readxl)
library(rgdal)
library(rgeos)
library(satellite)
library(satelliteTools)
library(sp)


# Other settings ---------------------------------------------------------------
rasterOptions(tmpdir = path_temp)

saga_cmd <- "C:/OSGeo4W64/apps/saga/saga_cmd.exe "
# initOTB("C:/OSGeo4W64/bin/")
initOTB("C:/PROGRA~1/OTB-5.8.0-win64/bin/")