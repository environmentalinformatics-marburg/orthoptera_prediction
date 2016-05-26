# Preprocess dataset by combining field and satellite observations
library(mapview)
library(raster)
library(rgdal)
library(sp)


if(Sys.info()["sysname"] == "Windows"){
  filepath_base <- "E:/analysis/orthoptera/"
} else {
  filepath_base <- "/media/tnauss/myWork/analysis/orthoptera/"
}

rasterOptions(tmpdir=paste0(filepath_base, "data/temp"))

filepath_landsat <- paste0(filepath_base, "data/landsat/")
filepath_modis <- paste0(filepath_base, "data/modis/processed/")
filepath_obsv <- paste0(filepath_base, "data/orthoptera/")


# Read observations
obsv <- read.table(paste0(filepath_obsv, "lvl0300_biodiversity_data.csv"),
                   sep = ";", dec = ",", header = TRUE)
obsv <- obsv[, -grep("greyval" , colnames(obsv))]
colnames(obsv)

shp <- obsv
coordinates(shp) <- ~coordW+coordN
projection(shp) <- CRS("+proj=utm +zone=37 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs ")
# mapview(shp)


# Read Landsat datasets and attach to observations
landsat_files <- list.files(filepath_landsat, 
                            pattern = glob2rx("gls2000_pca01_scale_hist*.tif"), 
                            full.names = TRUE)
lnds <- stack(landsat_files)
shp <- spTransform(shp, projection(lnds))
lnds_vals <- extract(lnds, shp)
colnames(lnds_vals) <- paste0("gls2000_", substr(colnames(lnds_vals), 21, nchar(colnames(lnds_vals))))
obsv <- cbind(obsv, lnds_vals)

# Read MODIS datasets and attach to observations

# Write dataset
write.table(as.data.frame(obsv), 
            file = paste0(filepath_obsv, "orthoptera_sat.csv"),
            sep = ",", dec = ".", row.names = FALSE)
# writeOGR(obsv, paste0(filepath_obsv, "orthoptera_sat.shp"), "orthoptera_sat", 
#          driver="ESRI Shapefile")
