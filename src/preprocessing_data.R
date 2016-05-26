# Preprocess dataset by combining field and satellite observations
library(mapview)
library(raster)
library(rgdal)
library(sp)


if(Sys.info()["sysname"] == "Windows"){
  filepath_base <- "G:/analysis/orthoptera/data/"
} else {
  filepath_base <- "/media/tnauss/myWork/analysis/orthoptera/data/"
}

rasterOptions(tmpdir=paste0(filepath_base, "temp"))

filepath_landsat <- paste0(filepath_base, "landsat/")
filepath_modis <- paste0(filepath_base, "modis/processed/")
filepath_obsv <- paste0(filepath_base, "orthoptera/")


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
modis_files <- list.files(paste0(filepath_modis, "MOD09GA_TZS"), 
                          pattern = glob2rx("*sur_refl_b0*tif"), 
                          full.names = TRUE)
modis <- stack(modis_files)
shp <- spTransform(shp, projection(modis))
modis_vals <- extract(modis, shp)

# Write dataset
write.table(as.data.frame(obsv), 
            file = paste0(filepath_obsv, "orthoptera_sat.csv"),
            sep = ",", dec = ".", row.names = FALSE)
# writeOGR(obsv, paste0(filepath_obsv, "orthoptera_sat.shp"), "orthoptera_sat", 
#          driver="ESRI Shapefile")
