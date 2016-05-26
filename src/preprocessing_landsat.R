# Preprocess Landsat data by computing mulit-scene mean values for each pixel
library(rgdal)
library(satellite)
library(sp)
library(glcm)

if(Sys.info()["sysname"] == "Windows"){
  filepath_base <- "E:/analysis/orthoptera/"
} else {
  filepath_base <- "/media/tnauss/myWork/analysis/orthoptera/"
}

rasterOptions(tmpdir=paste0(filepath_base, "data/temp"))

filepath_sat <- paste0(filepath_base, "data/landsat/GLS2000")
filepath_obsv <- paste0(filepath_base, "data/original/")

landsat_files <- list.files(filepath_sat, pattern = glob2rx("*.tif"), 
                            full.names = TRUE)

obsv <- readOGR(paste0(
  filepath_obsv, "lvl0300_biodiversity_data_utm37s_buffer_dissolve_single.shp"),
  "lvl0300_biodiversity_data_utm37s_buffer_dissolve_single")

landsat <- lapply(unique(substr(basename(landsat_files), 1, 4)), function(x){
  sat <- satellite(landsat_files[grepl(x, landsat_files)])
  sat <- crop(sat, extent(obsv))
  # sat <- alignGeometry(sat, getSatDataLayer(sat, bcde = "Bxx3"), method = "ngb")
  return(sat)
})

save(landsat, file = "E:/analysis/orthoptera/data/processed/landsat.RData")
# load("E:/analysis/orthoptera/data/processed/landsat.RData")


landsat_01_stack <- stack(subset(landsat[[1]], sid = c("Bxx4",
                                                       "Bxx5", "Bxx6", "Bxx7", "Bxx8", 
                                                       "Bxx9")))
landsat_02_stack <- stack(subset(landsat[[2]], sid = c("Bxx4",
                                                       "Bxx5", "Bxx6", "Bxx7", "Bxx8", 
                                                       "Bxx9")))
NAvalue(landsat_01_stack) <- 0
NAvalue(landsat_02_stack) <- 0

writeRaster(landsat_01_stack, paste0(filepath_base, "data/landsat/gls2000_01.tif"), overwrite = TRUE)
writeRaster(landsat_02_stack, paste0(filepath_base, "data/landsat/gls2000_02.tif"), overwrite = TRUE)



gls_scale <- raster(paste0(filepath_base, "data/landsat/gls2000_pca01_scale_hist.tif"))
gls_scale_glcm <- glcm(gls_scale)
names(gls_scale_glcm)
writeRaster(gls_scale_glcm, bylayer = "TRUE",
            filename = paste0(filepath_base, "data/landsat/gls2000_pca01_scale_hist_", names(gls_scale_glcm)), 
            format = "GTiff")
