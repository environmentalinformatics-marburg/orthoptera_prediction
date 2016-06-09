# Preprocess dataset by combining field and satellite observations
library(mapview)
library(raster)
library(rgdal)
library(sp)


if(Sys.info()["sysname"] == "Windows"){
  filepath_base <- "E:/analysis/orthoptera/data/"
} else {
  filepath_base <- "/media/tnauss/myWork/analysis/orthoptera/data/"
}

rasterOptions(tmpdir=paste0(filepath_base, "temp"))

filepath_landsat <- paste0(filepath_base, "landsat/")
filepath_modis <- paste0(filepath_base, "modis/processed/")
filepath_obsv <- paste0(filepath_base, "orthoptera/")
filepath_results <- paste0(filepath_base, "rdata/")

# Read observations
obsv <- read.table(paste0(filepath_obsv, "lvl0300_biodiversity_data.csv"),
                   sep = ";", dec = ",", header = TRUE)
obsv <- obsv[, -grep("greyval" , colnames(obsv))]
colnames(obsv)

shp <- obsv
coordinates(shp) <- ~coordW+coordN
projection(shp) <- CRS("+proj=utm +zone=37 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs ")
# mapview(shp)


# Read Landsat spatial datasets and attach to observations
landsat_files <- list.files(filepath_landsat, 
                            pattern = glob2rx("gls2000_pca01_scale_hist*.tif"), 
                            full.names = TRUE)
landsat <- stack(landsat_files)
shp_tf <- spTransform(shp, projection(landsat))
landsat_vals <- extract(landsat, shp_tf)
colnames(landsat_vals) <- paste0("gls2000_", substr(colnames(landsat_vals), 21, nchar(colnames(landsat_vals))))
obsv <- cbind(obsv, landsat_vals)
# save(landsat_vals, file = paste0(filepath_results, "preprocessing_data_landsat_vals.RData"))

# Read MODIS datasets and attach to observations
product <- c("MOD09GA", "MYD09GA")
mod_files <- list.files(paste0(filepath_modis, product[1], "_TZS"), 
                        pattern = glob2rx("*sur_refl_b0*tif"), 
                        full.names = TRUE)
myd_files <- list.files(paste0(filepath_modis, product[2], "_TZS"), 
                        pattern = glob2rx("*sur_refl_b0*tif"), 
                        full.names = TRUE)
modis_files <- c(mod_files, myd_files)

modis <- stack(modis_files)
shp_tf <- spTransform(shp, projection(modis))
modis_vals <- extract(modis, shp_tf)
modis_vals <- as.data.frame(modis_vals)
# save(modis_vals, file = paste0(filepath_results, "preprocessing_data_modis_vals.RData"))

bandnames <- c("sur_refl_b01_1", "sur_refl_b02_1", "sur_refl_b03_1", 
               "sur_refl_b04_1", "sur_refl_b05_1", "sur_refl_b06_1",
               "sur_refl_b07_1")
modis_obsv <- lapply(seq(nrow(obsv)), function(x){
  date <- gsub("-", "", substring(obsv$date_nocloud[x], 1, 8), "-")
  cnames <- c(paste0(product[1], ".A", date, ".", bandnames), 
              paste0(product[2], ".A", date, ".", bandnames))
  act_vals <- modis_vals[x, which(colnames(modis_vals) %in% cnames)]
  colnames(act_vals) <- sapply(strsplit(colnames(act_vals), "\\."), function(y){
    paste0(y[1], "_", y[3])
  })
  if(!any(grepl(product[1], colnames(act_vals)))){
    missing <- data.frame(t(rep(NA, 7)))
    colnames(missing) <- sapply(strsplit(cnames[1:7], "\\."), function(y){
      paste0(y[1], "_", y[3])
    })
    act_vals <- cbind(missing, act_vals)
    
  } else if(!any(grepl(product[2], colnames(act_vals)))){
    missing <- data.frame(t(rep(NA, 7)))
    colnames(missing) <- sapply(strsplit(cnames[8:14], "\\."), function(y){
      paste0(y[1], "_", y[3])
    })
    act_vals <- cbind(act_vals, missing)
  }
  return(cbind(data.frame(date_modis = date),act_vals))
})
modis_obsv <- do.call("rbind", modis_obsv)
obsv <- cbind(obsv, modis_obsv)

obsv <- obsv[, c(1:grep("diff_days_nocloud", colnames(obsv)),
                 grep("date_modis", colnames(obsv)),
                 (grep("diff_days_nocloud", colnames(obsv))+1):(grep("date_modis", colnames(obsv))-1),
                 (grep("date_modis", colnames(obsv))+1):length(colnames(obsv)))]
all.equal(gsub("-", "", substring(obsv$date_nocloud, 1, 8), "-"),
          as.character(obsv$date_modis))
# save(modis_obsv, file = paste0(filepath_results, "preprocessing_data_modis_obsv.RData"))

# Write dataset
write.table(obsv, 
            file = paste0(filepath_obsv, "orthoptera_sat.csv"),
            sep = ",", dec = ".", row.names = FALSE)
save(obsv, file = paste0(filepath_results, "preprocessing_data_obsv.RData"))

# Read Landsat spectral datasets and attach to observations (came in to a later analysis stage)
landsat_files <- list.files(filepath_landsat, 
                            pattern = glob2rx("gls2000.tif"), 
                            full.names = TRUE)
landsat <- stack(landsat_files)
shp_tf <- spTransform(shp, projection(landsat))
landsat_vals <- extract(landsat, shp_tf)
colnames(landsat_vals) <- paste0("gls2000_", c(seq(5), 7))
obsv <- cbind(obsv, landsat_vals)

# Write dataset
write.table(obsv, 
            file = paste0(filepath_obsv, "orthoptera_sat.csv"),
            sep = ",", dec = ".", row.names = FALSE)
save(obsv, file = paste0(filepath_results, "preprocessing_data_obsv.RData"))
