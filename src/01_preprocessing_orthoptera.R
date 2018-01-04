# Preprocess dataset by combining field and satellite observations
# Thomas Nauss

if(Sys.info()["sysname"] == "Windows"){
  source("F:/analysis/orthoptera/orthoptera_prediction/src/00_set_environment.R")
} else {
  source("/media/tnauss/myWork/analysis/orthoptera/orthoptera_prediction/src/00_set_environment.R")
}


# Prepare orthoptera observations ----------------------------------------------
# Read original data provided by CH as well as preprocessed data by YT and merge
# datasets.
obsv_ch <- read_excel(paste0(path_obsv, "abundance_matrix_hemp.xls"))  
colnames(obsv_ch)[which(colnames(obsv_ch) == "asl")] <- "elevation"

obsv_yt <- read_excel(paste0(path_obsv, "Grasshopper-Data.xlsx"))  
obsv_yt$plot <- as.numeric(obsv_yt$plot)
obsv_yt$date <- as.POSIXct(strptime(obsv_yt$date, "%d/%m/%Y"), tz = "UTC")

str(obsv_yt[, c("plot", "date", "coordN", "coordW")])
str(obsv_ch[, c("plot", "date", "coordN", "coordW")])

obsv <- merge(obsv_yt, obsv_ch, by = c("plot", "date", "coordN", "coordW",
                                       "elevation", "rich.plants", "rainfall"))

# obsv_shp_wgs <- obsv
# coordinates(obsv_shp_wgs) <- ~coordW+coordN
# projection(obsv_shp_wgs) <- CRS("+init=epsg:32737")
# saveRDS(obsv_shp_wgs, file = paste0(path_results, "obsv_shp_wgs.rds"))

obsv_shp_arc <- obsv
coordinates(obsv_shp_arc) <- ~coordW+coordN
projection(obsv_shp_arc) <- CRS("+init=epsg:21037")
saveRDS(obsv_shp_arc, file = paste0(path_results, "obsv_shp_arc.rds"))
