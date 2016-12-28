# Preprocess dataset by combining field and satellite observations

if(Sys.info()["sysname"] == "Windows"){
  source("F:/analysis/orthoptera/orthoptera_prediction/src/00_set_environment.R")
} else {
  source("/media/tnauss/myWork/analysis/orthoptera/orthoptera_prediction/src/00_set_environment.R")
}

# Prepare GLS2000 dataset ------------------------------------------------------
gls <- stack(paste0(filepath_landsat, "gls2000.tif"))
# mapview(gls) + obsv_shp_arc

ndvi <- (gls[[4]] - gls[[3]]) / (gls[[4]] + gls[[3]])


# Prepare orthoptera observations ----------------------------------------------
obsv <- read_excel(paste0(filepath_obsv, "Grasshopper-Data.xlsx"))  
obsv <- as.data.frame(obsv)
obsv$date_observation <- format(as.Date(obsv$date, "%d/%m/%Y"), "%Y-%j")

obsv_shp_wgs <- obsv
coordinates(obsv_shp_wgs) <- ~coordW+coordN
projection(obsv_shp_wgs) <- CRS("+init=epsg:32737")

obsv_shp_arc <- obsv
coordinates(obsv_shp_arc) <- ~coordW+coordN
projection(obsv_shp_arc) <- CRS("+init=epsg:21037")

obsv_shp_wgs <- spTransform(obsv_shp_wgs, crs(ndvi))
obsv_shp_arc <- spTransform(obsv_shp_arc, crs(ndvi))
# mapview(obsv_shp_wgs)+obsv_shp_arc


# Extract GLS2000 data ---------------------------------------------------------
ndvi_obs <- lapply(c(obsv_shp_wgs, obsv_shp_arc), function(obsv_shp){
  ndvi_plots <- extract(ndvi, obsv_shp, sp = TRUE)
  colnames(ndvi_plots@data)[ncol(ndvi_plots@data)] <- "NDVI"

  ndvi_plots_buffer <- extract(ndvi, obsv_shp, buffer = 60.0)
  ndvi_plots_buffer_stat <- lapply(seq(length(ndvi_plots_buffer)), function(i){
    data.frame(ID = obsv_shp@data[i,"ID"],
               NDVI_mean = mean(ndvi_plots_buffer[[i]]),
               NDVI_median = median(ndvi_plots_buffer[[i]]),
               NDVI_sd = sd(ndvi_plots_buffer[[i]]),
               NDVI_min = min(ndvi_plots_buffer[[i]]),
               NDVI_max = max(ndvi_plots_buffer[[i]]))
  })
  ndvi_plots_buffer_stat <- do.call("rbind", ndvi_plots_buffer_stat)
  merge(ndvi_plots, ndvi_plots_buffer_stat)
})

colnames(ndvi_obs[[1]]@data)[28:33] <- paste0(colnames(ndvi_obs[[1]]@data)[28:33], "_WGS")
colnames(ndvi_obs[[2]]@data)[28:33] <- paste0(colnames(ndvi_obs[[2]]@data)[28:33], "_ARC")
ndvi_plots_final <- merge(ndvi_obs[[1]], ndvi_obs[[2]]@data)
head(ndvi_plots_final@data)


saveRDS(ndvi_plots_final, 
        file = paste0(filepath_results, "ndvi_plots_final.RDS"))
saveRDS(as.data.frame(ndvi_plots_final), 
        file = paste0(filepath_results, "ndvi_plots_final_df.RDS"))

