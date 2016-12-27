# Preprocess dataset by combining field and satellite observations
# Thomas Nauss

if(Sys.info()["sysname"] == "Windows"){
  source("F:/analysis/orthoptera/orthoptera_prediction/src/00_set_environment.R")
} else {
  source("/media/tnauss/myWork/analysis/orthoptera/orthoptera_prediction/src/00_set_environment.R")
}

initOTB("C:/OSGeo4W64/bin/")
compute <- FALSE

# Prepare GLS2000 dataset ------------------------------------------------------
if(compute){
  gls <- stack(paste0(path_landsat, "gls2000.tif"))
  
  # Create one raster tile for each observation plot
  obsv_shp_wgs <- readRDS(file = paste0(path_results, "obsv_shp_wgs.rds"))
  obsv_shp_arc <- readRDS(file = paste0(path_results, "obsv_shp_arc.rds"))
  
  for(prj in c("wgs", "arc")){
    if(prj == "wgs"){
      obsv_shp <- obsv_shp_wgs
    } else {
      obsv_shp <- obsv_shp_arc
    }
    obsv_shp <- spTransform(obsv_shp, crs(gls))
    
    gls_snip <- snipRaster(raster=gls, spatial=obsv_shp, selector = "plot",
                           buffer=500, byid = TRUE)
    saveRDS(gls_snip, file = paste0(path_results, "gls_snip_", prj, ".rds"))
    
    # Compute pca for all plots
    gls_pca <- lapply(gls_snip, function(t){
      pca(stack(t), return_raster = TRUE)
    })
    saveRDS(gls_pca, file = paste0(path_results, "gls_pca_", prj, ".rds"))
    
    # Compute spectral indices for all plots
    mspec_indices <- lapply(seq(length(gls_snip)), function(i){
      if(i %% 10 == 0) print(paste0(i))
      mSpecIndices(blue = gls_snip[[i]][[1]], green = gls_snip[[i]][[2]], 
                   red = gls_snip[[i]][[3]], nir = gls_snip[[i]][[4]])
    })  
    saveRDS(mspec_indices, file = paste0(path_results, "mspec_indices_", prj, ".rds"))
    
    # Compute Haralick textures for all plots
    otb_txt <- lapply(seq(length(gls_pca)), function(i){
      if(i %% 10 == 0) print(paste0(i))
      otbTexturesHaralick(x=gls_pca[[i]], path_output = path_temp, 
                          return_raster = TRUE, texture="all",
                          channel = 1)
    })
    saveRDS(otb_txt, file = paste0(path_results, "otb_txt_", prj, ".rds"))
    
    # Compute glcm textures for all plots
    glcm_txt <- lapply(seq(length(gls_snip)), function(i){
      if(i %% 10 == 0) print(paste0(i))
      glcm(gls_pca[[i]][[1]], n_grey = 32, window = c(3,3),
           shift=list(c(0,1), c(1,1), c(1,0), c(1,-1)))
    })  
    saveRDS(glcm_txt, file = paste0(path_results, "glcm_txt_", prj, ".rds"))
  }
  
  # Combine results in one stack per tile
  gls_2000 <- lapply(seq(length(gls_snip)), function(i){
    stack(gls_snip[[i]], gls_pca[[i]], mspec_indices[[i]],
          otb_txt[[i]], glcm_txt[[i]])
  })
  saveRDS(gls_2000, file = paste0(path_results, "gls_2000_", prj, ".rds"))
  
} else {
  #   gls_snip_wgs <- readRDS(file = paste0(path_results, "gls_snip_", prj, ".rds"))
  #   gls_pca_wgs <- readRDS(file = paste0(path_results, "gls_pca_", prj, ".rds"))
  #   mspec_indices_wgs <- readRDS(file = paste0(path_results, "mspec_indices_", prj, ".rds"))
  #   otb_txt_wgs <- readRDS(file = paste0(path_results, "otb_txt_", prj, ".rds"))
  #   glcm_txt_wgs <- readRDS(file = paste0(path_results, "glcm_txt_", prj, ".rds"))
  gls_2000_wgs <- readRDS(file = paste0(path_results, "gls_2000_wgs.rds"))
  
  #   gls_snip_arc <- readRDS(file = paste0(path_results, "gls_snip_", prj, ".rds"))
  #   gls_pca_arc <- readRDS(file = paste0(path_results, "gls_pca_", prj, ".rds"))
  #   mspec_indices_arc <- readRDS(file = paste0(path_results, "mspec_indices_", prj, ".rds"))
  #   otb_txt_arc <- readRDS(file = paste0(path_results, "otb_txt_", prj, ".rds"))
  #   glcm_txt_arc <- readRDS(file = paste0(path_results, "glcm_txt_", prj, ".rds"))
  gls_2000_arc <- readRDS(file = paste0(path_results, "gls_2000_arc.rds"))
}


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


# ggplot(data = ndvi_plots_final@data, aes(x = NDVI_WGS, y = NDVI_ARC)) + 
#   geom_point() + 
#   geom_smooth()

saveRDS(ndvi_plots_final, 
        file = paste0(path_results, "ndvi_plots_final.RDS"))
saveRDS(as.data.frame(ndvi_plots_final), 
        file = paste0(path_results, "ndvi_plots_final_df.RDS"))

