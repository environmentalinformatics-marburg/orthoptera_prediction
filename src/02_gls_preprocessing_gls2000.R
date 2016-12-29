# Preprocess GLS200 dataset
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
                          return_raster = TRUE, 
                          parameters.xyrad=list(c(1,1)),
                          parameters.xyoff=list(c(1,1)),
                          texture="all",
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
    # Combine results in one stack per tile
    gls_2000 <- lapply(seq(length(gls_snip)), function(i){
      stack(gls_snip[[i]], gls_pca[[i]], mspec_indices[[i]],
            otb_txt[[i]], glcm_txt[[i]])
    })
    names(gls_2000) <- names(gls_snip)
    saveRDS(gls_2000, file = paste0(path_results, "gls_2000_", prj, ".rds"))
  }
  
} else {
  #     gls_snip_wgs <- readRDS(file = paste0(path_results, "gls_snip_", prj, ".rds"))
  #     gls_pca_wgs <- readRDS(file = paste0(path_results, "gls_pca_", prj, ".rds"))
  #     mspec_indices_wgs <- readRDS(file = paste0(path_results, "mspec_indices_", prj, ".rds"))
  #     otb_txt_wgs <- readRDS(file = paste0(path_results, "otb_txt_", prj, ".rds"))
  #     glcm_txt_wgs <- readRDS(file = paste0(path_results, "glcm_txt_", prj, ".rds"))
  gls_2000_wgs <- readRDS(file = paste0(path_results, "gls_2000_wgs.rds"))
  
  #   gls_snip_arc <- readRDS(file = paste0(path_results, "gls_snip_", prj, ".rds"))
  #   gls_pca_arc <- readRDS(file = paste0(path_results, "gls_pca_", prj, ".rds"))
  #   mspec_indices_arc <- readRDS(file = paste0(path_results, "mspec_indices_", prj, ".rds"))
  #   otb_txt_arc <- readRDS(file = paste0(path_results, "otb_txt_", prj, ".rds"))
  #   glcm_txt_arc <- readRDS(file = paste0(path_results, "glcm_txt_", prj, ".rds"))
  gls_2000_arc <- readRDS(file = paste0(path_results, "gls_2000_arc.rds"))
}



