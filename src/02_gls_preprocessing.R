# Preprocess GLS200 dataset
# Thomas Nauss

if(Sys.info()["sysname"] == "Windows"){
  source("F:/analysis/orthoptera/orthoptera_prediction/src/00_set_environment.R")
} else {
  source("/media/tnauss/myWork/analysis/orthoptera/orthoptera_prediction/src/00_gls_set_environment.R")
}

compute <- FALSE

# Prepare GLS2000 dataset ------------------------------------------------------
if(compute){
  gls <- stack(paste0(path_landsat, "gls2000.tif"))
  
  obsv_shp <- readRDS(file = paste0(path_results, "obsv_shp_arc.rds"))
  obsv_shp <- spTransform(obsv_shp, crs(gls))
  
  # Create one raster tile for each observation plot
  gls_snip <- snipRaster(raster=gls, spatial=obsv_shp, selector = "plot",
                         buffer=500, byid = TRUE)
  saveRDS(gls_snip, file = paste0(path_results, "gls_snip_arc.rds"))
  
  # Compute pca for all plots
  gls_pca <- pca(gls_snip, center = TRUE, scale = TRUE)
  saveRDS(gls_pca, file = paste0(path_results, "gls_pca_arc.rds"))
  
  # Compute spectral indices for all plots
  gls_mspec_indices <- lapply(seq(length(gls_snip)), function(i){
    if(i %% 10 == 0) print(paste0(i))
    mSpecIndices(blue = gls_snip[[i]][[1]], green = gls_snip[[i]][[2]], 
                 red = gls_snip[[i]][[3]], nir = gls_snip[[i]][[4]])
  })  
  names(gls_mspec_indices) <- names(gls_snip)
  saveRDS(gls_mspec_indices, file = paste0(path_results, "gls_mspec_indices_arc.rds"))
  
  # Compute Haralick textures for all plots based on PCA
  windows <- seq(1, 25, 2)

  minv <- min(unlist(lapply(gls_pca, function(s){minValue(s$PC1)})))
  maxv <- max(unlist(lapply(gls_pca, function(s){maxValue(s$PC1)})))
  
  gls_pca_otb_txt <- lapply(seq(length(gls_pca)), function(i){
    if(i %% 10 == 0) print(paste0(i))
    otb_txt_w <- lapply(windows, function(w){
      oth <- otbTexturesHaralick(x=gls_pca[[i]]$PC1, path_output = path_temp, 
                                 return_raster = TRUE, 
                                 parameters.xyrad=list(c(w,w)),
                                 parameters.xyoff=list(c(1,1)),
                                 parameters.minmax=c(minv, maxv),
                                 parameters.nbbin = 8,
                                 texture="all",
                                 channel = 1)
      names(oth) <- paste0("pca_", names(oth))
      return(oth)
    })
    return(otb_txt_w)
  })
  saveRDS(gls_pca_otb_txt, file = paste0(path_results, "gls_pca_otb_txt_arc_bu.rds"))
  names(gls_pca_otb_txt) <- names(gls_snip)
  saveRDS(gls_pca_otb_txt, file = paste0(path_results, "gls_pca_otb_txt_arc_new.rds"))
  
  # Compute Haralick textures for all plots based on NDVI
  minv <- -1
  maxv <- 1
  gls_ndvi_otb_txt <- lapply(seq(length(gls_mspec_indices)), function(i){
    if(i %% 10 == 0) print(paste0(i))
    otb_txt_w <- lapply(windows, function(w){
      oth <- otbTexturesHaralick(x=gls_mspec_indices[[i]]$NDVI, path_output = path_temp, 
                                 return_raster = TRUE, 
                                 parameters.xyrad=list(c(w,w)),
                                 parameters.xyoff=list(c(1,1)),
                                 parameters.minmax=c(minv, maxv),
                                 parameters.nbbin = 8,
                                 texture="all",
                                 channel = 1)
      names(oth) <- paste0("ndvi_", names(oth))
      return(oth)
    })
    return(otb_txt_w)
  })
  names(gls_ndvi_otb_txt) <- names(gls_snip)
  saveRDS(gls_ndvi_otb_txt, file = paste0(path_results, "gls_ndvi_otb_txt_arc.rds"))
  
  # Compute glcm textures for all plots based on PCA
  gls_pca_glcm_txt <- lapply(seq(length(gls_pca)), function(i){
    if(i %% 10 == 0) print(paste0(i))
    gt_w <- lapply(windows, function(w){
      gt <- glcm(gls_pca[[i]]$PC1, n_grey = 32, window = c(w,w),
                 shift=list(c(0,1), c(1,1), c(1,0), c(1,-1)))
      names(gt) <- paste0("pca_", names(gt), ".r", w)
      return(gt)
    })
    return(gt_w)
  })  
  names(gls_pca_glcm_txt) <- names(gls_snip)
  saveRDS(gls_pca_glcm_txt, file = paste0(path_results, "gls_pca_glcm_txt_arc.rds"))
  
  # Compute glcm textures for all plots based on NDVI
  gls_ndvi_glcm_txt <- lapply(seq(length(gls_mspec_indices)), function(i){
    if(i %% 10 == 0) print(paste0(i))
    gt_w <- lapply(windows, function(w){
      gt <- glcm(gls_mspec_indices[[i]]$NDVI, n_grey = 32, window = c(3,3),
               shift=list(c(0,1), c(1,1), c(1,0), c(1,-1)))
    names(gt) <- paste0("ndvi_", names(gt), ".r", w)
    return(gt)
    })
    return(gt_w)
  })  
  names(gls_ndvi_glcm_txt) <- names(gls_snip)
  saveRDS(gls_ndvi_glcm_txt, file = paste0(path_results, "gls_ndvi_glcm_txt_arc.rds"))
  
  # Combine results in one stack per tile
  gls_2000 <- lapply(seq(length(gls_snip)), function(i){
    stack(gls_snip[[i]], gls_pca[[i]], gls_mspec_indices[[i]],
          stack(gls_pca_otb_txt[[i]]), stack(gls_ndvi_otb_txt[[i]]),
          stack(gls_pca_glcm_txt[[i]]), stack(gls_ndvi_glcm_txt[[i]]))
  })
  names(gls_2000) <- names(gls_snip)
  saveRDS(gls_2000, file = paste0(path_results, "gls_2000_arc.rds"))
  
  
} else {
  gls_snip <- readRDS(file = paste0(path_results, "gls_snip_arc.rds"))
  gls_pca <- readRDS(file = paste0(path_results, "gls_pca_arc.rds"))
  gls_mspec_indices <- readRDS(file = paste0(path_results, "gls_mspec_indices_arc.rds"))
  gls_pca_otb_txt <- readRDS(file = paste0(path_results, "gls_pca_otb_txt_arc.rds"))
  gls_ndvi_otb_txt <- readRDS(file = paste0(path_results, "gls_ndvi_otb_txt_arc.rds"))
  gls_pca_glcm_txt <- readRDS(file = paste0(path_results, "gls_pca_glcm_txt_arc.rds"))
  gls_ndvi_glcm_txt <- readRDS(file = paste0(path_results, "gls_ndvi_glcm_txt_arc.rds"))
  gls_2000 <- readRDS(file = paste0(path_results, "gls_2000_arc.rds"))
}
