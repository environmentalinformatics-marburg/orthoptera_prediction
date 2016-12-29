# Preprocess MODIS dataset
# Thomas Nauss

if(Sys.info()["sysname"] == "Windows"){
  source("F:/analysis/orthoptera/orthoptera_prediction/src/00_set_environment.R")
} else {
  source("/media/tnauss/myWork/analysis/orthoptera/orthoptera_prediction/src/00_set_environment.R")
}

initOTB("C:/OSGeo4W64/bin/")
download <- FALSE
compute <- FALSE


# Download MODIS tiles ---------------------------------------------------------
if(download){
  # Read observation data for dates
  obsv <- read.table(paste0(path_obsv, "lvl0300_biodiversity_data.csv"), 
                     header = TRUE, sep = ";", dec = ",")
  dates <- gsub("-", "", substring(unique(obsv$date_nocloud), 1, 8), "-")
  
  # Set processing options
  MODISoptions(localArcPath = path_modis_local_arc, outProj = "+init=epsg:32737",
               outDirPath = path_modis_out_dir, MODISserverOrder = c("LAADS", "LPDAAC"),
               gdalPath = "C:/OSGeo4W64/bin")
  
  # Download MODIS tiles
  modis_product <- c("MOD09GA", "MOD09A1", "MOD09Q1", "MOD11A1",
                     "MYD09GA", "MYD09A1", "MYD09Q1", "MYD11A1")
  
  for(pdct in modis_product){
    for(dts in dates){
      job_name <- paste0(pdct, "_TZS")
      collection <- getCollection(pdct, forceCheck = TRUE)
      runGdal(pdct, collection = collection, job = job_name, 
              tileH = 21, tileV = 9,
              begin = dts, end = dts)
    }
  }
}

# Prepare MODIS dataset --------------------------------------------------------
# modis_product <- c("MOD09GA", "MOD09A1", "MOD09Q1",
#                    "MYD09GA", "MYD09A1", "MYD09Q1")
modis_sensors <- c("mod", "myd")

tmin <- as.POSIXct(strptime("2002-01-01", "%Y-%m-%d"), tz = "UTC")
tmax <- as.POSIXct(strptime("2013-01-01", "%Y-%m-%d"), tz = "UTC")

# Create one raster tile for each observation plot
obsv_shp_wgs <- readRDS(file = paste0(path_results, "obsv_shp_wgs.rds"))
obsv_shp_arc <- readRDS(file = paste0(path_results, "obsv_shp_arc.rds"))

obsv_shp_wgs_modis <- obsv_shp_wgs[obsv_shp_wgs@data$date >= tmin & 
                                     obsv_shp_wgs@data$date <= tmax,]
obsv_shp_arc_modis <- obsv_shp_wgs[obsv_shp_arc@data$date >= tmin & 
                                     obsv_shp_arc@data$date <= tmax,]

saveRDS(obsv_shp_wgs_modis, file = paste0(path_results, "obsv_shp_wgs_modis.rds"))
saveRDS(obsv_shp_arc_modis, file = paste0(path_results, "obsv_shp_arc_modis.rds"))

for(sensor in modis_sensors){
  
  if(sensor == "mod"){
    prdct <- "MOD09A1"
  } else {
    prdct <- "MYD09A1"
  }
  
  job_name <- paste0(pdct, "_TZS")
  modis_files <- list.files(paste0(path_modis_out_dir, job_name),
                            pattern = glob2rx("*sur_refl_b0*tif"), 
                            full.names = TRUE)
  
  for(prj in c("wgs", "arc")){
    if(prj == "wgs"){
      obsv_shp <- obsv_shp_wgs_modis
    } else {
      obsv_shp <- obsv_shp_arc_modis
    }
    obsv_shp <- spTransform(obsv_shp, crs(stack(modis_files[1])))
    
    time_match <- timeMatch(a = obsv_shp@data$date, 
                            b = unique(as.POSIXct(
                              strptime(substr(basename(modis_files), 10, 16),
                                       "%Y%j"), tz = "UTC")))
    
    modis_snip <- lapply(seq(nrow(time_match)), function(i){
      if(i %% 10 == 0) print(paste0("Processing ", i, " of ", nrow(time_match)))
      snipRaster(raster=stack(modis_files[grep(strftime(time_match$b[i], "%Y%j"), 
                                               modis_files)]), 
                 spatial=obsv_shp[i,], selector = NULL,
                 buffer=4500, byid = TRUE)
    })
    saveRDS(modis_snip, file = paste0(path_results, "modis_", sensor, "_snip_", prj, ".rds"))
  }
}

    
    # Compute pca for all plots
    modis_pca <- lapply(modis_snip, function(t){
      pca(stack(t), return_raster = TRUE)
    })
    saveRDS(modis_pca, file = paste0(path_results, "modis_pca_", prj, ".rds"))
    
    # Compute spectral indices for all plots
    mspec_indices <- lapply(seq(length(modis_snip)), function(i){
      if(i %% 10 == 0) print(paste0(i))
      mSpecIndices(blue = modis_snip[[i]][[1]], green = modis_snip[[i]][[2]], 
                   red = modis_snip[[i]][[3]], nir = modis_snip[[i]][[4]])
    })  
    saveRDS(mspec_indices, file = paste0(path_results, "mspec_indices_", prj, ".rds"))
    
    # Compute Haralick textures for all plots
    otb_txt <- lapply(seq(length(modis_pca)), function(i){
      if(i %% 10 == 0) print(paste0(i))
      otbTexturesHaralick(x=modis_pca[[i]], path_output = path_temp, 
                          return_raster = TRUE, 
                          parameters.xyrad=list(c(1,1)),
                          parameters.xyoff=list(c(1,1)),
                          texture="all",
                          channel = 1)
    })
    saveRDS(otb_txt, file = paste0(path_results, "otb_txt_", prj, ".rds"))
    
    # Compute glcm textures for all plots
    glcm_txt <- lapply(seq(length(modis_snip)), function(i){
      if(i %% 10 == 0) print(paste0(i))
      glcm(modis_pca[[i]][[1]], n_grey = 32, window = c(3,3),
           shift=list(c(0,1), c(1,1), c(1,0), c(1,-1)))
    })  
    
    saveRDS(glcm_txt, file = paste0(path_results, "glcm_txt_", prj, ".rds"))
    # Combine results in one stack per tile
    modis_2000 <- lapply(seq(length(modis_snip)), function(i){
      stack(modis_snip[[i]], modis_pca[[i]], mspec_indices[[i]],
            otb_txt[[i]], glcm_txt[[i]])
    })
    names(modis_2000) <- names(modis_snip)
    saveRDS(modis_2000, file = paste0(path_results, "modis_2000_", prj, ".rds"))
    
    
    
}


#TODO

# Prepare MODIS dataset ------------------------------------------------------
if(compute){
  gls <- stack(paste0(path_landsat, "gls2000.tif"))
  
    
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



