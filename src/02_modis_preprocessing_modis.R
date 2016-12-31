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
if(compute){
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
      pdct <- "MOD09A1"
    } else {
      pdct <- "MYD09A1"
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
      
      # Create one raster tile for each observation plot
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
      modis_snip <- unlist(modis_snip)
      names(modis_snip) <- time_match$a
      if(sensor == "myd"){
        modis_snip <- lapply(modis_snip, function(s){
          s[[-6]]
        })
      }
      saveRDS(modis_snip, file = paste0(path_results, "modis_", sensor, "_snip_", prj, ".rds"))
      
      # Compute pca for all plots
      modis_pca <- pca(modis_snip, ignore_names = TRUE, 
                       center = TRUE, scale = TRUE)
      saveRDS(modis_pca, file = paste0(path_results, "modis_", sensor, "_pca_", prj, ".rds"))
      
      # Compute spectral indices for all plots
      modis_mspec_indices <- lapply(seq(length(modis_snip)), function(i){
        if(i %% 10 == 0) print(paste0(i))
        mSpecIndices(blue = modis_snip[[i]][[3]], green = modis_snip[[i]][[4]], 
                     red = modis_snip[[i]][[1]], nir = modis_snip[[i]][[2]])
      })  
      names(modis_mspec_indices) <- names(modis_snip)
      saveRDS(modis_mspec_indices, file = paste0(path_results, "modis_", sensor, "_mspec_indices_", prj, ".rds"))
      
      # Compute Haralick textures for all plots based on PCA
      minv <- min(unlist(lapply(modis_pca, function(s){minValue(s$PC1)})))
      maxv <- max(unlist(lapply(modis_pca, function(s){maxValue(s$PC1)})))
      modis_pca_otb_txt <- lapply(seq(length(modis_pca)), function(i){
        if(i %% 10 == 0) print(paste0(i))
        oth <- otbTexturesHaralick(x=modis_pca[[i]]$PC1, path_output = path_temp, 
                                   return_raster = TRUE, 
                                   parameters.xyrad=list(c(1,1)),
                                   parameters.xyoff=list(c(1,1)),
                                   parameters.minmax=c(minv, maxv),
                                   parameters.nbbin = 8,
                                   texture="all",
                                   channel = 1)
        names(oth) <- paste0("pca_", names(oth))
        return(oth)
      })
      names(modis_pca_otb_txt) <- names(modis_snip)
      saveRDS(modis_pca_otb_txt, file = paste0(path_results, "modis_", sensor, "_pca_otb_txt_", prj, ".rds"))
      
      # Compute Haralick textures for all plots based on NDVI
      minv <- -1
      maxv <- 1
      modis_ndvi_otb_txt <- lapply(seq(length(modis_mspec_indices)), function(i){
        if(i %% 10 == 0) print(paste0(i))
        oth <- otbTexturesHaralick(x=modis_mspec_indices[[i]]$NDVI, path_output = path_temp, 
                                   return_raster = TRUE, 
                                   parameters.xyrad=list(c(1,1)),
                                   parameters.xyoff=list(c(1,1)),
                                   parameters.minmax=c(minv, maxv),
                                   parameters.nbbin = 8,
                                   texture="all",
                                   channel = 1)
        names(oth) <- paste0("ndvi_", names(oth))
        return(oth)
      })
      names(modis_ndvi_otb_txt) <- names(modis_snip)
      saveRDS(modis_ndvi_otb_txt, file = paste0(path_results, "modis_", sensor, "_ndvi_otb_txt_", prj, ".rds"))
      
      # Compute glcm textures for all plots based on PCA
      modis_pca_glcm_txt <- lapply(seq(length(modis_pca)), function(i){
        if(i %% 10 == 0) print(paste0(i))
        gt <- glcm(modis_pca[[i]]$PC1, n_grey = 32, window = c(3,3),
                   shift=list(c(0,1), c(1,1), c(1,0), c(1,-1)))
        names(gt) <- paste0("pca_", names(gt))
        return(gt)
      })  
      names(modis_pca_glcm_txt) <- names(modis_snip)
      saveRDS(modis_pca_glcm_txt, file = paste0(path_results, "modis_", sensor, "_pca_glcm_txt_", prj, ".rds"))
      
      # Compute glcm textures for all plots based on NDVI
      modis_ndvi_glcm_txt <- lapply(seq(length(modis_mspec_indices)), function(i){
        if(i %% 10 == 0) print(paste0(i))
        gt <- glcm(modis_mspec_indices[[i]]$NDVI, n_grey = 32, window = c(3,3),
                   shift=list(c(0,1), c(1,1), c(1,0), c(1,-1)))
        names(gt) <- paste0("ndvi_", names(gt))
        return(gt)
      })
      names(modis_ndvi_glcm_txt) <- names(modis_snip)
      saveRDS(modis_ndvi_glcm_txt, file = paste0(path_results, "modis_", sensor, "_ndvi_glcm_txt_", prj, ".rds"))
      
      # Combine results in one stack per tile
      modis <- lapply(seq(length(modis_snip)), function(i){
        stack(modis_snip[[i]], modis_pca[[i]], modis_mspec_indices[[i]],
              modis_pca_otb_txt[[i]], modis_ndvi_otb_txt[[i]],
              modis_pca_glcm_txt[[i]], modis_ndvi_glcm_txt[[i]])
      })
      names(modis) <- names(modis_snip)
      saveRDS(modis, file = paste0(path_results, "modis_", sensor, "_", prj, ".rds"))
    }
  }
} else {
  sensor = "myd"
  prj = "arc"
  modis_snip <- readRDS(file = paste0(path_results, "modis_", sensor, "_snip_", prj, ".rds"))
  modis_pca <- readRDS(file = paste0(path_results, "modis_", sensor, "_pca_", prj, ".rds"))
  modis_mspec_indices <- readRDS(file = paste0(path_results, "modis_", sensor, "_mspec_indices_", prj, ".rds"))
  modis_pca_otb_txt <- readRDS(file = paste0(path_results, "modis_", sensor, "_pca_otb_txt_", prj, ".rds"))
  modis_ndvi_otb_txt <- readRDS(file = paste0(path_results, "modis_", sensor, "_ndvi_otb_txt_", prj, ".rds"))
  modis_pca_glcm_txt <- readRDS(file = paste0(path_results, "modis_", sensor, "_pca_glcm_txt_", prj, ".rds"))
  modis_ndvi_glcm_txt <- readRDS(file = paste0(path_results, "modis_", sensor, "_ndvi_glcm_txt_", prj, ".rds"))
  
  modis_mod_wgs <- readRDS(file = paste0(path_results, "modis_mod_wgs.rds"))
  modis_mod_arc <- readRDS(file = paste0(path_results, "modis_mod_arc.rds"))
  
  modis_myd_wgs <- readRDS(file = paste0(path_results, "modis_myd_wgs.rds"))
  modis_myd_arc <- readRDS(file = paste0(path_results, "modis_myd_arc.rds"))
  
}
