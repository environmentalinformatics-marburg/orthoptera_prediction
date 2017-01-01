# Preprocess dataset by combining field and satellite observations
# Thomas Nauss

if(Sys.info()["sysname"] == "Windows"){
  source("F:/analysis/orthoptera/orthoptera_prediction/src/00_set_environment.R")
} else {
  source("/media/tnauss/myWork/analysis/orthoptera/orthoptera_prediction/src/00_set_environment.R")
}

compute <- TRUE

# Merge GLS2000 data with preprocessed orthoptera observations -----------------
if(compute){
  modis_sensors <- c("mod", "myd")
  
  obsv_shp_wgs <- readRDS(file = paste0(path_results, "obsv_shp_wgs_modis.rds"))
  obsv_shp_arc <- readRDS(file = paste0(path_results, "obsv_shp_arc_modis.rds"))
  
  for(sensor in modis_sensors){
    
    if(sensor == "mod"){
      modis_wgs <- readRDS(file = paste0(path_results, "modis_mod_wgs.rds"))
      modis_arc <- readRDS(file = paste0(path_results, "modis_mod_arc.rds"))
    } else {
      modis_wgs <- readRDS(file = paste0(path_results, "modis_myd_wgs.rds"))
      modis_arc <- readRDS(file = paste0(path_results, "modis_myd_arc.rds"))
    }
    cprj <- c("wgs", "arc")
    modis_plots <- lapply(cprj, function(prj){
      if(prj == "wgs"){
        modis <- modis_wgs
        obsv_shp <- spTransform(obsv_shp_wgs, crs(modis[[1]]))
      } else {
        modis <- modis_arc
        obsv_shp <- spTransform(obsv_shp_arc, crs(modis[[1]]))
      }
      
#       obsv_shp@data$date_double <- as.double(obsv_shp@data$date)
#       names(modis) <- as.double(as.POSIXct(strptime(names(modis), "%Y-%m-%d"), tz = "UTC"))

      plots <- extractFromRasterSnips(raster = modis, spatial = obsv_shp, 
                                      selector = NULL, buffer = 850)
    })
    names(modis_plots) <- cprj
    saveRDS(modis_plots, file = paste0(path_results, "modis", sensor, "plots.rds"))
  } 
} else {
  modis_mod_plots <- readRDS(file = paste0(path_results, "modis_mod_plots.rds"))
  modis_myd_plots <- readRDS(file = paste0(path_results, "modis_myd_plots.rds"))
}
  

# Prepare gpm data set used for remote sensing prediction study ----------------
if(compute){

  modis_sensors <- c("mod", "myd")
  for(sensor in modis_sensors){
    
    if(sensor == "mod"){
      modis_plots <- readRDS(file = paste0(path_results, "modis_mod_plots.rds"))
    } else {
      modis_plots <- readRDS(file = paste0(path_results, "modis_myd_plots.rds"))
    }
    cprj <- c("wgs", "arc")
    obsv_modis <- lapply(cprj, function(prj){
      obsv <- as.data.frame(modis_plots[[prj]])
      col_selector <- which(names(obsv) == "plot")
      
      col_meta <- c(seq(which(names(obsv) == "date"), 
                        which(names(obsv) == "elevation")),
                    seq(which(names(obsv) == "rainfall"), 
                        which(names(obsv) == "year")),
                    seq(which(names(obsv) == "ID"), 
                        which(names(obsv) == "veg.nr")),
                    seq(which(names(obsv) == "coordW"), 
                        which(names(obsv) == "coordN")))
      
      col_diversity <- c(which(names(obsv) == "rich.plants"),
                         seq(which(names(obsv) ==  "specrich"), 
                             which(names(obsv) == "shannon.nw")))
      
      col_species <- seq(which(names(obsv) == "Abisares.depressus"),
                         which(names(obsv) == "Zonocerus.elegans"))
      
      col_precitors <- seq(which(names(obsv) == "gls2000.1"),
                           which(names(obsv) == "ndvi_glcm_correlation_sd"))
      
      if((length(unique(c(col_selector, col_meta, col_diversity, 
                          col_species, col_precitors))) == ncol(obsv))){
        obsv[, col_species][obsv[, col_species] > 0] <- "yes"
        obsv[, col_species][obsv[, col_species] == "0"] <- "no"
        for(i in col_species){
          obsv[, i] <- as.factor(obsv[, i])
        }
        meta <- createGPMMeta(obsv, type = "input",
                              selector = col_selector, 
                              response = col_species, 
                              independent = col_precitors, 
                              meta = c(col_meta, col_diversity))
        obsv <- gpm(obsv, meta, scale = TRUE)
      } else {
        obsv <- NULL
      }
      return(obsv)
    })
    names(obsv_modis) <- cprj
    
    saveRDS(obsv_modis, file = paste0(path_results, "obsv_modis_", sensor, ".rds"))
    saveRDS(obsv_modis[["wgs"]], file = paste0(path_results, "obsv_modis_", sensor, "_wgs.rds"))
    saveRDS(obsv_modis[["arc"]], file = paste0(path_results, "obsv_modis_", sensor, "_arc.rds"))
  }    
} else {
  obsv_modis_mod <- readRDS(file = paste0(path_results, "obsv_modis_mod.rds"))
  obsv_modis_myd <- readRDS(file = paste0(path_results, "obsv_modis_myd.rds"))
}


# Select responses occuring at least 20 plots on average -----------------------
if(compute){
  
  modis_sensors <- c("mod", "myd")
  for(sensor in modis_sensors){
    
    if(sensor == "mod"){
      obsv_modis <- readRDS(file = paste0(path_results, "obsv_modis_mod.rds"))
    } else {
      obsv_modis <- readRDS(file = paste0(path_results, "obsv_modis_myd.rds"))
    }

    cprj <- c("wgs", "arc")
    obsv_modis <- lapply(cprj, function(prj){
      obsv <- obsv_modis[[prj]]
      obsv <- minimumOccurence(x = obsv,
                               occurence = "yes", 
                               resample = 100, 
                               thv = 25)
      obsv@meta$input$RESPONSE <- obsv@meta$input$MIN_OCCURENCE$names
      return(obsv)
    })
    names(obsv_modis) <- cprj
    saveRDS(obsv_modis, file = paste0(path_results, "obsv_modis_", sensor, "_minimumOccurence.rds"))
  }
  } else {
    obsv_modis_mod <- readRDS(file = paste0(path_results, "obsv_modis_mod_minimumOccurence.rds"))
    obsv_modis_myd <- readRDS(file = paste0(path_results, "obsv_modis_myd_minimumOccurence.rds"))
}



# Compile model training and evaluation dataset --------------------------------
if(compute){
  
  modis_sensors <- c("mod", "myd")
  for(sensor in modis_sensors){
    
    if(sensor == "mod"){
      obsv_modis <- readRDS(file = paste0(path_results, "obsv_modis_mod_minimumOccurence.rds"))
    } else {
      obsv_modis <- readRDS(file = paste0(path_results, "obsv_modis_myd_minimumOccurence.rds"))
    }
    cprj <- c("wgs", "arc")
    obsv_modis <- lapply(cprj, function(prj){
      obsv <- obsv_modis[[prj]]
      
      # Compute resamples
      obsv <- resamplingsByVariable(x = obsv,
                                    use_selector = TRUE,
                                    grabs = 1,
                                    resample = 100)
      
      # Split resamples into training and testing samples
      obsv <- splitMultResp(x = obsv, 
                            p = 0.75, 
                            use_selector = TRUE)
      
    })
    names(obsv_modis) <- cprj
    saveRDS(obsv_modis, file = paste0(path_results, "obsv_modis_", sensor, "_traintest.rds"))
  }
  
} else {
  obsv_modis_mod <- readRDS(file = paste0(path_results, "obsv_modis_mod_traintest.rds"))
  obsv_modis_myd <- readRDS(file = paste0(path_results, "obsv_modis_myd_traintest.rds"))
}
