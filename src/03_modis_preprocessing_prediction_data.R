# Preprocess dataset by combining field and satellite observations
# Thomas Nauss

if(Sys.info()["sysname"] == "Windows"){
  source("F:/analysis/orthoptera/orthoptera_prediction/src/00_set_environment.R")
} else {
  source("/media/tnauss/myWork/analysis/orthoptera/orthoptera_prediction/src/00_set_environment.R")
}

compute <- TRUE

# Prepare MODIS data ----------------------------------------------------------
if(compute){
  modis_sensors <- c("mod", "myd")
  
  obsv_shp_arc <- readRDS(file = paste0(path_results, "obsv_shp_arc_modis.rds"))
  
  for(sensor in modis_sensors){
    
    if(sensor == "mod"){
      modis <- readRDS(file = paste0(path_results, "modis_mod_arc.rds"))
    } else {
      modis <- readRDS(file = paste0(path_results, "modis_myd_arc.rds"))
    }
    obsv_shp <- spTransform(obsv_shp_arc, crs(modis[[1]]))

    modis_plots <- extractFromRasterSnips(raster = modis, spatial = obsv_shp, 
                                          selector = NULL, buffer = 850)
    saveRDS(modis_plots, file = paste0(path_results, "modis_", sensor, "_plots.rds"))
  } 
} else {
  modis_mod_plots <- readRDS(file = paste0(path_results, "modis_mod_plots.rds"))
  modis_myd_plots <- readRDS(file = paste0(path_results, "modis_myd_plots.rds"))
}
  

# Prepare MODIS data set used for remote sensing prediction study -------------
if(compute){

  modis_sensors <- c("mod", "myd")
  for(sensor in modis_sensors){
    
    if(sensor == "mod"){
      modis_plots <- readRDS(file = paste0(path_results, "modis_mod_plots.rds"))
    } else {
      modis_plots <- readRDS(file = paste0(path_results, "modis_myd_plots.rds"))
    }
    obsv_modis <- as.data.frame(modis_plots)
      col_selector <- which(names(obsv_modis) == "plot")
      
      col_meta <- c(seq(which(names(obsv_modis) == "date"), 
                        which(names(obsv_modis) == "elevation")),
                    seq(which(names(obsv_modis) == "rainfall"), 
                        which(names(obsv_modis) == "year")),
                    seq(which(names(obsv_modis) == "ID"), 
                        which(names(obsv_modis) == "veg.nr")),
                    seq(which(names(obsv_modis) == "coordW"), 
                        which(names(obsv_modis) == "coordN")))
      
      col_diversity <- c(which(names(obsv_modis) == "rich.plants"),
                         seq(which(names(obsv_modis) ==  "specrich"), 
                             which(names(obsv_modis) == "shannon.nw")))
      
      col_species <- seq(which(names(obsv_modis) == "Abisares.depressus"),
                         which(names(obsv_modis) == "Zonocerus.elegans"))
      
      col_precitors <- seq(which(names(obsv_modis) == "modis_sur_refl_b01_mean"),
                           which(names(obsv_modis) == "ndvi_glcm_correlation_var"))
      
      if((length(unique(c(col_selector, col_meta, col_diversity, 
                          col_species, col_precitors))) == ncol(obsv_modis))){
        obsv_modis[, col_species][obsv_modis[, col_species] > 0] <- "yes"
        obsv_modis[, col_species][obsv_modis[, col_species] == "0"] <- "no"
        for(i in col_species){
          obsv_modis[, i] <- as.factor(obsv_modis[, i])
        }
        meta <- createGPMMeta(obsv_modis, type = "input",
                              selector = col_selector, 
                              response = col_species, 
                              predictor = col_precitors, 
                              meta = c(col_meta, col_diversity))
        obsv_modis <- gpm(obsv_modis, meta, scale = TRUE)
      } else {
        obsv_modis <- NULL
      }

    saveRDS(obsv_modis, file = paste0(path_results, "obsv_modis_gpm_", sensor, ".rds"))
  }    
} else {
  obsv_modis_mod <- readRDS(file = paste0(path_results, "obsv_modis_gpm_mod.rds"))
  obsv_modis_myd <- readRDS(file = paste0(path_results, "obsv_modis_gpm_myd.rds"))
}


# Select responses occuring at least 25 plots on average -----------------------
if(compute){
  modis_sensors <- c("mod", "myd")
  for(sensor in modis_sensors){
    
    if(sensor == "mod"){
      obsv_modis <- readRDS(file = paste0(path_results, "obsv_modis_gpm_mod.rds"))
    } else {
      obsv_modis <- readRDS(file = paste0(path_results, "obsv_modis_gpm_myd.rds"))
    }

    obsv_modis <- minimumOccurence(x = obsv_modis,
                               occurence = "yes", 
                               resample = 100, 
                               thv = 25)
    obsv_modis@meta$input$RESPONSE <- obsv_modis@meta$input$MIN_OCCURENCE$names

    saveRDS(obsv_modis, file = paste0(path_results, "obsv_modis_gpm_", sensor, "_minimumOccurence.rds"))
  }
  } else {
    obsv_modis_mod <- readRDS(file = paste0(path_results, "obsv_modis_gpm_mod_minimumOccurence.rds"))
    obsv_modis_myd <- readRDS(file = paste0(path_results, "obsv_modis_gpm_myd_minimumOccurence.rds"))
}


# Clean predictor variables ----------------------------------------------------
if(compute){
  
  modis_sensors <- c("mod", "myd")
  for(sensor in modis_sensors){
    
    if(sensor == "mod"){
      obsv_modis <- readRDS(file = paste0(path_results, "obsv_modis_gpm_mod_minimumOccurence.rds"))
    } else {
      obsv_modis <- readRDS(file = paste0(path_results, "obsv_modis_gpm_myd_minimumOccurence.rds"))
    }
    
    obsv_modis <- cleanPredictors(x = obsv_modis, nzv = TRUE, 
                                  highcor = TRUE, cutoff = 0.90)

    
    # Get some predictors back even if removed:
    # - all original bands
    # - all PCA band
    # - all VI bands
    # - all predictors computed with a windows radius of 1 back on track
    org_pca_vi_bands <- obsv_modis@meta$input$PREDICTOR[1:17][
      which(!(obsv_modis@meta$input$PREDICTOR[1:17] %in% obsv_modis@meta$input$PREDICTOR_FINAL))]
    obsv_modis@meta$input$PREDICTOR_FINAL <- c(obsv_modis@meta$input$PREDICTOR_FINAL, org_pca_vi_bands)
    
    r1 <- obsv_modis@meta$input$PREDICTOR[grep("r1_|r1o", obsv_modis@meta$input$PREDICTOR)]
    r1 <- r1[which(!(r1 %in% obsv_modis@meta$input$PREDICTOR_FINAL))]
    
    obsv_modis@meta$input$PREDICTOR_FINAL <- c(obsv_modis@meta$input$PREDICTOR_FINAL, r1)
    
    # Remove all additionally added predictors including NAs
    na_check <- colSums(is.na(obsv_modis@data$input[, obsv_modis@meta$input$PREDICTOR_FINAL]))
    narm <- which(obsv_modis@meta$input$PREDICTOR_FINAL %in%  names(na_check[na_check>0]))
    obsv_modis@meta$input$PREDICTOR_FINAL <- obsv_modis@meta$input$PREDICTOR_FINAL[-narm]
    
    saveRDS(obsv_modis, file = paste0(path_results, "obsv_modis_gpm_", sensor, "_cleanPredictors.rds"))
  }
} else {
  obsv_modis_mod <- readRDS(file = paste0(path_results, "obsv_modis_gpm_mod_cleanPredictors.rds"))
  obsv_modis_myd <- readRDS(file = paste0(path_results, "obsv_modis_gpm_myd_cleanPredictors.rds"))
}

# Compile model training and evaluation dataset --------------------------------
if(compute){
  
  modis_sensors <- c("mod", "myd")
  for(sensor in modis_sensors){
    
    if(sensor == "mod"){
      obsv_modis <- readRDS(file = paste0(path_results, "obsv_modis_gpm_mod_cleanPredictors.rds"))
    } else {
      obsv_modis <- readRDS(file = paste0(path_results, "obsv_modis_gpm_myd_cleanPredictors.rds"))
    }

    # Compute resamples
    obsv_modis <- resamplingsByVariable(x = obsv_modis,
                                        use_selector = TRUE,
                                        grabs = 1,
                                        resample = 100)
    
    # Split resamples into training and testing samples
    obsv_modis <- splitMultResp(x = obsv_modis, 
                                p = 0.85, 
                                use_selector = FALSE)
      
    saveRDS(obsv_modis, file = paste0(path_results, "obsv_modis_gpm_", sensor, "_traintest.rds"))
  }
  
} else {
  obsv_modis_mod <- readRDS(file = paste0(path_results, "obsv_modis_gpm_mod_traintest.rds"))
  obsv_modis_myd <- readRDS(file = paste0(path_results, "obsv_modis_gpm_myd_traintest.rds"))
}
