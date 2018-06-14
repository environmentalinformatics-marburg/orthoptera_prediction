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
      modis <- readRDS(file = paste0(path_results, "mod_arc.rds"))
    } else {
      modis <- readRDS(file = paste0(path_results, "myd_arc.rds"))
    }
    obsv_shp <- spTransform(obsv_shp_arc, crs(modis[[1]]))
    
    modis_plots <- extractFromRasterSnips(raster = modis, spatial = obsv_shp, 
                                          selector = NULL, buffer = 850)
    saveRDS(modis_plots, file = paste0(path_results, sensor, "_plots.rds"))
  } 
} else {
  mod_plots <- readRDS(file = paste0(path_results, "mod_plots.rds"))
  myd_plots <- readRDS(file = paste0(path_results, "myd_plots.rds"))
}


# Prepare MODIS data set used for remote sensing prediction study -------------
if(compute){
  
  modis_sensors <- c("mod", "myd")
  for(sensor in modis_sensors){
    
    if(sensor == "mod"){
      modis_plots <- readRDS(file = paste0(path_results, "mod_plots.rds"))
    } else {
      modis_plots <- readRDS(file = paste0(path_results, "myd_plots.rds"))
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
      
      rm_pred = lapply(col_precitors, function(i){
        if(any(is.na(obsv_modis[, i]))){
          rm = i
        } else {
          rm = NULL
        }
        return(rm)
      })
      rm_pred = unlist(rm_pred)
      colnames(obsv_modis)[rm_pred]
      col_precitors = col_precitors[-which(col_precitors %in% rm_pred)]
      
      meta <- createGPMMeta(obsv_modis, type = "input",
                            selector = col_selector, 
                            response = col_species, 
                            predictor = col_precitors, 
                            meta = c(col_meta, col_diversity))
      obsv_modis <- gpm(obsv_modis, meta, scale = TRUE)
    } else {
      obsv_modis <- NULL
    }
    
    saveRDS(obsv_modis, file = paste0(path_results, sensor, "_gpm.rds"))
  }    
} else {
  obsv_mod <- readRDS(file = paste0(path_results, "mod_gpm.rds"))
  obsv_myd <- readRDS(file = paste0(path_results, "myd_gpm.rds"))
}


# Select responses occuring at least 25 plots on average -----------------------
if(compute){
  modis_sensors <- c("mod", "myd")
  for(sensor in modis_sensors){
    
    if(sensor == "mod"){
      obsv_modis <- readRDS(file = paste0(path_results, "mod_gpm.rds"))
    } else {
      obsv_modis <- readRDS(file = paste0(path_results, "myd_gpm.rds"))
    }
    
    obsv_modis <- minimumOccurence(x = obsv_modis,
                                   occurence = "yes", 
                                   resample = 100, 
                                   thv = 25)
    obsv_modis@meta$input$RESPONSE <- obsv_modis@meta$input$MIN_OCCURENCE$names
    
    saveRDS(obsv_modis, file = paste0(path_results, sensor, "_gpm_minimumOccurence.rds"))
  }
} else {
  obsv_mod <- readRDS(file = paste0(path_results, "mod_gpm_minimumOccurence.rds"))
  obsv_myd <- readRDS(file = paste0(path_results, "myd_gpm_minimumOccurence.rds"))
}


# Clean predictor variables ----------------------------------------------------
if(compute){
  
  modis_sensors <- c("mod", "myd")
  for(sensor in modis_sensors){
    
    if(sensor == "mod"){
      obsv_modis <- readRDS(file = paste0(path_results, "mod_gpm_minimumOccurence.rds"))
    } else {
      obsv_modis <- readRDS(file = paste0(path_results, "myd_gpm_minimumOccurence.rds"))
    }
    
    # Remove highly correlated or near zero variance predictors
    sel = c(grep(glob2rx("pca*b1r1o1m*mean"), 
                 obsv_modis@meta$input$PREDICTOR, value = TRUE),
            grep(glob2rx("pca*b1r1o1m*sd"), 
                 obsv_modis@meta$input$PREDICTOR, value = TRUE))
    # sel =
    #   c("gls2000.1_mean", "gls2000.2_mean", "gls2000.3_mean", "gls2000.4_mean",
    #     "gls2000.5_mean", "gls2000.6_mean", "NDVI_mean", "GNDVI_mean", "SR_mean",
    #     sel)
    
    obsv_modis@meta$input$PREDICTOR_FINAL =  sel
    
    obsv_modis <- cleanPredictors(x = obsv_modis, nzv = TRUE, 
                                  highcor = TRUE, cutoff = 0.80,
                                  rmvna = TRUE)
    obsv_modis@meta$input$PREDICTOR_FINAL =    
      c("modis_sur_refl_b01_mean", "modis_sur_refl_b02_mean", "modis_sur_refl_b03_mean", "modis_sur_refl_b04_mean",
        "modis_sur_refl_b05_mean", "modis_sur_refl_b07_mean", "NDVI_mean", "GNDVI_mean", "SR_mean",
        obsv_modis@meta$input$PREDICTOR_FINAL)
    
    obsv_modis@meta$input$PREDICTOR_FINAL
    
    # obsv_modis <- cleanPredictors(x = obsv_modis, nzv = TRUE, 
    #                               highcor = TRUE, cutoff = 0.75)
    
    
    # Get some predictors back even if removed:
    # - all original bands
    # - all PCA band
    # - all VI bands
    # - all predictors computed with a windows radius of 1 back on track
    # org_pca_vi_bands <- obsv_modis@meta$input$PREDICTOR[1:17][
    #   which(!(obsv_modis@meta$input$PREDICTOR[1:17] %in% obsv_modis@meta$input$PREDICTOR_FINAL))]
    # obsv_modis@meta$input$PREDICTOR_FINAL <- c(obsv_modis@meta$input$PREDICTOR_FINAL, org_pca_vi_bands)
    # 
    # r1 <- obsv_modis@meta$input$PREDICTOR[grep("r1_|r1o", obsv_modis@meta$input$PREDICTOR)]
    # r1 <- r1[which(!(r1 %in% obsv_modis@meta$input$PREDICTOR_FINAL))]
    # 
    # obsv_modis@meta$input$PREDICTOR_FINAL <- c(obsv_modis@meta$input$PREDICTOR_FINAL, r1)
    # 
    # # Remove all additionally added predictors including NAs
    # na_check <- colSums(is.na(obsv_modis@data$input[, obsv_modis@meta$input$PREDICTOR_FINAL]))
    # narm <- which(obsv_modis@meta$input$PREDICTOR_FINAL %in%  names(na_check[na_check>0]))
    # obsv_modis@meta$input$PREDICTOR_FINAL <- obsv_modis@meta$input$PREDICTOR_FINAL[-narm]
    
    saveRDS(obsv_modis, file = paste0(path_results, sensor, "_gpm_cleanPredictors.rds"))
  }
} else {
  obsv_mod <- readRDS(file = paste0(path_results, "mod_gpm_cleanPredictors.rds"))
  obsv_myd <- readRDS(file = paste0(path_results, "myd_gpm_cleanPredictors.rds"))
}

# Compile model training and evaluation dataset --------------------------------
if(compute){
  
  modis_sensors <- c("mod", "myd")
  for(sensor in modis_sensors){
    
    if(sensor == "mod"){
      obsv_modis <- readRDS(file = paste0(path_results, "mod_gpm_cleanPredictors.rds"))
    } else {
      obsv_modis <- readRDS(file = paste0(path_results, "myd_gpm_cleanPredictors.rds"))
    }
    
    # Compute resamples
    obsv_modis <- resamplingsByVariable(x = obsv_modis,
                                        use_selector = TRUE,
                                        grabs = 1,
                                        resample = 25)
    
    # Split resamples into training and testing samples
    obsv_modis <- splitMultResp(x = obsv_modis, 
                                p = 0.85, 
                                use_selector = FALSE)
    
    saveRDS(obsv_modis, file = paste0(path_results, sensor, "_gpm_traintest.rds"))
  }
  
} else {
  obsv_mod <- readRDS(file = paste0(path_results, "mod_gpm_traintest.rds"))
  obsv_myd <- readRDS(file = paste0(path_results, "myd_gpm_traintest.rds"))
}
