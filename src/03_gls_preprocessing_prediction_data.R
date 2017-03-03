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
  gls_2000 <- readRDS(file = paste0(path_results, "gls_2000_arc.rds"))
  
  obsv_shp <- readRDS(file = paste0(path_results, "obsv_shp_arc.rds"))
  obsv_shp <- spTransform(obsv_shp, CRS(projection(gls_2000[[1]])))

  gls_2000_plots <- extractFromRasterSnips(raster = gls_2000, spatial = obsv_shp, 
                                           selector = "plot", buffer = 50)
  saveRDS(gls_2000_plots, file = paste0(path_results, "gls_2000_plots.rds"))
} else {
  gls_2000_plots <- readRDS(file = paste0(path_results, "gls_2000_plots.rds"))
}


# Prepare gpm data set used for remote sensing prediction study ----------------
if(compute){
  obsv_gls <- as.data.frame(gls_2000_plots)
  col_selector <- which(names(obsv_gls) == "plot")
    
  col_meta <- c(seq(which(names(obsv_gls) == "date"), 
                    which(names(obsv_gls) == "elevation")),
                seq(which(names(obsv_gls) == "rainfall"), 
                    which(names(obsv_gls) == "year")),
                seq(which(names(obsv_gls) == "ID"), 
                    which(names(obsv_gls) == "veg.nr")),
                seq(which(names(obsv_gls) == "coordW"), 
                    which(names(obsv_gls) == "coordN")))
    
  col_diversity <- c(which(names(obsv_gls) == "rich.plants"),
                     seq(which(names(obsv_gls) ==  "specrich"), 
                         which(names(obsv_gls) == "shannon.nw")))
    
  col_species <- seq(which(names(obsv_gls) == "Abisares.depressus"),
                     which(names(obsv_gls) == "Zonocerus.elegans"))
    
  col_precitors <- seq(which(names(obsv_gls) == "gls2000.1_mean"),
                       which(names(obsv_gls) == "ndvi_glcm_correlation.r25_var"))
    
  if((length(unique(c(col_selector, col_meta, col_diversity, 
                      col_species, col_precitors))) == ncol(obsv_gls))){
    obsv_gls[, col_species][obsv_gls[, col_species] > 0] <- "yes"
    obsv_gls[, col_species][obsv_gls[, col_species] == "0"] <- "no"
    for(i in col_species){
      obsv_gls[, i] <- as.factor(obsv_gls[, i])
    }
    meta <- createGPMMeta(obsv_gls, type = "input",
                          selector = col_selector, 
                          response = col_species, 
                          predictor = col_precitors, 
                          meta = c(col_meta, col_diversity))
    obsv_gls <- gpm(obsv_gls, meta, scale = TRUE)
  } else {
    obsv_gls <- NULL
  }


  saveRDS(obsv_gls, file = paste0(path_results, "obsv_gls_gpm.rds"))
} else {
  obsv_gls <- readRDS(file = paste0(path_results, "obsv_gls_gpm.rds"))
}


# Select responses occuring at least 25 plots on average -----------------------
if(compute){
  obsv_gls <- minimumOccurence(x = obsv_gls,
                               occurence = "yes", 
                               resample = 100, 
                               thv = 25)
  saveRDS(obsv_gls, file = paste0(path_results, "obsv_gls_gpm_minimumOccurence.rds"))
} else {
  obsv_gls <- readRDS(file = paste0(path_results, "obsv_gls_gpm_minimumOccurence.rds"))
}


# Clean predictor variables ----------------------------------------------------
if(compute){
  # Remove highly correlated or near zero variance predictors
  obsv_gls <- cleanPredictors(x = obsv_gls, nzv = TRUE, 
                              highcor = TRUE, cutoff = 0.80,
                              rmvna = TRUE)
  
  # Get some predictors back even if removed:
  # - all original bands
  # - all PCA band
  # - all VI bands
  # - all predictors computed with a windows radius of 1 back on track
  org_pca_vi_bands <- obsv_gls@meta$input$PREDICTOR[1:15][
    which(!(obsv_gls@meta$input$PREDICTOR[1:15] %in% obsv_gls@meta$input$PREDICTOR_FINAL))]
  obsv_gls@meta$input$PREDICTOR_FINAL <- c(obsv_gls@meta$input$PREDICTOR_FINAL, org_pca_vi_bands)
  
  r1 <- obsv_gls@meta$input$PREDICTOR[grep("r1_|r1o", obsv_gls@meta$input$PREDICTOR)]
  r1 <- r1[which(!(r1 %in% obsv_gls@meta$input$PREDICTOR_FINAL))]
  
  obsv_gls@meta$input$PREDICTOR_FINAL <- c(obsv_gls@meta$input$PREDICTOR_FINAL, r1)
  
  # Remove all additionally added predictors including NAs
  na_check <- colSums(is.na(obsv_gls@data$input[, obsv_gls@meta$input$PREDICTOR_FINAL]))
  narm <- which(obsv_gls@meta$input$PREDICTOR_FINAL %in%  names(na_check[na_check>0]))
  obsv_gls@meta$input$PREDICTOR_FINAL <- obsv_gls@meta$input$PREDICTOR_FINAL[-narm]

  saveRDS(obsv_gls, file = paste0(path_results, "obsv_gls_gpm_cleanPredictors.rds"))
} else {
  obsv_gls <- readRDS(file = paste0(path_results, "obsv_gls_gpm_cleanPredictors.rds"))
}


# Compile model training and evaluation dataset --------------------------------
if(compute){
    # Compute resamples
  obsv_gls <- resamplingsByVariable(x = obsv_gls,
                                    use_selector = TRUE,
                                    grabs = 1,
                                    resample = 100)
  
  # Split resamples into training and testing samples
  obsv_gls <- splitMultResp(x = obsv_gls, 
                            p = 0.85, 
                            use_selector = FALSE)
  saveRDS(obsv_gls, file = paste0(path_results, "obsv_gls_gpm_traintest.rds"))
} else {
  obsv_gls <- readRDS(file = paste0(path_results, "obsv_gls_gpm_traintest.rds"))
}
