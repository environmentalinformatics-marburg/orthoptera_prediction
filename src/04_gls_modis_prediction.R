# Predict orthoptera based on satellite observations
# Thomas Nauss

if(Sys.info()["sysname"] == "Windows"){
  source("F:/analysis/orthoptera/orthoptera_prediction/src/00_set_environment.R")
} else {
  source("/media/tnauss/myWork/analysis/orthoptera/orthoptera_prediction/src/00_set_environment.R")
}

compute <- TRUE

# Predict dataset --------------------------------------------------------------
if(compute){
  obsv_gpm <- readRDS(file = paste0(path_results, "obsv_gls_gpm_traintest.rds"))
  n_var <- c(seq(1, 29), seq(40, length(obsv_gpm@meta$input$PREDICTOR_FINAL), 40))
  
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  
  obsv_gpm <- trainModel(x = obsv_gpm,
                         n_var = n_var, 
                         mthd = "rf",
                         mode = "rfe",
                         seed_nbr = 11, 
                         cv_nbr = 5,
                         var_selection = "indv", 
                         metric = "Kappa",
                         response_nbr = c(1, 5, 10),
                         resample_nbr = c(1:10),
                         filepath_tmp = path_temp)
  saveRDS(obsv_gpm, file = paste0(path_results, "obsv_gls_gpm_trainmodel.rds"))
  
  
  
  obsv_modis_mod <- readRDS(file = paste0(path_results, "obsv_modis_gpm_mod_traintest.rds"))
  
  n_var <- c(seq(1, 29), seq(40, length(obsv_modis_mod@meta$input$PREDICTOR_FINAL), 40))
  
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  
  obsv_modis_mod <- trainModel(x = obsv_modis_mod,
                               n_var = NULL, 
                               mthd = "rf",
                               mode = "rfe",
                               seed_nbr = 11, 
                               cv_nbr = 5,
                               var_selection = "indv", 
                               filepath_tmp = NULL)
  saveRDS(obsv_modis_mod, file = paste0(path_results, "obsv_mod_gpm_trainmodel.rds"))
  
  
  
  obsv_modis_myd <- readRDS(file = paste0(path_results, "obsv_modis_gpm_myd_traintest.rds"))
  
  n_var <- c(seq(1, 29), seq(40, length(obsv_modis_myd@meta$input$PREDICTOR_FINAL), 40))
  
  obsv_modis_myd <- trainModel(x = obsv_modis_myd,
                               n_var = NULL, 
                               mthd = "rf",
                               mode = "rfe",
                               seed_nbr = 11, 
                               cv_nbr = 5,
                               var_selection = "indv", 
                               filepath_tmp = NULL)
  saveRDS(obsv_modis_myd, file = paste0(path_results, "obsv_myd_gpm_trainmodel.rds"))
}
