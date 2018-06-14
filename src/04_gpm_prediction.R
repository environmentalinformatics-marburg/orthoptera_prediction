# Predict orthoptera based on satellite observations
# Thomas Nauss

if(Sys.info()["sysname"] == "Windows"){
  source("D:/orthoptera/orthoptera_prediction/src/00_set_environment.R")
} else {
  source("/media/tnauss/myWork/analysis/orthoptera/orthoptera_prediction/src/00_set_environment.R")
}

compute = TRUE
cores = 4

# Predict dataset --------------------------------------------------------------
if(compute){
  
  
  # GLS
  obsv_gpm = readRDS(file = paste0(path_results, "obsv_gls_gpm_traintest.rds"))
  n_var = c(seq(1, 20), length(obsv_gpm@meta$input$PREDICTOR_FINAL))
  
  if(length(showConnections()) == 0){
    cl = parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
  }
  
  obsv_gpm = trainModel(x = obsv_gpm,
                        metric = "Accuracy",
                        n_var = NULL, 
                        mthd = "pls",
                        mode = "ffs",
                        seed_nbr = 11, 
                        cv_nbr = NULL,
                        var_selection = "indv",
                        filepath_tmp = path_temp)
  saveRDS(obsv_gpm, file = paste0(path_results, "obsv_gls_gpm_trainmodel.rds"))
  
  
  # Terra-MODIS
  obsv_modis_mod = readRDS(file = paste0(path_results, "obsv_gpm_mod_traintest.rds"))
  
  n_var = c(seq(1, 20), length(obsv_modis_mod@meta$input$PREDICTOR_FINAL))
  
  if(length(showConnections()) == 0){
    cl = parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
  }
  
  obsv_modis_mod = trainModel(x = obsv_modis_mod,
                               metric = "Accuracy",
                               n_var = NULL, 
                               mthd = "pls",
                               mode = "ffs",
                               seed_nbr = 11, 
                               cv_nbr = NULL,
                               var_selection = "indv",
                               filepath_tmp = path_temp)
  saveRDS(obsv_modis_mod, file = paste0(path_results, "obsv_mod_gpm_trainmodel.rds"))
  
  
  # Terra-MODIS
  if(length(showConnections()) == 0){
    cl = parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
  }
  
  obsv_modis_myd = readRDS(file = paste0(path_results, "obsv_gpm_myd_traintest.rds"))
  
  n_var = c(seq(1, 20), length(obsv_modis_myd@meta$input$PREDICTOR_FINAL))
  obsv_modis_myd = trainModel(x = obsv_modis_myd,
                               metric = "Accuracy",
                               n_var = NULL, 
                               mthd = "pls",
                               mode = "ffs",
                               seed_nbr = 11, 
                               cv_nbr = NULL,
                               var_selection = "indv",
                               filepath_tmp = path_temp)
  saveRDS(obsv_modis_myd, file = paste0(path_results, "obsv_myd_gpm_trainmodel.rds"))
}
