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
  
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  
  obsv_gpm_wgs <- trainModel(x = obsv_gpm,
                         n_var = NULL, 
                         mthd = "rf",
                         mode = "rfe",
                         seed_nbr = 11, 
                         cv_nbr = 5,
                         var_selection = "sd", 
                         filepath_tmp = path_temp)
  saveRDS(obsv_gpm_wgs, file = paste0(path_results, "obsv_gls_gpm_wgs_trainModel.rds"))
  
  obsv_gpm_arc <- trainModel(x = obsv_gpm[["arc"]],
                             n_var = NULL, 
                             mthd = "rf",
                             mode = "rfe",
                             seed_nbr = 11, 
                             cv_nbr = 5,
                             var_selection = "sd", 
                             response_nbr = seq(2), 
                             resample_nbr = seq(2),
                             filepath_tmp = path_temp)
  saveRDS(obsv_gpm_arc, file = paste0(path_results, "obsv_gls_gpm_arc_trainModel.rds"))

} else {
  obsv_gpm_wgs <- readRDS(file = paste0(path_results, "obsv_gls_gpm_wgs_trainModel.rds"))
  obsv_gpm_arc <- readRDS(file = paste0(path_results, "obsv_gls_gpm_arc_trainModel.rds"))
}

# var_imp <- compVarImp(models@model$rf_rfe, scale = FALSE)
# 
# var_imp_scale <- compVarImp(models@model$rf_rfe, scale = TRUE)
# 
# var_imp_plot <- plotVarImp(var_imp)
# 
# var_imp_heat <- plotVarImpHeatmap(var_imp_scale, xlab = "Species", ylab = "Band")
# 
# tstat <- compContTests(models@model$rf_rfe, mean = TRUE)
# 
# tstat_mean <- merge(tstat[[1]], obsv_gpm[[prj]]@meta$input$MIN_OCCURENCE, 
#                     by.x = "Response", by.y="names")
# 
# tstat_mean[order(tstat_mean$Kappa_mean, decreasing = TRUE),]
# 
# ggplot(data = tstat_mean, aes(x = mo_mean, y = Kappa_mean)) + geom_point() + geom_smooth()
