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

  obsv_modis_mod <- readRDS(file = paste0(path_results, "obsv_modis_mod_traintest.rds"))
  
  obsv_modis_mod_wgs <- trainModel(x = obsv_modis_mod[["wgs"]],
                         n_var = NULL, 
                         mthd = "rf",
                         mode = "rfe",
                         seed_nbr = 11, 
                         cv_nbr = 5,
                         var_selection = "sd", 
                         response_nbr = seq(3), 
                         resample_nbr = seq(2),
                         filepath_tmp = NULL)
  saveRDS(obsv_modis_mod_wgs, file = paste0(path_results, "obsv_modis_mod_wgs_trainModel.rds"))
  
  obsv_modis_mod_arc <- trainModel(x = obsv_modis_mod[["arc"]],
                                   n_var = NULL, 
                                   mthd = "rf",
                                   mode = "rfe",
                                   seed_nbr = 11, 
                                   cv_nbr = 5,
                                   var_selection = "sd", 
                                   response_nbr = seq(3), 
                                   resample_nbr = seq(2),
                                   filepath_tmp = NULL)
  saveRDS(obsv_modis_mod_arc, file = paste0(path_results, "obsv_modis_mod_arc_trainModel.rds"))
  
  
  obsv_modis_myd <- readRDS(file = paste0(path_results, "obsv_modis_myd_traintest.rds"))
  
  obsv_modis_myd_wgs <- trainModel(x = obsv_modis_myd[["wgs"]],
                                   n_var = NULL, 
                                   mthd = "rf",
                                   mode = "rfe",
                                   seed_nbr = 11, 
                                   cv_nbr = 5,
                                   var_selection = "sd", 
                                   response_nbr = seq(3), 
                                   resample_nbr = seq(2),
                                   filepath_tmp = NULL)
  saveRDS(obsv_modis_myd_wgs, file = paste0(path_results, "obsv_modis_myd_wgs_trainModel.rds"))
  
  obsv_modis_myd_arc <- trainModel(x = obsv_modis_myd[["arc"]],
                                   n_var = NULL, 
                                   mthd = "rf",
                                   mode = "rfe",
                                   seed_nbr = 11, 
                                   cv_nbr = 5,
                                   var_selection = "sd", 
                                   response_nbr = seq(3), 
                                   resample_nbr = seq(2),
                                   filepath_tmp = NULL)
  saveRDS(obsv_modis_myd_arc, file = paste0(path_results, "obsv_modis_myd_arc_trainModel.rds"))
  
} else {
  obsv_modis_mod_wgs <- readRDS(file = paste0(path_results, "obsv_modis_mod_wgs_trainModel.rds"))
  obsv_modis_mod_arc <- readRDS(file = paste0(path_results, "obsv_modis_mod_arc_trainModel.rds"))
  obsv_modis_myd_wgs - readRDS(file = paste0(path_results, "obsv_modis_myd_wgs_trainModel.rds"))
  obsv_modis_myd_arc <- readRDS(file = paste0(path_results, "obsv_modis_myd_arc_trainModel.rds"))
}


# var_imp <- compVarImp(models@model, scale = FALSE)
# 
# var_imp_scale <- compVarImp(models@model, scale = TRUE)
# 
# var_imp_plot <- plotVarImp(var_imp)
# 
# var_imp_heat <- plotVarImpHeatmap(var_imp_scale, xlab = "Species", ylab = "Band")
# 
# tstat <- compContTests(models@model, mean = TRUE)
# 
# tstat_mean <- merge(tstat[[1]], obsv_modis[[prj]]@meta$input$MIN_OCCURENCE, 
#                     by.x = "Response", by.y="names")
# 
# tstat_mean[order(tstat_mean$Kappa_mean, decreasing = TRUE),]
# 
# ggplot(data = tstat_mean, aes(x = mo_mean, y = Kappa_mean)) + geom_point() + geom_smooth()
