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

  obsv_gpm <- readRDS(file = paste0(path_results, "gls_2000_gpm_traintest.rds"))
  
  # n_var <- c(seq(1, 29), seq(40, length(obsv_gpm@meta$input$PREDICTOR_FINAL), 40))
  n_var <- c(seq(1, length(obsv_gpm@meta$input$PREDICTOR_FINAL)))
  
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  
  obsv_gpm <- trainModel(x = obsv_gpm,
                         n_var = n_var, 
                         mthd = "rf",
                         mode = "rfe",
                         response_nbr = 1,
                         resample_nbr = 1,
                         seed_nbr = 11, 
                         cv_nbr = 5,
                         var_selection = "indv", 
                         filepath_tmp = path_temp)
  
  
  saveRDS(obsv_gpm, file = paste0(path_results, "gls_2000_gpm_trainmodel_pls.rds"))
  
} else {
  obsv_gpm <- readRDS(file = paste0(path_results, "gls_2000_gpm_trainmodel.rds"))
}

var_imp <- compVarImp(obsv_gpm@model$rf_rfe, scale = FALSE)

var_imp_scale <- compVarImp(obsv_gpm@model$rf_rfe, scale = TRUE)

var_imp_plot <- plotVarImp(var_imp)

var_imp_heat <- plotVarImpHeatmap(var_imp_scale, xlab = "Species", ylab = "Band")

tstat <- compContTests(obsv_gpm@model$rf_rfe, mean = TRUE)

tstat_mean <- merge(tstat[[1]], obsv_gpm@meta$input$MIN_OCCURENCE,
                    by.x = "Response", by.y="names")

tstat_mean[order(tstat_mean$Kappa_mean, decreasing = TRUE),]

ggplot(data = tstat_mean, aes(x = mo_mean, y = Kappa_mean)) + geom_point() + geom_smooth()
