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

  obsv_mod <- readRDS(file = paste0(path_results, "myd_gpm_traintest.rds"))
  
  n_var <- seq(length(obsv_mod@meta$input$PREDICTOR_FINAL))
  
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  
  obsv_mod <- trainModel(x = obsv_mod,
                         n_var = NULL, 
                         mthd = "rf",
                         mode = "rfe",
                         seed_nbr = 11, 
                         cv_nbr = 2,
                         response_nbr = c(2),
                         var_selection = "indv", 
                         metric = "Kappa",
                         tune_length = 5,
                         filepath_tmp = path_temp)
  len = length(obsv_mod@model$rf_rfe[[1]])
  t = lapply(seq(len), function(i){
    data.frame(r = obsv_mod@model$rf_rfe[[1]][[i]]$testing$RESPONSE,
    p = obsv_mod@model$rf_rfe[[1]][[i]]$testing$PREDICTED)
  })
  t = do.call("rbind", t)
  obsv_mod@model$rf_rfe[[1]][[5]]$model
  
  obsv_myd = readRDS(file = paste0(path_results, "bu_2017-12-31/mod_gpm_trainmodel.rds"))
  len = length(obsv_myd@model$rf_rfe[[1]])
  t = lapply(seq(len), function(i){
    data.frame(r = obsv_myd@model$rf_rfe[[3]][[i]]$testing$RESPONSE,
               p = obsv_myd@model$rf_rfe[[3]][[i]]$testing$PREDICTED)
  })
  t = do.call("rbind", t)
  obsv_myd@model$rf_rfe[[3]][[i]]$model
  t$p[t$r]
  
  saveRDS(obsv_mod, file = paste0(path_results, "mod_gpm_trainmodel.rds"))
  

  obsv_myd <- readRDS(file = paste0(path_results, "myd_gpm_traintest.rds"))
  
  n_var <- c(seq(1, 29), seq(40, length(obsv_myd@meta$input$PREDICTOR_FINAL), 40))
  
  obsv_myd <- trainModel(x = obsv_myd,
                                   n_var = NULL, 
                                   mthd = "rf",
                                   mode = "rfe",
                                   seed_nbr = 11, 
                                   cv_nbr = 2,
                                   var_selection = "indv", 
                                   filepath_tmp = NULL)
  saveRDS(obsv_myd, file = paste0(path_results, "myd_gpm_trainmodel.rds"))
  

} else {
  obsv_mod_wgs <- readRDS(file = paste0(path_results, "obsv_mod_wgs_trainModel.rds"))
  obsv_mod_arc <- readRDS(file = paste0(path_results, "obsv_mod_arc_trainModel.rds"))
  obsv_myd_wgs - readRDS(file = paste0(path_results, "obsv_myd_wgs_trainModel.rds"))
  obsv_myd_arc <- readRDS(file = paste0(path_results, "obsv_myd_arc_trainModel.rds"))
}


var_imp <- compVarImp(obsv_mod@model$rf_rfe, scale = FALSE)

var_imp_scale <- compVarImp(obsv_mod@model$rf_rfe, scale = TRUE)

var_imp_plot <- plotVarImp(var_imp)

var_imp_heat <- plotVarImpHeatmap(var_imp_scale, xlab = "Species", ylab = "Band")

tstat_mod <- compContTests(obsv_mod@model$rf_rfe, mean = TRUE)
tstat_myd <- compContTests(obsv_myd@model$rf_rfe, mean = TRUE)

tstat_mean <- merge(tstat[[1]], obsv_gpm@meta$input$MIN_OCCURENCE,
                    by.x = "Response", by.y="names")

tstat_mean[order(tstat_mean$Kappa_mean, decreasing = TRUE),]

ggplot(data = tstat_mean, aes(x = mo_mean, y = Kappa_mean)) + geom_point() + geom_smooth()
