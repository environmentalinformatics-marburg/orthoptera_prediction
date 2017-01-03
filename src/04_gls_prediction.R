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

  obsv_gpm <- readRDS(file = paste0(path_results, "gls_obsv_gpm_traintest.rds"))
  
  cprj <- c("wgs", "arc")
  obsv_gpm <- lapply(cprj, function(prj){
    models <- trainModel(x = obsv_gpm[[prj]],
                         n_var = NULL, 
                         mthd = "rf",
                         mode = "rfe",
                         seed_nbr = 11, 
                         cv_nbr = 5,
                         var_selection = "sd", 
                         response_nbr = 2, 
                         resample_nbr = seq(9),
                         filepath_tmp = NULL)
    })
  names(obsv_gpm) <- cprj
  saveRDS(obsv_gpm, file = paste0(path_results, "gls_obsv_gpm_trainModel.rds"))
} else {
  obsv_gpm <- readRDS(file = paste0(path_results, "gls_obsv_gpm_trainModel.rds"))
}


var_imp <- compVarImp(models@model, scale = FALSE)

var_imp_scale <- compVarImp(models@model, scale = TRUE)

var_imp_plot <- plotVarImp(var_imp)

var_imp_heat <- plotVarImpHeatmap(var_imp_scale, xlab = "Species", ylab = "Band")

tstat <- compContTests(models@model, mean = TRUE)

tstat_mean <- merge(tstat[[1]], obsv_gpm[[prj]]@meta$input$MIN_OCCURENCE, 
                    by.x = "Response", by.y="names")

tstat_mean[order(tstat_mean$Kappa_mean, decreasing = TRUE),]

ggplot(data = tstat_mean, aes(x = mo_mean, y = Kappa_mean)) + geom_point() + geom_smooth()
