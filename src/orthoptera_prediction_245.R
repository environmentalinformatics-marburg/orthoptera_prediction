# Predcit orthoptera occurence
library(gpm)

if(Sys.info()["sysname"] == "Windows"){
  filepath_base <- "E:/analysis/orthoptera/data/"
} else {
  filepath_base <- "/home/dogbert/orthoptera/"
}

filepath_results <- paste0(filepath_base, "rdata/")


# Evaluate prediction models ---------------------------------------------------
load(file = paste0(filepath_results, "orthoptera_prediction_obsv_gpm.RData"))
load(file = paste0(filepath_results, "orthoptera_prediction_prevalence.RData"))
load(paste0(filepath_results, "orthoptera_prediction_orthoptera_resamples.RData"))
load(paste0(filepath_results, "orthoptera_prediction_orthoptera_trte.RData"))

# Check for NA and remove those columns
independent <- obsv_gpm@meta$input$INDEPENDENT
independent <- c(independent, "asl")
independent <- independent[sapply(independent, function(x){!any(is.na(obsv_gpm@data$input[,x]))})]

n_vars <- c(seq(length(independent)))
models <- trainModel(x = obsv_gpm,
                     response = prevalence$RESPONSE, independent = independent,
                     resamples = orthoptera_trte, n_var = n_vars,
                     mthd = "rf", seed_nbr = 11, cv_nbr = 5,
                     var_selection = "sd",
                     filepath_tmp = filepath_results)
save(models, file = paste0(filepath_results, "orthoptera_prediction_models_rf_2016-06-06_rfe_sd_245.RData"))
load(paste0(filepath_results, "gpm_trainModel_model_instances_001.RData"))
