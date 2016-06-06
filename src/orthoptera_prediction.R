# Predcit orthoptera occurence
library(gpm)
library(ggplot2)

if(Sys.info()["sysname"] == "Windows"){
  filepath_base <- "E:/analysis/orthoptera/data/"
} else {
  filepath_base <- "/media/tnauss/myWork/analysis/orthoptera/data/"
}

filepath_obsv <- paste0(filepath_base, "orthoptera/")
filepath_results <- paste0(filepath_base, "rdata/")



# Prepare observations ---------------------------------------------------------
# Read observation dataset
load(paste0(filepath_results, "preprocessing_data_obsv.RData"))

# Replace number of observations and NAs to 1/0
colnames(obsv)

col_meta <- seq(1, 14)
col_species <- seq(15, 179)
col_sat <- seq(180, ncol(obsv))

obsv[, col_species][!is.na(obsv[, col_species])] <- "yes"
obsv[, col_species][is.na(obsv[, col_species])] <- "no"
for(i in col_species){
  obsv[, i] <- as.factor(obsv[, i])
}

obsv[, grep("MOD", colnames(obsv))] <- obsv[, grep("MOD", colnames(obsv))]/10000.0
obsv[, grep("MYD", colnames(obsv))] <- obsv[, grep("MYD", colnames(obsv))]/10000.0

# Compile dataset containing complete cases only
summary(obsv)

meta <- createGPMMeta(obsv, type = "input",
                      selector = 1, response = col_species, 
                      independent = col_sat, meta = col_meta)
obsv_gpm <- gpm(obsv, meta, scale = TRUE)
# summary(obsv_gpm@data$input)
# save(obsv_gpm, file = paste0(filepath_results, "orthoptera_prediction_obsv_gpm.RData"))


# Sampling ---------------------------------------------------------------------
# Select responses occuring at least across 20 unique selector values on average
# load(paste0(filepath_results, "orthoptera_prediction_obsv_gpm.RData"))
plotid <- obsv_gpm@data$input$plot
observations <- obsv_gpm@data$input[, col_species]
min_occurence <- minimumOccurence(x = observations, selector = plotid,
                                  occurence = "yes", 
                                  resample = 100, thv = 20)
prevalence <- data.frame(min_occurence[[2]])
prevalence$RESPONSE <- rownames(prevalence)
names(prevalence)[1] <- "OCCURENCE"
rownames(prevalence) <- NULL
# save(prevalence, file = paste0(filepath_results, "orthoptera_prediction_prevalence.RData"))


# Compile model evaluation dataset ---------------------------------------------
orthoptera_resamples <- resamplingsByVariable(x = obsv_gpm@data$input, 
                                              selector = plotid, 
                                              grabs = 1,
                                              resample = 100)
# save(orthoptera_resamples, file = paste0(filepath_results, "orthoptera_prediction_orthoptera_resamples.RData"))


# Split dataset into testing and training samples for each individual species --
# load(file = paste0(filepath_results, "orthoptera_prediction_obsv_gpm.RData"))
# load(paste0(filepath_results, "orthoptera_prediction_orthoptera_resamples.RData"))
# load(paste0(filepath_results, "orthoptera_prediction_prevalence.RData"))
orthoptera_trte <- splitMultResp(x = obsv_gpm@data$input, 
                                 response = prevalence$RESPONSE,
                                 resamples = orthoptera_resamples)
# save(orthoptera_trte, file = paste0(filepath_results, "orthoptera_prediction_orthoptera_trte.RData"))


# Evaluate prediction models ---------------------------------------------------
# load(file = paste0(filepath_results, "orthoptera_prediction_obsv_gpm.RData"))
# load(file = paste0(filepath_results, "orthoptera_prediction_prevalence.RData"))
# load(paste0(filepath_results, "orthoptera_prediction_orthoptera_resamples.RData"))
# load(paste0(filepath_results, "orthoptera_prediction_orthoptera_trte.RData"))
# Check for NA and remove those columns
independent <- obsv_gpm@meta$input$INDEPENDENT
# independent <- c(independent, "asl")
independent <- independent[sapply(independent, function(x){!any(is.na(obsv_gpm@data$input[,x]))})]

n_vars <- c(seq(1,15,1), seq(16, length(independent), 3))
models <- trainModel(x = obsv_gpm, 
                     response = prevalence$RESPONSE, independent = independent,
                     resamples = orthoptera_trte,  mode = "ffs",
                     n_var = n_vars, 
                     mthd = "rf", seed_nbr = 11, cv_nbr = 2)
save(models, file = paste0(filepath_results, "orthoptera_prediction_models_rf_2016-05-29_ffs.RData"))

# n_vars <- c(seq(1,15,1), seq(16, length(independent), 3))
# models <- trainModel(x = obsv_gpm, 
#                      response = prevalence$RESPONSE, independent = independent,
#                      resamples = orthoptera_trte, n_var = n_vars,
#                      mthd = "rf", seed_nbr = 11, cv_nbr = 2,
#                      filepath_tmp = filepath_results)
# save(models, file = paste0(filepath_results, "orthoptera_prediction_models_rf_2016-05-29_rfe_als.RData"))

# model_files <- list.files(filepath_results, pattern = glob2rx("gpm_trainModel_model_instances_*"),
#                           full.names = TRUE)
# models <- lapply(model_files, function(x){
#   load(x)
#   return(model_instances)
# })


# models <- trainModel(x = orthoptera@data$input, 
#                      response = response, independent = independent,
#                      resamples = orthoptera_trte, n_var = seq(1,30,2),
#                      response_nbr = seq(5),
#                      mthd = "avNNet")
# 
# 
# models <- trainModel(x = orthoptera@data$input, 
#                      response = response, independent = independent,
#                      resamples = orthoptera_trte, n_var = seq(1,30,2),
#                      mthd = "avNNet")
# save(models, file = "processed/models_avnnet.RData")
# load("processed/models_rf-2015-11-26.RData")
# load("processed/models_rf.RData")

var_imp <- compVarImp(models, scale = FALSE)

var_imp_scale <- compVarImp(models, scale = TRUE)

var_imp_plot <- plotVarImp(var_imp)

var_imp_heat <- plotVarImpHeatmap(var_imp_scale, xlab = "Species", ylab = "Band")

tstat <- compContTests(models, mean = TRUE)

tstat_mean <- merge(tstat[[1]], prevalence, by.x = "Response", by.y = "RESPONSE")

tstat_mean[order(tstat_mean$Kappa_mean, decreasing = TRUE),]

ggplot(data = tstat_mean, aes(x = OCCURENCE, y = Kappa_mean)) + geom_point() + geom_smooth()

# save(tstat, tstat_mean, file = "processed/tstat_mean_rf.RData")


# 
# 
# response_nbr <- seq(length(response))
# resample_nbr <- seq(length(orthoptera_trte[[1]]))
# response_instances <- lapply(response_nbr, function(i){
#   model_instances <- lapply(resample_nbr, function(j){
#     print(paste0("Computing resample instance ", j, 
#                  " of response instance ", i, "..."))
#     act_resample <- orthoptera_trte[[i]][[j]]
#     resp <- obsv_gpm@data$input[act_resample$training$SAMPLES, 
#                          act_resample$training$RESPONSE]
#     test <- grep("yes", resp)
#     return(test)
#     })
#   })    
# response_instances
# 
# do.call("rbind", response_instances)
# 
# t <- lapply(seq(5), function(y){
#   print(y)
#   res <- try(x>5)
# #    if(inherits(res, "try-error"))
# #    {
# #      res <- res
# #    }
# })
