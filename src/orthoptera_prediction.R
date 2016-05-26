# Predcit orthoptera occurence
library(gpm)
library(ggplot2)

if(Sys.info()["sysname"] == "Windows"){
  filepath_base <- "E:/analysis/orthoptera/"
} else {
  filepath_base <- "/media/tnauss/myWork/analysis/orthoptera/"
}

filepath_obsv <- paste0(filepath_base, "data/orthoptera/")
filepath_results <- paste0(filepath_base, "data/gpm_results/")


# Prepare observations ---------------------------------------------------------
# Read observation dataset
obsv <- read.table(paste0(filepath_obsv, "orthoptera_sat.csv"),
                   sep = ",", dec = ".", header = TRUE)

# Replace number of observations and NAs to 1/0
colnames(obsv)

obsv[, 14:178][!is.na(obsv[, 14:178])] <- "yes"
obsv[, 14:178][is.na(obsv[, 14:178])] <- "no"
for(i in seq(14, 178)){
  obsv[, i] <- as.factor(obsv[, i])
}

# Compile dataset containing complete cases only
obsv <- obsv[, -which(colnames(obsv) == "gls2000_hist_glcm_correlation")]
any(is.na(obsv[, -7]))

col_meta <- seq(1, 13)
col_species <- seq(14, 178)
col_sat <- seq(179, ncol(obsv))

meta <- createGPMMeta(obsv, type = "input",
                      selector = 1, response = col_species, 
                      independent = col_sat, meta = col_meta)
obsv <- gpm(obsv, meta)
# save(obsv, file = paste0(filepath_results, "obsv.rda"))


# Sampling ---------------------------------------------------------------------
# Select responses occuring at least across 20 unique selector values on average
# load(paste0(filepath_results, "obsv.rda"))
plotid <- obsv@data$input$plot
observations <- obsv@data$input[, col_species]
min_occurence <- minimumOccurence(x = observations, selector = plotid,
                                  occurence = "yes", 
                                  resample = 100, thv = 20)
prevalent_species <- min_occurence[[1]]
prevalence <- data.frame(min_occurence[[2]])
prevalence$RESPONSE <- rownames(prevalence)
names(prevalence)[1] <- "OCCURENCE"
rownames(prevalence) <- NULL
# save(prevalent_species, file = paste0(filepath_results, "prevalent_species.rda"))
# save(prevalence, file = paste0(filepath_results, "prevalence.rda"))


# Compile model evaluation dataset ---------------------------------------------
orthoptera_resamples <- resamplingsByVariable(x = obsv@data$input, 
                                              selector = plotid, 
                                              grabs = 1,
                                              resample = 100)
# save(orthoptera_resamples, file = paste0(filepath_results, "orthoptera_resamples.rda"))


# Split dataset into testing and training samples for each individual species --
# load(paste0(filepath_results, "obsv.rda"))
# load("processed/prevalent_species.rda")
# load("processed/orthoptera_resamples.rda")
col_response <- prevalent_species
orthoptera_trte <- splitMultResp(x = obsv@data$input, 
                                 response = col_response,
                                 resamples = orthoptera_resamples)
# save(orthoptera_trte, file = paste0(filepath_results, "orthoptera_trte.rda"))




# Evaluate prediction models ---------------------------------------------------
# load(paste0(filepath_results, "orthoptera.rda"))
# load(paste0(filepath_results, "prevalent_species.rda"))
# load(paste0(filepath_results, "orthoptera_trte.rda"))
# load(paste0(filepath_results, "obsv.rda"))
response <- prevalent_species
independent <- obsv@meta$input$INDEPENDENT

models <- trainModel(x = obsv, 
                     response = response, independent = independent,
                     resamples = orthoptera_trte,  mode = "ffs",
                     n_var = seq(1,30,5), resample_nbr = 2, response_nbr = seq(2),
                     mthd = "rf", seed_nbr = 11, cv_nbr = 2)

models <- trainModel(x = obsv, 
                     response = response, independent = independent,
                     resamples = orthoptera_trte, n_var = seq(1,30,5),
                     mthd = "rf", seed_nbr = 11, cv_nbr = 2)
# save(models, file = "processed/models_rf-2016-03-04.rda")


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
# save(models, file = "processed/models_avnnet.rda")
# load("processed/models_rf-2015-11-26.rda")
# load("processed/models_rf.rda")
load("processed/models_rf-2016-03-04.rda")

var_imp <- compVarImp(models, scale = FALSE)

var_imp_scale <- compVarImp(models, scale = TRUE)

var_imp_plot <- plotVarImp(var_imp)

var_imp_heat <- plotVarImpHeatmap(var_imp_scale, xlab = "Species", ylab = "Band")

tstat <- compContTests(models, mean = TRUE)

tstat_mean <- merge(tstat[[1]], prevalence, by.x = "Response", by.y = "RESPONSE")


tstat_mean[order(tstat_mean$Kappa_mean, decreasing = TRUE),]

ggplot(data = tstat_mean, aes(x = OCCURENCE, y = Kappa_mean)) + geom_point() + geom_smooth()

# save(tstat, tstat_mean, file = "processed/tstat_mean_rf.rda")



