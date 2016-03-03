setwd("D:/active/orthoptera/data")

# Libraries --------------------------------------------------------------------
# library(gpm)




# Read and adjust data from S. Schlauss, level 300 -----------------------------
orthoptera <- read.table("original/lvl0300_biodiversity_data.csv", 
                         header = TRUE, sep = ";", dec = ",")

# Replace number of observations and NAs to 1/0
orthoptera[, 14:178][!is.na(orthoptera[, 14:178])] <- "yes"
orthoptera[, 14:178][is.na(orthoptera[, 14:178])] <- "no"
for(i in seq(14, 178)){
  orthoptera[, i] <- as.factor(orthoptera[, i])
}

# Compile dataset containing complete cases only
orthoptera <- 
  orthoptera[, -(which(colnames(orthoptera) == "greyval_band_11") : 
                   which(colnames(orthoptera) == "greyval_band_16"))]
any(is.na(orthoptera[, -7]))

col_meta <- seq(1, 13)
col_species <- seq(14, 178)
col_modis <- seq(179, 208)

meta <- createGPMMeta(orthoptera, type = "input",
                      selector = 1, response = seq(14, 178), 
                      independent = seq(179, 208), meta = c(2: 13))
orthoptera <- gpm(orthoptera, meta)
# save(orthoptera, file = "processed/orthoptera.rda")




# Select responses occuring at least across 20 unique selector values on average
# load("processed/orthoptera.rda")
plotid <- orthoptera@data$input$plot
observations <- orthoptera@data$input[, col_species]
min_occurence <- minimumOccurence(x = observations, selector = plotid,
                                  occurence = "yes", 
                                  resample = 100, thv = 20)
prevalent_species <- min_occurence[[1]]
prevalence <- data.frame(min_occurence[[2]])
prevalence$RESPONSE <- rownames(prevalence)
names(prevalence)[1] <- "OCCURENCE"
rownames(prevalence) <- NULL
# save(prevalent_species, file = "processed/prevalent_species.rda")
# save(prevalence, file = "processed/prevalence.rda")


# Compile model evaluation dataset ---------------------------------------------
orthoptera_resamples <- resamplingsByVariable(x = orthoptera@data$input, 
                                              selector = plotid, 
                                              grabs = 1,
                                              resample = 100)
# save(orthoptera_resamples, file = "processed/orthoptera_resamples.rda")




# Split dataset into testing and training samples for each individual species --
# load("processed/orthoptera.rda")
# load("processed/prevalent_species.rda")
# load("processed/orthoptera_resamples.rda")
col_response <- prevalent_species
orthoptera_trte <- splitMultResp(x = orthoptera@data$input, 
                                 response = col_response,
                                 resamples = orthoptera_resamples)
# save(orthoptera_trte, file = "processed/orthoptera_trte.rda")




# Evaluate prediction models ---------------------------------------------------
load("processed/orthoptera.rda")
load("processed/prevalent_species.rda")
load("processed/orthoptera_trte.rda")
response <- prevalent_species
independent <- orthoptera@meta$input$INDEPENDENT

# models <- trainModel(x = orthoptera@data$input, 
#                      response = response, independent = independent,
#                      resamples = orthoptera_trte, n_var = seq(1,30,5),
#                      mthd = "rf")
# save(models, file = "processed/models_rf-2015-11-26.rda")


models <- trainModel(x = orthoptera@data$input, 
                     response = response, independent = independent,
                     resamples = orthoptera_trte, n_var = seq(1,30,2),
                     response_nbr = seq(5),
                     mthd = "avNNet")


models <- trainModel(x = orthoptera@data$input, 
                     response = response, independent = independent,
                     resamples = orthoptera_trte, n_var = seq(1,30,2),
                     mthd = "avNNet")
# save(models, file = "processed/models_avnnet.rda")
# load("processed/models_rf.rda")

var_imp <- compVarImp(models)

var_imp_plot <- plotVarImp(var_imp)

var_imp_heat <- plotVarImpHeatmap(var_imp, xlab = "Species", ylab = "Band")

# png("var_imp_heat.png")
# var_imp_heat
# dev.off()

tests <- compContTests(models)

tstat_mean <- lapply(tests, function(x){
  data.frame(RESPONSE = x$RESPONSE[1], 
             KAPPA_MEAN = mean(x$Kappa, na.rm = TRUE),
             POD_MEAN = mean(x$POD, na.rm = TRUE),
             FAR_MEAN = mean(x$FAR, na.rm = TRUE), 
             POFD_MEAN = mean(x$POFD, na.rm = TRUE),
             ACCURACY_MEAN = mean(x$ACCURACY, na.rm = TRUE),
             SR_MEAN = mean(x$SR, na.rm = TRUE),
             TS_MEAN = mean(x$TS, na.rm = TRUE),
             ETS_MEAN = mean(x$ETS, na.rm = TRUE),
             HK_MEAN = mean(x$HK, na.rm = TRUE))
})
tstat_mean <- do.call("rbind", tstat_mean)
tstat_mean <- merge(tstat_mean, prevalence, by = "RESPONSE")
tstat_mean[order(tstat_mean$KAPPA_MEAN, decreasing = TRUE),]
# save(tstat_mean, file = "processed/tstat_mean_rf.rda")
library(corrplot)
corrplot(cor(tstat_mean[, -1]))

