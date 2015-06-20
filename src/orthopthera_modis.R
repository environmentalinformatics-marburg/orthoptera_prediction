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
any(is.na(orthoptera_cplt[, -7]))

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
prevalent_species <- minimumOccurence(x = observations, selector = plotid,
                                      occurence = "yes", 
                                      resample = 100, thv = 20)
# save(prevalent_species, file = "processed/prevalent_species.rda")




# Compile model evaluation dataset ---------------------------------------------
orthoptera_resamples <- resamplingsByVariable(x = orthoptera@data$input, 
                                              selector = plotid, 
                                              resample = 100)
save(orthoptera_resamples, file = "processed/orthoptera_resamples.rda")




# Split dataset into testing and training samples for each individual species --
# load(""processed/orthoptera_mdl.rda")
col_response <- seq(14,15)
orthoptera_trte <- splitMultResp(x = orthoptera@data$input, 
                                 response = col_response,
                                 resamples = orthoptera_resamples)
save(orthoptera_mdl_trte, file = "processed/orthoptera_mdl_trte.rda")




# Evaluate prediction models ---------------------------------------------------
# load("processed/orthoptera_mdl_trte.rda")
response <- 2
independent <- seq(3, 32)


models <- trainModel(orthoptera_mdl_trte, 
                response_column = 2, independent_columns = seq(3, 32),
                response_nbr = c(1,2), model_instc = 1)
# devtools::use_data(models, overwrite = TRUE)
# load(models)
# load(orthoptera_mdl_trte)


plot(varImp(rfe_model$fit,scale=TRUE))

test <- lapply(seq(100), function(x){
  act_test <- predict(models[[x]], orthoptera_mdl_trte[[1]][[x]]$test[, independent])
  act_obs <- orthoptera_mdl_trte[[1]][[x]]$test[, response]
  calcKappa(ftable(data.frame(PREDICT = act_test$pred, OBSERVERD = act_obs)))[1]
})
summary(unlist(test))