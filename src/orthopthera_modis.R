setwd("D:/active/orthoptera/data")

# Libraries --------------------------------------------------------------------
# library(gpm)




# Read and adjust data from S. Schlauss, level 300 -----------------------------
orthoptera <- read.table("original/lvl0300_biodiversity_data.csv", 
                         header = TRUE, sep = ";", dec = ",")

# Replace number of observations and NAs to 1/0
orthoptera[, 14:178][is.na(orthoptera[, 14:178])] <- 0
orthoptera[, 14:178][orthoptera[, 14:178] > 0] <- 1

str(orthoptera)
col_meta <- seq(1, 13)
col_species <- seq(14, 178)
col_modis <- seq(179, 216)

save(orthoptera, file = "processed/orthoptera.rda")




# Select species which occure at least on 20 unique plots on average -----------
# load("processed/orthoptera.rda")  
plotid <- orthoptera$plot
observations <- orthoptera[, col_species]
prevalent_species <- meanOccupancy(plotid, observations, 
                                   resample = 100, thv = 20)
save(prevalent_species, file = "processed/prevalent_species.rda")




# Compile dataset containing complete cases only -------------------------------
# load("processed/orthoptera.rda")
# load("processed/prevalent_species.rda")
orthoptera_cplt <- orthoptera[, c(col_meta, 
                           which(colnames(orthoptera) %in% prevalent_species),
                           col_modis)]
summary(orthoptera_cplt)
orthoptera_cplt <- 
  orthoptera_cplt[, -(which(colnames(orthoptera_cplt) == "greyval_band_11") : 
                        which(colnames(orthoptera_cplt) == "greyval_band_16"))]
any(is.na(orthoptera_cplt[, -7]))
save(orthoptera_cplt, file = "processed/orthoptera_cplt.rda")




# Compile model evaluation dataset ---------------------------------------------
# Since caret's rfe function does not only allow letter characters as factor,
# change 0/1 to no/yes for the species occurence columns.
# load("processed/orthoptera_cplt.rda")
for(i in seq(14, 34)){
  temp <- as.character(orthoptera_cplt[, i])
  temp[temp == "0"] <- "no"
  temp[temp == "1"] <- "yes"
  orthoptera_cplt[, i] <- as.factor(temp)
}

orthoptera_mdl <- conditionalSample(orthoptera_cplt[, -(2:13)], orthoptera_cplt$plot, 
                               resample = 100)
save(orthoptera_cplt, file = "processed/orthoptera_mdl.rda")




# Split dataset into testing and training samples for each individual species --
# load(""processed/orthoptera_mdl.rda")
col_response <- seq(2,22)
orthoptera_mdl_trte <- splitByFrequency(orthoptera_mdl, col_response)
save(orthoptera_mdl_trte, file = "processed/orthoptera_mdl_trte.rda")




# Evaluate prediction models ---------------------------------------------------
# load("processed/orthoptera_mdl_trte.rda")
response <- 2
independent <- seq(3, 32)


models <- lapply(seq(100), function(x){
  act <- orthoptera_mdl_trte[[1]][[x]]
  act_train <- act$training
  act_test <- act$test
  
  # # Convert response variables to factors
  # act_train[, response] <- as.factor(act_train[, response])
  # act_test[, response] <- as.factor(act_test[, response])
  # 
  # resp2 <- act_train[, response]
  # resp <- rep("no", length(resp2))
  # resp[resp2 == "1"] <- "yes"
  # resp <- as.factor(resp)
  
  resp <- act_train[, response]
  indp <- act_train[, independent]
  
  set.seed(10)
  cv_splits <- createFolds(resp, k=5, returnTrain = TRUE)
  
  # thresholds=c(seq(0.0, 0.40, 0.02),seq(0.50,1,0.1))
  # summaryFunction = "fourStats"
  
  rfeCntrl <- rfeControl(functions = caretFuncs,
                         method="cv", index = cv_splits,
                         returnResamp = "all",
                         verbose = FALSE,
                         rerank=FALSE)
  
  trCntr <- trainControl(method="cv", number = 5, repeats = 1, verbose = FALSE)
  
  n_var <- seq(2, ncol(indp), 8)
  
  
  # method = rf_thvs
  # 
  #   ctrl <- trainControl(index=cvSplits,
  #                        method="cv",
  #                        summaryFunction = eval(parse(text=summaryFunction)),
  #                        classProbs = classProbs)
  
  rfe_model <- rfe(indp, resp,
                   metric = "Accuracy", method = "rf", 
                   sizes = n_var,
                   rfeControl = rfeCntrl,
                   trControl = trCntr, verbose = FALSE,
                   tuneGrid = expand.grid(mtry = n_var))
  
  return(rfe_model)
})

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