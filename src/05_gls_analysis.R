# Predict orthoptera based on satellite observations
# Thomas Nauss

if(Sys.info()["sysname"] == "Windows"){
  source("D:/orthoptera/orthoptera_prediction/src/00_set_environment.R")
} else {
  source("/media/tnauss/myWork/analysis/orthoptera/orthoptera_prediction/src/00_set_environment.R")
}

compute <- TRUE

# Analyse model results --------------------------------------------------------
load(paste0(path_results, "gpm_trainModel_model_instances_004.RData"))
i001 <- model_instances
i002 <- model_instances
i003 <- model_instances
i004 <- model_instances


vi <- lapply(i003, function(i){
  vi <- caret::varImp(i$model)
  variables <- rownames(vi)
  data.frame(VARIABLE = variables,
             IMPORTANCE = vi$Overall)
})
vi_species <- do.call("rbind", vi)
vi_count <- vi_species %>% count(VARIABLE)
vi_mean <- vi_species %>% group_by(VARIABLE) %>% summarise(mean = mean(IMPORTANCE))
vi <- merge(vi_count, vi_mean)
vi$RESPONSE <- vi_species$RESPONSE[1]
vi <- vi[order(vi$n, decreasing = TRUE), ,drop = FALSE]

vi_length <- lapply(i003, function(i){
  length(i$model$finalModel$importance)
})
summary(unlist(vi_length))



caret::varImp(i001[[100]]$model)


# var_imp <- compVarImp(models@model$rf_rfe, scale = FALSE)
# 
# var_imp_scale <- compVarImp(models@model$rf_rfe, scale = TRUE)
# 
# var_imp_plot <- plotVarImp(var_imp)
# 
# var_imp_heat <- plotVarImpHeatmap(var_imp_scale, xlab = "Species", ylab = "Band")
# 
# tstat <- compContTests(models@model$rf_rfe, mean = TRUE)
# 
# tstat_mean <- merge(tstat[[1]], obsv_gpm[[prj]]@meta$input$MIN_OCCURENCE, 
#                     by.x = "Response", by.y="names")
# 
# tstat_mean[order(tstat_mean$Kappa_mean, decreasing = TRUE),]
# 
# ggplot(data = tstat_mean, aes(x = mo_mean, y = Kappa_mean)) + geom_point() + geom_smooth()
