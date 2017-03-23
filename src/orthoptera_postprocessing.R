# Libraries --------------------------------------------------------------------
library(gpm)
# library(ggplot2)
# library(ggmap)
# library(mapview)
library(rgdal)
library(raster)
library(sp)
library(reshape2)
# library(quantreg)
# library(latticeExtra)


# Set path ---------------------------------------------------------------------
if(Sys.info()["sysname"] == "Windows"){
  filepath_base <- "F:/analysis/orthoptera/data/"
} else {
  filepath_base <- "/media/tnauss/myWork/analysis/orthoptera/data/"
}

filepath_results <- paste0(filepath_base, "rdata/")
filepath_obsv <- paste0(filepath_base, "orthoptera/")
filepath_traits <- paste0(filepath_base, "traits/")
filepath_figures <- paste0(filepath_base, "figures/")


# Read data --------------------------------------------------------------------
load(paste0(filepath_results, "orthoptera_prediction_obsv_gpm.RData"))
load(paste0(filepath_results, "orthoptera_prediction_prevalence.RData"))
load(paste0(filepath_results, "orthoptera_prediction_orthoptera_resamples.RData"))
load(paste0(filepath_results, "orthoptera_prediction_orthoptera_trte.RData"))

dem <- raster(paste0(filepath_obsv, "DEM_UTM37S_WGS84_30m_Hemp.tif"))
traits <- read.table(paste0(filepath_traits, "Orthoptera_10_11_2015.csv"),
                     header = TRUE, sep = ";", dec = ".")





load(paste0(filepath_results, "orthoptera_prediction_models_rf_2016-06-06_rfe_sd_mspc.RData"))
mspc <- models

load(paste0(filepath_results, "orthoptera_prediction_models_rf_2016-06-06_rfe_sd_lspt.RData"))
lspt <- models

load(paste0(filepath_results, "orthoptera_prediction_models_rf_2016-06-06_rfe_sd_mspc_lspt.RData"))
mspc_lspt <- models

load(paste0(filepath_results, "orthoptera_prediction_models_rf_2016-06-06_rfe_sd_mspc_lspt_asl.RData"))
mspc_lspt_asl <- models

rm(models)

models_01 <- list(mspc = mspc, lspt = lspt, mspc_lspt = mspc_lspt, 
               mspc_lspt_asl = mspc_lspt_asl)
# save(models_01, file = paste0(filepath_results, "orthoptera_postprocessing_models_01.RData"))
# rm(mspc, lspt, mspc_lspt, mspc_lspt_asl, models_01)
load(file = paste0(filepath_results, "orthoptera_postprocessing_models_01.RData"))

load(paste0(filepath_results, "orthoptera_prediction_models_rf_2016-06-07_rfe_sd_lspc.RData"))
lspc <- models

load(paste0(filepath_results, "orthoptera_prediction_models_rf_2016-06-06_rfe_sd_lspc_lspt.RData"))
lspc_lspt <- models

load(paste0(filepath_results, "orthoptera_prediction_models_rf_2016-06-07_rfe_sd_lspc_lspt_als.RData"))
lspc_lspt_als <- models

load(paste0(filepath_results, "orthoptera_prediction_models_rf_2016-06-06_rfe_sd_lspc_lspt_mspc.RData"))
lspc_lspt_mspc <- models

load(paste0(filepath_results, "orthoptera_prediction_models_rf_2016-06-06_rfe_sd_lspc_lspt_mspc_als.RData"))
lspc_lspt_mspc_als <- models

rm(models)

models_02 <- list(lspc = lspc, lspc_lspt = lspc_lspt, lspc_lspt_als = lspc_lspt_als, 
               lspc_lspt_mspc = lspc_lspt_mspc, lspc_lspt_mspc_als = lspc_lspt_mspc_als)
# save(models_02, file = paste0(filepath_results, "orthoptera_postprocessing_models_02.RData"))
# rm(lspc, lspc_lspt, lspc_lspt_als, lspc_lspt_mspc, lspc_lspt_mspc_als)
load(file = paste0(filepath_results, "orthoptera_postprocessing_models_02.RData"))


# Compute test statistics ------------------------------------------------------
tstat_mspc <- compContTests(models_01$mspc, mean = TRUE)
tstat_mspc_mean <- merge(tstat_mspc[[1]], prevalence, by.x = "Response", by.y = "RESPONSE")

tstat_lspt <- compContTests(models_01$lspt, mean = TRUE)
tstat_lspt_mean <- merge(tstat_lspt[[1]], prevalence, by.x = "Response", by.y = "RESPONSE")

tstat_mspc_lspt <- compContTests(models_01$mspc_lspt, mean = TRUE)
tstat_mspc_lspt_mean <- merge(tstat_mspc_lspt[[1]], prevalence, by.x = "Response", by.y = "RESPONSE")

tstat_mspc_lspt_asl <- compContTests(models_01$mspc_lspt_asl, mean = TRUE)
tstat_mspc_lspt_asl_mean <- merge(tstat_mspc_lspt_asl[[1]], prevalence, by.x = "Response", by.y = "RESPONSE")

tstat_merged_01 <- merge(tstat_mspc_mean, tstat_lspt_mean, by = "Response")
colnames(tstat_merged_01)[2:ncol(tstat_merged_01)] <- c(paste0(colnames(tstat_mspc_mean)[2:ncol(tstat_mspc_mean)], "_mspc"), 
                                                  paste0(colnames(tstat_mspc_mean)[2:ncol(tstat_mspc_mean)], "_lspt"))

tstat_merged_01 <- merge(tstat_merged_01, tstat_mspc_lspt_mean, by = "Response")
colnames(tstat_merged_01)[(ncol(tstat_merged_01)-ncol(tstat_mspc_mean)+2):ncol(tstat_merged_01)] <- 
  paste0(colnames(tstat_mspc_mean)[2:ncol(tstat_mspc_mean)], "_mspc_lspt")

tstat_merged_01 <- merge(tstat_merged_01, tstat_mspc_lspt_asl_mean, by = "Response")
colnames(tstat_merged_01)[(ncol(tstat_merged_01)-ncol(tstat_mspc_mean)+2):ncol(tstat_merged_01)] <- 
  paste0(colnames(tstat_mspc_mean)[2:ncol(tstat_mspc_mean)], "_mspc_lspt_asl")

tstat_merged_01 <- melt(tstat_merged_01, id.vars = "Response")

tstat_01 <- list(tstat_mspc = tstat_mspc, 
              tstat_lspt = tstat_lspt, 
              tstat_mspc_lspt = tstat_mspc_lspt, 
              tstat_mspc_lspt_asl = tstat_mspc_lspt_asl)
tstat_mean_01 <- list(tstat_mspc_mean = tstat_mspc_mean, 
                   tstat_lspt_mean = tstat_lspt_mean, 
                   tstat_mspc_lspt_mean = tstat_mspc_lspt_mean, 
                   tstat_mspc_lspt_asl_mean = tstat_mspc_lspt_asl_mean)
save(tstat_01, tstat_mean_01, tstat_merged_01,
     file = paste0(filepath_results, "orthoptera_postprocessing_tstat_01.RData"))

# ggplot(data = tstat_merged[grepl("POFD", tstat_merged$variable),], 
#        aes(x = Response, y = value, color = variable, group = variable)) + 
#   geom_line() + 
#   coord_flip()



tstat_lspc <- compContTests(models_02$lspc, mean = TRUE)
tstat_lspc_mean <- merge(tstat_lspc[[1]], prevalence, by.x = "Response", by.y = "RESPONSE")

tstat_lspc_lspt <- compContTests(models_02$lspc_lspt, mean = TRUE)
tstat_lspc_lspt_mean <- merge(tstat_lspc_lspt[[1]], prevalence, by.x = "Response", by.y = "RESPONSE")

tstat_lspc_lspt_als <- compContTests(models_02$lspc_lspt_als, mean = TRUE)
tstat_lspc_lspt_als_mean <- merge(tstat_lspc_lspt_als[[1]], prevalence, by.x = "Response", by.y = "RESPONSE")

tstat_lspc_lspt_mspc <- compContTests(models_02$lspc_lspt_mspc, mean = TRUE)
tstat_lspc_lspt_mspc_mean <- merge(tstat_lspc_lspt_mspc[[1]], prevalence, by.x = "Response", by.y = "RESPONSE")

tstat_lspc_lspt_mspc_als <- compContTests(models_02$lspc_lspt_mspc_als, mean = TRUE)
tstat_lspc_lspt_mspc_als_mean <- merge(tstat_lspc_lspt_mspc_als[[1]], prevalence, by.x = "Response", by.y = "RESPONSE")


tstat_merged_02 <- merge(tstat_lspc_mean, tstat_lspc_lspt_mean, by = "Response")
colnames(tstat_merged_02)[2:ncol(tstat_merged_02)] <- c(paste0(colnames(tstat_lspc_mean)[2:ncol(tstat_lspc_mean)], "_lspc"), 
                                                  paste0(colnames(tstat_lspc_lspt_mean)[2:ncol(tstat_lspc_lspt_mean)], "_lspc_lspt"))

tstat_merged_02 <- merge(tstat_merged_02, tstat_lspc_lspt_als_mean, by = "Response")
colnames(tstat_merged_02)[(ncol(tstat_merged_02)-ncol(tstat_lspc_lspt_mean)+2):ncol(tstat_merged_02)] <- 
  paste0(colnames(tstat_lspc_lspt_mean)[2:ncol(tstat_lspc_lspt_mean)], "_lspc_lspt_als")

tstat_merged_02 <- merge(tstat_merged_02, tstat_lspc_lspt_mspc_mean, by = "Response")
colnames(tstat_merged_02)[(ncol(tstat_merged_02)-ncol(tstat_lspc_lspt_mean)+2):ncol(tstat_merged_02)] <- 
  paste0(colnames(tstat_lspc_lspt_mean)[2:ncol(tstat_lspc_lspt_mean)], "_lspc_lspt_mspc")

tstat_merged_02 <- merge(tstat_merged_02, tstat_lspc_lspt_mspc_als_mean, by = "Response")
colnames(tstat_merged_02)[(ncol(tstat_merged_02)-ncol(tstat_lspc_lspt_mean)+2):ncol(tstat_merged_02)] <- 
  paste0(colnames(tstat_lspc_lspt_mean)[2:ncol(tstat_lspc_lspt_mean)], "_lspc_lspt_mspc_als")

tstat_merged_02 <- melt(tstat_merged_02, id.vars = "Response")

tstat_02 <- list(tstat_lspc = tstat_lspc, 
              tstat_lspc_lspt = tstat_lspc_lspt, 
              tstat_lspc_lspt_als = tstat_lspc_lspt_als, 
              tstat_lspc_lspt_mspc = tstat_lspc_lspt_mspc,
              tstat_lspc_lspt_mspc_als = tstat_lspc_lspt_mspc_als)

tstat_mean_02 <- list(tstat_lspc_mean = tstat_lspc_mean, 
                   tstat_lspc_lspt_mean = tstat_lspc_lspt_mean, 
                   tstat_lspc_lspt_als_mean = tstat_lspc_lspt_als_mean, 
                   tstat_lspc_lspt_mspc_mean = tstat_lspc_lspt_mspc_mean,
                   tstat_lspc_lspt_mspc_als_mean = tstat_lspc_lspt_mspc_als_mean)

save(tstat_02, tstat_mean_02, tstat_merged_02,
     file = paste0(filepath_results, "orthoptera_postprocessing_tstat_02.RData"))





# Reduce model dataset to species for which traits are available----------------
reduce2traits_models <- function(models, trait_names){
  del_nbr <- lapply(seq(length(models)), function(x){
    species <- models[[x]][[1]]$response
    if(species %in% trait_names){
      del_nbr <- x
    } else {
      del_nbr <- NA
    }
    return(del_nbr)
  })
  del_nbr <- unlist(del_nbr)
  del_nbr <- del_nbr[!is.na(del_nbr)]
  models_subset <- models[del_nbr]
}

reduce2traits_tstat <- function(stat, trait_names){
  tstat_subset_01 <- stat[[1]][stat[[1]]$Response %in% trait_names, ]
  del_nbr <- lapply(seq(length(stat[[2]])), function(x){
    species <- unique(stat[[2]][[x]]$Response)
    if(species %in% trait_names){
      del_nbr <- x
    } else {
      del_nbr <- NA
    }
    return(del_nbr)
  })
  del_nbr <- unlist(del_nbr)
  del_nbr <- del_nbr[!is.na(del_nbr)]
  tstat_subset_02 <- stat[[2]][del_nbr]
  list(tstat_subset_01, tstat_subset_02)
}

# Species for which trait data is available
trait_names <- gsub("\\ ", ".", traits$Name)

unique(trait_names[!(trait_names %in% prevalence$RESPONSE)])
trait_names[trait_names == "Cyrtacanthacris.tatarica.tatarica"] <- "Cyrtacanthacris.tatarica"
traits$Name <- trait_names

# Models
models_sub_01 <- lapply(models_01, function(x){
  reduce2traits_models(x, trait_names)
})
final_species <- lapply(seq(length(models_sub_01[[1]])), function(x){
  models_sub_01[[1]][[x]][[1]]$response
})
final_species <- unlist(final_species)

models_sub_02 <- lapply(models_02, function(x){
  reduce2traits_models(x, trait_names)
})

final_species <- lapply(seq(length(models_sub_02[[1]])), function(x){
  models_sub_02[[1]][[x]][[1]]$response
})
final_species <- unlist(final_species)

# TRTE
del_nbr <- lapply(seq(length(orthoptera_trte)), function(x){
  species <- orthoptera_trte[[x]][[1]]$training$RESPONSE
  if(species %in% final_species){
    del_nbr <- x
  } else {
    del_nbr <- NA
  }
  return(del_nbr)
})
del_nbr <- unlist(del_nbr)
del_nbr <- del_nbr[!is.na(del_nbr)]
orthoptera_trte_sub <- orthoptera_trte[del_nbr]

# OBSV
obsv_gpm_subset <- obsv_gpm
obsv_gpm_subset@meta$input$RESPONSE <- obsv_gpm_subset@meta$input$RESPONSE[obsv_gpm_subset@meta$input$RESPONSE %in% final_species]
obsv_gpm_subset@data$input <- obsv_gpm_subset@data$input[colnames(obsv_gpm_subset@data$input) %in% unlist(obsv_gpm_subset@meta$input)]

# TSTAT
tstat_sub_01 <- lapply(tstat_01, function(x){
  reduce2traits_tstat(x, final_species)
  })

tstat_mean_sub_01 <- lapply(tstat_mean_01, function(x){
  tstat_mean <- x[x$Response %in% final_species, ]
})

tstat_merged_sub_01 <- tstat_merged_01[tstat_merged_01$Response %in% final_species, ]
tstat_merged_sub_01 <- droplevels(tstat_merged_sub_01)

tstat_merged_sub_inv_01 <- tstat_merged_01[!tstat_merged_01$Response %in% final_species, ]
tstat_merged_sub_inv_01 <- droplevels(tstat_merged_sub_inv_01)
# unique(tstat_merged_sub_inv_01$Response)

tstat_sub_02 <- lapply(tstat_02, function(x){
  reduce2traits_tstat(x, final_species)
})

tstat_mean_sub_02 <- lapply(tstat_mean_02, function(x){
  tstat_mean <- x[x$Response %in% final_species, ]
})

tstat_merged_sub_02 <- tstat_merged_02[tstat_merged_02$Response %in% final_species, ]
tstat_merged_sub_02 <- droplevels(tstat_merged_sub_02)

tstat_merged_sub_inv_02 <- tstat_merged_02[!tstat_merged_02$Response %in% final_species, ]
tstat_merged_sub_inv_02 <- droplevels(tstat_merged_sub_inv_02)


# Create observational dataset
dfx <- lapply(orthoptera_trte_sub, function(x){
  dfy <- lapply(x, function(y){
    c(y$training$SAMPLES, y$testing$SAMPLES)
  })
  unique(unlist(dfy))
})
rows <- unique(unlist(dfx))

obsv_samples <- obsv_gpm_subset@data$input[rows, ]

obsv_samples$date_observation <- as.Date(obsv_samples$date_observation, format="%Y-%j")
obsv_samples$year <- as.numeric(format(obsv_samples$date_observation, "%Y"))

# save(models_sub_01, orthoptera_trte_sub, obsv_gpm_subset,
#      tstat_sub_01, tstat_mean_sub_01, tstat_merged_sub_01, tstat_merged_sub_inv_01, obsv_samples,
#      trait_names, final_species,
#      file = paste0(filepath_results, "orthoptera_figures_trait_subsets_01.RData"))

# save(models_sub_02, orthoptera_trte_sub, obsv_gpm_subset,
#      tstat_sub_02, tstat_mean_sub_02, tstat_merged_sub_02, tstat_merged_sub_inv_02, obsv_samples,
#      trait_names, final_species,
#      file = paste0(filepath_results, "orthoptera_figures_trait_subsets_02.RData"))
