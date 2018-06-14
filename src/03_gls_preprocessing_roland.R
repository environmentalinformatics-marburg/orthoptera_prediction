# Preprocess dataset by combining field and satellite observations
# Thomas Nauss

if(Sys.info()["sysname"] == "Windows"){
  source("F:/analysis/orthoptera/orthoptera_prediction/src/00_set_environment.R")
} else {
  source("/media/tnauss/myWork/analysis/orthoptera/orthoptera_prediction/src/00_set_environment.R")
}

# Preparation
# gls_2000_plots <- readRDS(file = paste0(path_results, "gls_2000_plots.rds"))
# obsv_gls <- as.data.frame(gls_2000_plots)
# saveRDS(obsv_gls, file = paste0(path_results, "obsv_gls_df.rds"))

# Dataset for Roland
obsv_gls <- readRDS(obsv_gls, file = paste0(path_results, "obsv_gls_df.rds"))

col_selector <- which(names(obsv_gls) == "plot")
    
col_meta <- c(seq(which(names(obsv_gls) == "date"), 
                  which(names(obsv_gls) == "elevation")),
              seq(which(names(obsv_gls) == "rainfall"), 
                  which(names(obsv_gls) == "year")),
              seq(which(names(obsv_gls) == "ID"), 
                  which(names(obsv_gls) == "veg.nr")),
              seq(which(names(obsv_gls) == "coordW"), 
                  which(names(obsv_gls) == "coordN")))

col_diversity <- c(which(names(obsv_gls) == "rich.plants"),
                   seq(which(names(obsv_gls) ==  "specrich"), 
                       which(names(obsv_gls) == "shannon.nw")))
    
col_species <- seq(which(names(obsv_gls) == "Abisares.depressus"),
                   which(names(obsv_gls) == "Zonocerus.elegans"))
    
col_precitors <- seq(which(names(obsv_gls) == "gls2000.1_mean"),
                     which(names(obsv_gls) == "ndvi_glcm_correlation.r25_var"))
    
