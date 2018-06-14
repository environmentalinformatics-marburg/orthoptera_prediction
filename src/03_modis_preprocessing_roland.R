# Preprocess dataset by combining field and satellite observations
# Thomas Nauss

if(Sys.info()["sysname"] == "Windows"){
  source("F:/analysis/orthoptera/orthoptera_prediction/src/00_set_environment.R")
} else {
  source("/media/tnauss/myWork/analysis/orthoptera/orthoptera_prediction/src/00_set_environment.R")
}


# Preparation
# modis_mod_plots <- readRDS(file = paste0(path_results, "modis_mod_plots.rds"))
# modis_myd_plots <- readRDS(file = paste0(path_results, "modis_myd_plots.rds"))
# modis_mod <- as.data.frame(modis_mod_plots)
# modis_myd <- as.data.frame(modis_myd_plots)
# saveRDS(modis_mod, file = paste0(path_results, "modis_mod_df.rds"))
# saveRDS(modis_myd, file = paste0(path_results, "modis_myd_df.rds"))


modis_mod <- readRDS(file = paste0(path_results, "modis_mod_df.rds"))

# Columns in Terra-MODIS dataset (modis_mod)
col_selector <- which(names(modis_mod) == "plot")

col_meta <- c(seq(which(names(modis_mod) == "date"), 
                  which(names(modis_mod) == "elevation")),
              seq(which(names(modis_mod) == "rainfall"), 
                  which(names(modis_mod) == "year")),
              seq(which(names(modis_mod) == "ID"), 
                  which(names(modis_mod) == "veg.nr")),
              seq(which(names(modis_mod) == "coordW"), 
                  which(names(modis_mod) == "coordN")))

col_diversity <- c(which(names(modis_mod) == "rich.plants"),
                   seq(which(names(modis_mod) ==  "specrich"), 
                       which(names(modis_mod) == "shannon.nw")))

col_species <- seq(which(names(modis_mod) == "Abisares.depressus"),
                   which(names(modis_mod) == "Zonocerus.elegans"))

col_precitors <- seq(which(names(modis_mod) == "modis_sur_refl_b01_mean"),
                     which(names(modis_mod) == "ndvi_glcm_correlation_var"))

