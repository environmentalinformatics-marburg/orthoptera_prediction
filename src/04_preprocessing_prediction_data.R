# Preprocess dataset by combining field and satellite observations
# Thomas Nauss

if(Sys.info()["sysname"] == "Windows"){
  source("F:/analysis/orthoptera/orthoptera_prediction/src/00_set_environment.R")
} else {
  source("/media/tnauss/myWork/analysis/orthoptera/orthoptera_prediction/src/00_set_environment.R")
}

compute <- FALSE

# Merge GLS2000 data with preprocessed orthoptera observations -----------------
if(compute){
  gls_2000_wgs <- readRDS(file = paste0(path_results, "gls_2000_wgs.rds"))
  gls_2000_arc <- readRDS(file = paste0(path_results, "gls_2000_arc.rds"))
  
  obsv_shp_wgs <- readRDS(file = paste0(path_results, "obsv_shp_wgs.rds"))
  obsv_shp_arc <- readRDS(file = paste0(path_results, "obsv_shp_arc.rds"))
  
  cprj <- c("wgs", "arc")
  gls_2000_plots <- lapply(cprj, function(prj){
    if(prj == "wgs"){
      gls_2000 <- gls_2000_wgs
      obsv_shp <- spTransform(obsv_shp_wgs, crs(gls_2000[[1]]))
    } else {
      gls_2000 <- gls_2000_arc
      obsv_shp <- spTransform(obsv_shp_arc, crs(gls_2000[[1]]))
    }
    plots <- extractFromRasterSnips(raster = gls_2000, spatial = obsv_shp, 
                                    selector = "plot", buffer = 50)
  })
  names(gls_2000_plots) <- cprj
  saveRDS(gls_2000_plots, file = paste0(path_results, "gls_2000_plots.rds"))
} else {
  gls_2000_plots <- readRDS(file = paste0(path_results, "gls_2000_plots.rds"))
}

# Prepare data set used for remote sensing prediction study --------------------
cprj <- c("wgs", "arc")
obsv_gpm <- lapply(cprj, function(prj){
  obsv <- as.data.frame(gls_2000_plots[[prj]])
  col_selector <- which(names(obsv) == "plot")
  
  col_meta <- c(seq(which(names(obsv) == "date"), 
                    which(names(obsv) == "elevation")),
                seq(which(names(obsv) == "rainfall"), 
                    which(names(obsv) == "year")),
                seq(which(names(obsv) == "ID"), 
                    which(names(obsv) == "veg.nr")),
                seq(which(names(obsv) == "coordW"), 
                    which(names(obsv) == "coordN")))
  
  col_diversity <- c(which(names(obsv) == "rich.plants"),
                     seq(which(names(obsv) ==  "specrich"), 
                         which(names(obsv) == "shannon.nw")))
  
  col_species <- seq(which(names(obsv) == "Abisares.depressus"),
                     which(names(obsv) == "Zonocerus.elegans"))
  
  col_precitors <- seq(which(names(obsv) == "gls2000.1"),
                       which(names(obsv) == "glcm_correlation_sd"))
  
  if((length(unique(c(col_selector, col_meta, col_diversity, 
                        col_species, col_precitors))) == ncol(obsv))){
    obsv[, col_species][obsv[, col_species] > 0] <- "yes"
    obsv[, col_species][obsv[, col_species] == "0"] <- "no"
    for(i in col_species){
      obsv[, i] <- as.factor(obsv[, i])
    }
    meta <- createGPMMeta(obsv, type = "input",
                          selector = col_selector, 
                          response = col_species, 
                          independent = col_precitors, 
                          meta = c(col_meta, col_diversity))
    obsv <- gpm(obsv, meta, scale = TRUE)
  } else {
    obsv <- NULL
  }
  return(obsv)
})
names(obsv_gpm) <- cprj

saveRDS(obsv_gpm[["wgs"]], file = paste0(path_results, "obsv_gpm_wgs.rds"))
saveRDS(obsv_gpm[["arc"]], file = paste0(path_results, "obsv_gpm_arc.rds"))




