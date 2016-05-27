# Preprocess MODIS data
# devtools::install_github("MatMatt/MODIS", ref = "develop")
# setRepositories() 
# install.packages(c(' mapdata', 'ptw '),dependencies=TRUE)
library(rgdal)
library(MODIS)

if(Sys.info()["sysname"] == "Windows"){
  filepath_base <- "G:/analysis/orthoptera/data/"
} else {
  filepath_base <- "/media/tnauss/myWork/analysis/orthoptera/data/"
}

localArcPath <- paste0(filepath_base, "modis/modis_arc/")
outDirPath <- paste0(filepath_base, "modis/processed/")
cropPath <-paste0(filepath_base, "modis/croped/")

obsv_filepath <- paste0(filepath_base, "orthoptera/")

rasterOptions(tmpdir= paste0(filepath_base, "temp/"))

# Read observation data for dates
obsv <- read.table(paste0(obsv_filepath, "lvl0300_biodiversity_data.csv"), 
                   header = TRUE, sep = ";", dec = ",")
# head(obsv)
dates <- gsub("-", "", substring(unique(obsv$date_nocloud), 1, 8), "-")

# Set processing options
MODISoptions(localArcPath = localArcPath, outProj = "+init=epsg:32737",
             outDirPath = outDirPath, MODISserverOrder = c("LAADS", "LPDAAC"),
             gdalPath = "C:/OSGeo4W64/bin")


# Download MODIS tiles
modis_product <- c("MOD09GA", "MOD09A1", "MOD09Q1", "MOD11A1",
                   "MYD09GA", "MYD09A1", "MYD09Q1", "MYD11A1")

for(pdct in modis_product){
  for(dts in dates){
    job_name <- paste0(pdct, "_TZS")
    collection <- getCollection(pdct, forceCheck = TRUE)
    runGdal(pdct, collection = collection, job = job_name, 
            tileH = 21, tileV = 9,
            begin = dts, end = dts)
  }
}

  
# Post-processing
obsv_buffer <- readOGR(paste0(
  obsv_filepath, "lvl0300_biodiversity_data_utm37s_buffer_dissolve_single.shp"),
  "lvl0300_biodiversity_data_utm37s_buffer_dissolve_single")

modis_product <- c("MOD09GA", "MOD09A1", "MOD09Q1",
                   "MYD09GA", "MYD09A1", "MYD09Q1")
for(pdct in modis_product){
  job_name <- paste0(pdct, "_TZS")
  modis_files <- list.files(paste0(outDirPath, job_name),
                            pattern = glob2rx("*sur_refl_b0*tif"), 
                            full.names = TRUE)
  modis_data <- stack(modis_files)
  obsv_buffer_tf <- spTransform(obsv_buffer, projection(modis_data))
  modis_data_crop <- crop(modis_data, obsv_buffer_tf)     
  
  modis_outpath <- paste0(outDirPath, job_name, "_croped")
  dir.create(modis_outpath,  recursive = FALSE)
  modis_outfiles <- paste0(modis_outpath, "/", basename(modis_files))
  
  writeRaster(modis_data_crop, filename = modis_outfiles, 
              format = "GTiff", bylayer = TRUE)
  
}
