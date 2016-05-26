# Preprocess MODIS data
# devtools::install_github("MatMatt/MODIS", ref = "develop")
# setRepositories() 
# install.packages(c(' mapdata', 'ptw '),dependencies=TRUE)
library(rgdal)
library(MODIS)

if(Sys.info()["sysname"] == "Windows"){
  filepath_base <- "E:/analysis/orthoptera/"
} else {
  filepath_base <- "/media/tnauss/myWork/analysis/orthoptera/"
}

localArcPath <- paste0(filepath_base, "modis/modis_arc/")
outDirPath <- paste0(filepath_base, "modis/processed/")
cropPath <-paste0(filepath_base, "modis/croped/")

obsv_filepath <- paste0(filepath_base, "data/orthoptera/")

rasterOptions(tmpdir= paste0(filepath_base, "temp/"))

# Read observation data for dates
obsv <- read.table(paste0(obsv_filepath, "lvl0300_biodiversity_data.csv"), 
                   header = TRUE, sep = ";", dec = ",")
# head(obsv)
dates <- gsub("-", "", substring(unique(obsv$date_nocloud), 1, 8), "-")

# Set processing options
MODISoptions(localArcPath = localArcPath, outProj = "+init=epsg:32737",
             outDirPath = outDirPath, MODISserverOrder = c("LAADS", "LPDAAC"))
# , gdalPath = "C:/OSGeo4W64/bin")


# Download MODIS tiles
modis_product <- c("MOD09GA", "MOD09A1", "MOD09Q1", "MOD11A1",
                   "MYD09GA", "MYD09A1", "MYD09Q1", "MYD11A1")
collection <- getCollection("MYD09GA", forceCheck = TRUE)

for(pdct in modis_product){
  for(dts in dates){
    job_name <- paste0(pdct, "_TZS")
    collection <- getCollection(pdct, forceCheck = FALSE)
    collection[[1]] <- "006"
    runGdal(pdct, collection = collection, job = job_name, 
            tileH = 15, tileV = 7,
            begin = dts, end = dts)
  }
}


modis_product <- "MOD09A1"
job_name <- paste0(modis_product, "_CPV")

collection <- getCollection(modis_product, forceCheck = FALSE)
collection[[1]] <- "006"

runGdal(modis_product, collection = collection, job = job_name, 
        tileH = 15, tileV = 7,
        begin = "2004274", end = "2004276")






# Post-processing
# Get country boundaries
# country_shp <- getData(country = "CPV", level = 0, path = shpPath)
country_shp <- readOGR(paste0(shpPath, "CAP_admin_UTM26_SHP/CAP_UTM26.shp"), 
                       layer = "CAP_UTM26")
projection(country_shp) <- 
  CRS("+proj=utm +zone=26 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ")
country_shp <- spTransform(country_shp, CRS = CRS("+init=epsg:32626"))
# plot(country_shp)

modis_files <- list.files(paste0(getOption("MODIS_outDirPath"), "/", job_name),
                          pattern = glob2rx("*.tif"), full.names = TRUE)
modis_data <- stack(modis_files)
modis_data_crop <- crop(modis_data, country_shp)     
writeRaster(modis_data_crop, paste0(cropPath, "MODIS.tif"), bylayer = TRUE)


# lst <- lapply(c("NDVI", "VI_Quality"), function(i) {
#   fls <- list.files(paste0(getOption("MODIS_outDirPath"), "/MOD09A1_CPV"), 
#                     pattern = paste0(i, ".tif$"), full.names = TRUE)
#   rst <- stack(fls)
#   rst_crp <- crop(rst, shp)     
#   if (i == "NDVI") 
#     rst_crp <- rst_crp * 0.0001
#   return(rst_crp)
# })

# ## install and load 'MODIS' version 0.10-18 (assuming that the older version of 
# ## 'MODIS' is located inside your current working directory)
# detach("package:MODIS", unload = TRUE)
# install.packages("MODIS_0.10-18.tar.gz", repos = NULL, type = "source")
# library(MODIS)