# Preprocess dataset by combining field and satellite observations
# Thomas Nauss

if(Sys.info()["sysname"] == "Windows"){
  source("D:/orthoptera/orthoptera_prediction/src/00_set_environment.R")
} else {
  source("/media/tnauss/myWork/analysis/orthoptera/orthoptera_prediction/src/00_set_environment.R")
}


# Prepare orthoptera traits ----------------------------------------------------
obsv_traits <- read_excel(paste0(path_obsv, "Orthoptera_traits_10_11_2015.xls"),
                          sheet = "Daten", skip = 0)  
cn <- c("Page", "Number", "Name", "Sex", "Body_length", "Pronot_length",
        "Pronot_width", "Eye_crosssection", "Eye_absolute", "Ovipositor", 
        "Wings", "Femur_posterior", "Tibia_posterior", "Femur_medial", 
        "Femur_anterior", "Month", "Year", "HS", "Elevation")
cnt <- c("Body_length", "Pronot_length",
         "Pronot_width", "Eye_crosssection", "Eye_absolute", "Ovipositor", 
         "Wings", "Femur_posterior", "Tibia_posterior", "Femur_medial", 
         "Femur_anterior")
colnames(obsv_traits) <- cn
obsv_traits <- as.data.frame(obsv_traits)

for(n in cnt){
  obsv_traits[, n] <- as.numeric(obsv_traits[, n])
}
obsv_traits <- obsv_traits[complete.cases(obsv_traits[, cnt]), ]
obsv_traits$Name <- gsub(x = obsv_traits$Name, pattern = " ", 
                         replacement = ".")


obsv_traits$Name[which(obsv_traits$Name == "Acorypha.lacticosta.")] <- "Acorypha.laticosta"
obsv_traits$Name[which(obsv_traits$Name == "Acanthacris.ruficornis")] <- "Acanthacris.ruficornis.ruficornis"
obsv_traits$Name[which(obsv_traits$Name == "Cataloipus.oberthuri.oberthuri" )] <- "Cataloipus.oberthuri"
obsv_traits$Name[which(obsv_traits$Name == "Cyrtacanthacris.tatarica.tatarica")] <- "Cyrtacanthacris.tatarica"
obsv_traits$Name[which(obsv_traits$Name == "Gastrimargus.africanus.africanus")] <- "Gastrimargus.africanus"
obsv_traits$Name[which(obsv_traits$Name == "Gastrimargus.verticalis..")] <- "Gastrimargus.verticalis"
obsv_traits$Name[which(obsv_traits$Name == "Gymnobothrus.levipes..levipes")] <- "Gymnobothroides.levipes"
obsv_traits$Name[which(obsv_traits$Name == "Gymnobothrus.levipes.levipes")] <- "Gymnobothroides.levipes" 
obsv_traits$Name[which(obsv_traits$Name == "Gymnobothrus.temporalis.flexuosus")] <- "Gymnobothrus.flexuosus"
obsv_traits$Name[which(obsv_traits$Name == "Gymnobothrus.temporalis.temporalis")] <- "Gymnobothrus.temporalis"
obsv_traits$Name[which(obsv_traits$Name == "Oedaleus.sengalensis")] <- "Oedaleus.senegalensis"
obsv_traits$Name[which(obsv_traits$Name == "Pycnodictya.galinieri.galineri")] <- "Pycnodictya.galinieri"

obsv_traits <- aggregate(obsv_traits[, cnt], by = list(obsv_traits$Name), 
                              FUN = mean)
names(obsv_traits)[1] <- "Name"

# Compute PCA
obsv_traits_pca <- prcomp(obsv_traits[, cnt], center = TRUE, scale = FALSE)
biplot(obsv_traits_pca,  choices = c(1,2))

obsv_traits <- cbind(obsv_traits, obsv_traits_pca$x)

saveRDS(obsv_traits, paste0(path_results, "obsv_traits.rds"))
