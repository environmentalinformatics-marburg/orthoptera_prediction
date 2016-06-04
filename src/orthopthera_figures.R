# Libraries --------------------------------------------------------------------
# library(gpm)
library(ggplot2)
library(ggmap)
library(mapview)
library(rgdal)
library(raster)
library(sp)
library(reshape2)
library(quantreg)
library(latticeExtra)

# Read data --------------------------------------------------------------------
if(Sys.info()["sysname"] == "Windows"){
  filepath_base <- "E:/analysis/orthoptera/data/"
} else {
  filepath_base <- "/media/tnauss/myWork/analysis/orthoptera/data/"
}

filepath_results <- paste0(filepath_base, "rdata/")
filepath_obsv <- paste0(filepath_base, "orthoptera/")
filepath_traits <- paste0(filepath_base, "traits/")
filepath_figures <- paste0(filepath_base, "figures/")


load(paste0(filepath_results, "orthoptera_prediction_obsv_gpm.RData"))
load(paste0(filepath_results, "orthoptera_prediction_prevalence.RData"))
load(paste0(filepath_results, "orthoptera_prediction_orthoptera_resamples.RData"))
load(paste0(filepath_results, "orthoptera_prediction_orthoptera_trte.RData"))
load(paste0(filepath_results, "orthoptera_prediction_models_rf_2016-05-28_rfe.RData"))
dem <- raster(paste0(filepath_obsv, "DEM_UTM37S_WGS84_30m_Hemp.tif"))
traits <- read.table(paste0(filepath_traits, "Orthoptera_10_11_2015.csv"),
                     header = TRUE, sep = ";", dec = ".")


# Reduce model dataset to species for which traits are available----------------
trait_names <- gsub("\\ ", ".", traits$Name)

trait_names[!(trait_names %in% prevalence$RESPONSE)]
trait_names[trait_names == "Cyrtacanthacris.tatarica.tatarica"] <- "Cyrtacanthacris.tatarica"
traits$Name <- trait_names

trait_names <- unique(traits$Name)
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

final_species <- lapply(seq(length(models_subset)), function(x){
  models_subset[[x]][[1]]$response
})
final_species <- unlist(final_species)

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
orthoptera_trte_subset <- orthoptera_trte[del_nbr]

obsv_gpm_subset <- obsv_gpm
obsv_gpm_subset@meta$input$RESPONSE <- obsv_gpm_subset@meta$input$RESPONSE[obsv_gpm_subset@meta$input$RESPONSE %in% final_species]
obsv_gpm_subset@data$input <- obsv_gpm_subset@data$input[colnames(obsv_gpm_subset@data$input) %in% unlist(obsv_gpm_subset@meta$input)]


# Calculate basic model statistics ---------------------------------------------
var_imp <- compVarImp(models_subset, scale = FALSE)
var_imp_scale <- compVarImp(models_subset, scale = TRUE)
var_imp_plot <- plotVarImp(var_imp)
var_imp_heat <- plotVarImpHeatmap(var_imp_scale, xlab = "Species", ylab = "Band")
tstat <- compContTests(models_subset, mean = TRUE)
tstat_mean <- merge(tstat[[1]], prevalence, by.x = "Response", by.y = "RESPONSE")
tstat_mean[order(tstat_mean$Kappa_mean, decreasing = TRUE),]
# save(var_imp, var_imp_scale, tstat, tstat_mean, 
#      file = paste0(filepath_results, "orthoptera_figures_var_imp_tstat.RData"))
# load(paste0(filepath_results, "orthoptera_figures_var_imp_tstat.RData"))


dfx <- lapply(orthoptera_trte_subset, function(x){
  dfy <- lapply(x, function(y){
    c(y$training$SAMPLES, y$testing$SAMPLES)
  })
  unique(unlist(dfy))
})
rows <- unique(unlist(dfx))

sample <- obsv_gpm_subset@data$input[rows, ]

sample$date_observation <- as.Date(sample$date_observation, format="%Y-%j")
sample$year <- as.numeric(format(sample$date_observation, "%Y"))



# sample_names <- gsub("\\.", " ", colnames(sample)) 
# sample_names[sample_names == "Cyrtacanthacris tatarica"] <-
#   "Cyrtacanthacris tatarica tatarica"
# sample_names[sample_names == "Gymnobothrus flexuosus"] <-
#   "Gymnobothrus temporalis flexuosus"

# sample[, colnames(sample) %in% obsv_gpm_subset@meta$input$META]
# which(sample_names %in% unique(traits$Name))

# tstat_mean$Response <- gsub("\\.", " ", tstat_mean$Response) 
# tstat_mean$Response[tstat_mean$Response == "Cyrtacanthacris tatarica"] <-
#   "Cyrtacanthacris tatarica tatarica"
# tstat_mean$Response[tstat_mean$Response == "Gymnobothrus flexuosus"] <-
#   "Gymnobothrus temporalis flexuosus"
# tstat_mean$Response_SHORT <- substr(tstat_mean$Response, 1, 4)

# Adjust prevalence information
# prevalence$RESPONSE <- gsub("\\.", " ", prevalence$RESPONSE) 
# prevalence$RESPONSE[prevalence$RESPONSE == "Cyrtacanthacris tatarica"] <-
#   "Cyrtacanthacris tatarica tatarica"
# prevalence$RESPONSE[prevalence$RESPONSE == "Gymnobothrus flexuosus"] <-
#   "Gymnobothrus temporalis flexuosus"



# Create maps ------------------------------------------------------------------
map_extent <- get_map(location = c(36.93865,
                                   -3.454621,
                                   37.76235,
                                   -2.775392),
                      scale = "auto",
                      maptype = "satellite",
                      color = "bw",
                      source = "google")

map <- ggmap(map_extent, 
             extent = "normal",
             maprange = TRUE)

obs <- geom_point(aes(x = lon,
                      y = lat,
                      size = 2,
                      colour = nr.of.species),
                  show.legend = FALSE,
                  data = sample)

colourscale <- scale_colour_gradient(low = "white", 
                                     high = "darkgreen", 
                                     name = "Species number",
                                     limits=c(0, 30))

style_map_01 <- theme(legend.background = element_rect(colour = "black"),
                      plot.title = element_text(size = 20))

labels_map_01 <- labs(title = "Orthoptera observations 2002 - 2012")

png(paste0(filepath_figures, "fig_01_map.png"), 
    width = 1024 * 6, 
    height = 748 * 6, 
    units = "px", 
    res = 600)
map + obs + colourscale + labels_map_01 + style_map_01 + theme_bw()
graphics.off()


obs <- geom_point(aes(x = lon,
                      y = lat,
                      size = 1,
                      colour = nr.of.species),
                  show.legend = FALSE,
                  data = sample)

png(paste0(filepath_figures, "fig_a01_maps.png"), 
    width = 1024 * 6, 
    height = 748 * 6, 
    units = "px", 
    res = 600)
map + obs + colourscale + labels_map_01 + style_map_01 + theme_bw() + facet_wrap(~ year)
graphics.off()


# Elevation statistics ---------------------------------------------------------
sample_utm <- sample
coordinates(sample_utm) <- ~lon+lat
projection(sample_utm) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs ")
sample_utm <- spTransform(sample_utm, projection(dem))

sample_utm@data$masl <- extract(dem, sample_utm)

summary(sample_utm@data$masl)
quantile(sample_utm@data$masl, probs = seq(0, 1, 0.1), na.rm = TRUE)


# NDVI statistics --------------------------------------------------------------
sample_utm@data$NDVI <- 
  (sample_utm@data$MYD09GA_sur_refl_b02_1 - sample_utm@data$MYD09GA_sur_refl_b01_1) / 
  (sample_utm@data$MYD09GA_sur_refl_b02_1 + sample_utm@data$MYD09GA_sur_refl_b01_1)

summary(sample_utm@data$NDVI)
quantile(sample_utm@data$NDVI, probs = seq(0, 1, 0.1), na.rm = TRUE)

gsub("\\.", " ", prevalence$RESPONSE) 

quantile(sample_utm@data$diff_days_nocloud, seq(0, 1, 0.05))


# Prediction performance -------------------------------------------------------
plot_scores <- plotClassPerformance(tstat, scores = c("Kappa", "ETS"))

png(paste0(filepath_figures, "fig_02_scores.png"), 
    width = 1024 * 6, 
    height = 748 * 6, 
    units = "px", 
    res = 600)
plot_scores
graphics.off()


# Test dependency of prediction accuracy ---------------------------------------
test <- lapply(models_subset, function(x){
  act_test <- lapply(x, function(y){
    if(inherits(y$model, "try-error")){
      NULL
    } else {
      df <- data.frame(TARGET = y$response,
                       RESPONSE = y$testing$RESPONSE,
                       PREDICTION = y$testing$PREDICTED$pred,
                       NDVI = (y$testing$INDEPENDENT$MYD09GA_sur_refl_b02_1 - y$testing$INDEPENDENT$MYD09GA_sur_refl_b01_1) / 
                         (y$testing$INDEPENDENT$MYD09GA_sur_refl_b02_1 + y$testing$INDEPENDENT$MYD09GA_sur_refl_b01_1),
                       y$testing$META)
    }  
  })
  do.call("rbind", act_test)
})
test <- do.call("rbind", test)    

test_utm <- test
coordinates(test_utm) <- ~lon+lat
projection(test_utm) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs ")
test_utm <- spTransform(test_utm, projection(dem))

test_utm@data$masl <- extract(dem, test_utm)
test_utm@data$performance <- test_utm@data$RESPONSE == test_utm@data$PREDICTION
test_utm@data$performance_num <- as.numeric(test_utm@data$performance)

test_utm@data <- merge(test_utm@data, prevalence, by.x = "TARGET", by.y = "RESPONSE")

ndvi_mean <- aggregate(test_utm@data$NDVI, by = list(test_utm@data$TARGET), FUN = "mean")
colnames(ndvi_mean) <- c("Response", "NDVI")

asl_mean <- aggregate(test_utm@data$asl, by = list(test_utm@data$TARGET), FUN = "mean")
colnames(asl_mean) <- c("Response", "asl")

rainfall_mean <- aggregate(test_utm@data$rainfall, by = list(test_utm@data$TARGET), FUN = "mean")
colnames(rainfall_mean) <- c("Response", "rainfall")

species_mean <- aggregate(test_utm@data$nr.of.species, by = list(test_utm@data$TARGET), FUN = "mean")
colnames(species_mean) <- c("Response", "nr.of.species")

tstat_mean <- merge(tstat[[1]], prevalence, by.x = "Response", by.y = "RESPONSE")
tstat_mean <- merge(tstat_mean, ndvi_mean, by = "Response")
tstat_mean <- merge(tstat_mean, asl_mean, by = "Response")
tstat_mean <- merge(tstat_mean, rainfall_mean, by = "Response")
tstat_mean <- merge(tstat_mean, species_mean, by = "Response")

nv <- function(x){
  (x - min(x)) * (1 - 0) / (max(x) - min(x)) + 0
}

tstat_mean$OCCURENCE_scaled <- nv(tstat_mean$OCCURENCE)
tstat_mean$NDVI_scaled <- nv(tstat_mean$NDVI)
tstat_mean$asl_scaled <- nv(tstat_mean$asl)
tstat_mean$rainfall_scaled <- nv(tstat_mean$rainfall)
tstat_mean$nr.of.species_scaled <- nv(tstat_mean$nr.of.species)

tstat_mean_melt <- melt(tstat_mean[, colnames(tstat_mean) %in%
                                     c("Kappa_mean", "OCCURENCE_scaled", 
                                       "NDVI_scaled", "asl_scaled", 
                                       "rainfall_scaled", 
                                       "nr.of.species_scaled")], 
                        id.vars = "Kappa_mean")

png(paste0(filepath_figures, "fig_03_scores_dependency.png"), 
    width = 1024 * 6, 
    height = 748 * 6, 
    units = "px", 
    res = 600)
ggplot(data = tstat_mean_melt, aes(x = value, y = Kappa_mean)) + 
  geom_point() + 
  geom_smooth() + 
  facet_wrap(~variable,  scales = "free_y")
graphics.off()

# ggplot(data = tstat_mean, aes(x = OCCURENCE, y = Kappa_mean)) + geom_point() + geom_smooth()
# ggplot(data = tstat_mean, aes(x = NDVI, y = Kappa_mean)) + geom_point() + geom_smooth()
# ggplot(data = tstat_mean, aes(x = asl, y = Kappa_mean)) + geom_point() + geom_smooth()
# ggplot(data = tstat_mean, aes(x = rainfall, y = Kappa_mean)) + geom_point() + geom_smooth()
# ggplot(data = tstat_mean, aes(x = nr.of.species,  y= Kappa_mean)) + geom_point() + geom_smooth()

test_utm_melt <- 
  melt(test_utm@data[, colnames(test_utm@data) %in% 
                       c("TARGET", "performance", "NDVI", "asl", "rainfall", "nr.of.species")],
       id.var = c("TARGET", "performance"))

png("figures/fig_04_scores_dependency_boxplots.png", 
    width = 1024 * 6, 
    height = 748 * 6, 
    units = "px", 
    res = 600)
ggplot(data = test_utm_melt, aes(x = TARGET, y = value, fill = performance)) + 
  geom_boxplot(notch = TRUE) + facet_wrap(~variable,  scales = "free_y") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0))
graphics.off()


# png("figures/fig_03a_scores_dependency_01.png", 
#     width = 1024 * 6, 
#     height = 748 * 6, 
#     units = "px", 
#     res = 600)
# ggplot(data = test_utm@data, aes(x = TARGET, y = nr.of.species, fill = performance)) + 
#   geom_boxplot(notch = TRUE) +
#   theme(axis.text.x=element_text(angle = -90, hjust = 0))
# graphics.off()
# 
# png("figures/fig_03b_scores_dependency_01.png", 
#     width = 1024 * 6, 
#     height = 748 * 6, 
#     units = "px", 
#     res = 600)
# ggplot(data = test_utm@data, aes(x = TARGET, y = masl, fill = performance)) + 
#   geom_boxplot(notch = TRUE) +
#   theme(axis.text.x=element_text(angle = -90, hjust = 0))
# graphics.off()


# Variable importance ----------------------------------------------------------
png("figures/fig_05_variable_importance.png", 
    width = 1024 * 6, 
    height = 748 * 6, 
    units = "px", 
    res = 600)
plotVarImpHeatmap(var_imp, xlab = "Species", ylab = "Band")
graphics.off()


# Prediction dependency from functional diversity ------------------------------
clrs_hcl <- function(n) {
  hcl(h = seq(0, 260, length.out = n), 
      c = 60, l = 50, 
      fixup = TRUE)
}

scatter3D <- function(x, y, z, 
                      color = c("blue", "lightgreen", "gold", "red"),
                      alpha = 0.5,
                      ...) {
  stopifnot(require("latticeExtra"))
  limits <- range(na.exclude(x), na.exclude(y))
  limits[1] <- limits[1] - 0.1 * limits[2]
  limits[2] <- limits[2] + 0.1 * limits[2]
  zmin <- min(na.exclude(z)) - 0.1 * min(na.exclude(z))
  zmax <- max(na.exclude(z)) + 0.1 * max(na.exclude(z))
  color <- colorRampPalette(color)
  plot <- xyplot(y ~ x, z = z, xlim = limits, ylim = limits,     
                 panel = function() {
                   panel.grid(h = -1, v = -1)
                   panel.levelplot.points(x = x, y = y, z = z)                                 
                   panel.abline(a = 0, b = 1, col = "red", lty = 2, 
                                lwd = 1.5)
                 },
                 aspect = "iso")
  plot <- update(plot,
                 par.settings = custom.theme(region = color(200), alpha = alpha),
                 legend = list(top = list(fun = "draw.colorkey", 
                                          args = list(list(space = "top",
                                                           at = zmin:zmax,
                                                           width = 1,
                                                           col = color(200))))),
                 ...)
  return(plot)
}

# Adjust test statistics
tstat_mean$Response <- gsub("\\.", " ", tstat_mean$Response) 
tstat_mean$Response[tstat_mean$Response == "Cyrtacanthacris tatarica"] <-
  "Cyrtacanthacris tatarica tatarica"
tstat_mean$Response[tstat_mean$Response == "Gymnobothrus flexuosus"] <-
  "Gymnobothrus temporalis flexuosus"
tstat_mean$Response_SHORT <- substr(tstat_mean$Response, 1, 4)

# Adjust prevalence information
prevalence$RESPONSE <- gsub("\\.", " ", prevalence$RESPONSE) 
prevalence$RESPONSE[prevalence$RESPONSE == "Cyrtacanthacris tatarica"] <-
  "Cyrtacanthacris tatarica tatarica"
prevalence$RESPONSE[prevalence$RESPONSE == "Gymnobothrus flexuosus"] <-
  "Gymnobothrus temporalis flexuosus"

traits <- read.table(paste0(filepath_traits, "Orthoptera_10_11_2015.csv"),
                     header = TRUE, sep = ";", dec = ".")
traits <- traits[traits$Name %in% tstat_mean$Response, ]
traits[, 5:15] <- log(traits[, 5:15])
traits_mean <- aggregate(traits[, 5:15], by = list(traits$Name), FUN = mean,
                         na.rm = TRUE)
# Combine model statistics and traits and analyse information
expl <- cbind(tstat_mean[tstat_mean$Response %in% unique(traits_mean$Group.1),],
              traits_mean, by.x = "Response", by.y = "Group.1")
expl <- cbind(expl, prevalence[prevalence$RESPONSE %in% unique(expl$Response),],
              by.x = "Response", by.y = "RESPONSE")

# save(expl, file = paste0(filepath_results, "orthoptera_figures_expl.RData"))
# load(paste0(filepath_results, "orthoptera_figures_expl.RData"))

# PCA over traits
trait_names <- c("Gesamtlaenge", "PronotLaenge", "PronotBreite",
                 "AugenDurchm", "AugenAbs", "Ovipositor", "Femurhinten",
                 "Tibiahinten", "Femurmitte", "Femurvorne")
trait_names <- c("Gesamtlaenge", "PronotLaenge", "PronotBreite",
                 "AugenDurchm", "AugenAbs", "Femurhinten",
                 "Tibiahinten", "Femurmitte", "Femurvorne")
pca <- prcomp(expl[, trait_names], center = TRUE, scale = FALSE)

biplot(pca,  choices = c(1,2))


expl_pca <- cbind(expl, pca$x)
# save(expl_pca, file = paste0(filepath_results, "orthoptera_figures_expl_pca.RData"))

# Correlation
values <- c("KAPPA_MEAN", "POD_MEAN", "FAR_MEAN", "POFD_MEAN",
            "ACCURACY_MEAN", "SR_MEAN", "TS_MEAN", "ETS_MEAN", 
            "HK_MEAN",
            "Gesamtlaenge", "PronotLaenge", "PronotBreite",
            "AugenDurchm", "AugenAbs", "Ovipositor", "Femurhinten",
            "Tibiahinten", "Femurmitte", "Femurvorne",
            "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", 
            "PC8", "PC9", "PC10")

# expl_cor <- cor(expl_pca[, values])
# corrplot(expl_cor, type ="lower")

# Visual analysis
colors <- clrs_hcl(nrow(expl_pca))

rq <- rq(expl_pca[, "Kappa_mean"] ~ expl_pca[, "PC1"], tau = 0.5)

png("figures/fig_05_kappa_pc3.png", 
    width = 1024 * 6, 
    height = 748 * 6, 
    units = "px", 
    res = 600)
xyplot(expl_pca[, "Kappa_mean"] ~ expl_pca[, "PC1"], data = expl_pca,
       fill.color = colors, xlab = "PC #1", ylab = "Kappa mean",
       panel = function(x, y,fill.color,...) {
         fill = fill.color
         panel.xyplot(x, y, pch=19, cex = 2, col=fill, grid = TRUE)
         panel.text(x, y-0.02, cex = 0.75, labels = expl_pca$Response)
         panel.lmline(x, y, ...)
         panel.abline(rq$coefficients, col = "red")
       },
       par.settings = custom.theme(region = colors, alpha = 1.0),
       legend = list(
         top = list(fun = "draw.colorkey", 
                    args = list(
                      list(space = "top",
                           at = seq(min(expl_pca$OCCURENCE),
                                    max(expl_pca$OCCURENCE), 
                                    length.out = nrow(expl_pca)),
                           width = 1, col = colors)))))

graphics.off()

  expl_pca_rq <- rq(expl_pca[, "Kappa_mean"] ~ expl_pca[, "PC1"], 
                    tau = seq(0.05, 0.95, 0.05))
  plot(expl_pca_rq)
  plot(summary(expl_pca_rq), ylim = c(-0.5, 0.5))

  expl_pca_rq <- rq(expl_pca[, "Kappa_mean"] ~ expl_pca[, "PC1"], 
                    tau = seq(0.05, 0.95, 0.05))
  summary(expl_pca_rq)

  
  