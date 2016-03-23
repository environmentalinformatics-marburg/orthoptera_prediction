setwd("D:/active/orthoptera/data")

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

# Read data --------------------------------------------------------------------
load("processed/orthoptera.rda")
load("processed/prevalent_species.rda")
load("processed/orthoptera_trte.rda")
load("processed/prevalence.rda")
load("processed/orthoptera_trte.rda")
load("processed/tstat_mean_rf.rda")
load("processed/var_imp.rda")
load("processed/models_rf-2016-03-04.rda")

dem <- raster("original/DEM_UTM37S_WGS84_30m_Hemp.tif")


dfx <- lapply(orthoptera_trte, function(x){
  dfy <- lapply(x, function(y){
    c(y$training$SAMPLES, y$testing$SAMPLES)
  })
  unique(unlist(dfy))
})
rows <- unique(unlist(dfx))

sample <- orthoptera@data$input[rows, ]

sample$date_observation <- as.Date(sample$date_observation, format="%Y-%j")
sample$year <- as.numeric(format(sample$date_observation, "%Y"))


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

png("figures/fig_01_map.png", 
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

png("figures/fig_a01_maps.png", 
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
  (sample_utm@data$greyval_band_02 - sample_utm@data$greyval_band_01) / 
  (sample_utm@data$greyval_band_02 + sample_utm@data$greyval_band_01)

summary(sample_utm@data$NDVI)
quantile(sample_utm@data$NDVI, probs = seq(0, 1, 0.1), na.rm = TRUE)

gsub("\\.", " ", prevalent_species) 

quantile(sample_utm@data$diff_days_nocloud, seq(0, 1, 0.05))


# Prediction performance -------------------------------------------------------
plot_scores <- plotClassPerformance(tstat, scores = c("Kappa", "ETS"))

png("figures/fig_02_scores.png", 
    width = 1024 * 6, 
    height = 748 * 6, 
    units = "px", 
    res = 600)
plot_scores
graphics.off()


# Test dependency of prediction accuracy ---------------------------------------
test <- lapply(models, function(x){
  act_test <- lapply(x, function(y){
    df <- data.frame(TARGET = y$response,
                     RESPONSE = y$testing$RESPONSE,
                     PREDICTION = y$testing$PREDICTED$pred,
                     NDVI = (y$testing$INDEPENDENT$greyval_band_02 - y$testing$INDEPENDENT$greyval_band_01) / 
                       (y$testing$INDEPENDENT$greyval_band_02 + y$testing$INDEPENDENT$greyval_band_01),
                     y$testing$META)
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

png("figures/fig_03_scores_dependency.png", 
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

traits <- read.table("traits/Orthoptera_10_11_2015.csv",
                     header = TRUE, sep = ";", dec = ".")
traits <- traits[traits$Name %in% tstat_mean$Response, ]
traits[, 5:15] <- log(traits[, 5:15])
traits_mean <- aggregate(traits[, 5:15], by = list(traits$Name), FUN = mean,
                         na.rm = TRUE)
# Combine model statistics and traits and analyse information
expl <- cbind(tstat_mean, traits_mean, by.x = "Response", by.y = "Group.1")
expl <- cbind(expl, prevalence, by.x = "Response", by.y = "RESPONSE")

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
# save(expl_pca, file = "D:/active/orthoptera/orthoptera_prediction/src/expl_pca.rda")

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

rq <- rq(expl_pca[, "Kappa_mean"] ~ expl_pca[, "PC2"], tau = 0.3)

png("figures/fig_05_kappa_pc3.png", 
    width = 1024 * 6, 
    height = 748 * 6, 
    units = "px", 
    res = 600)
xyplot(expl_pca[, "Kappa_mean"] ~ expl_pca[, "PC2"], data = expl_pca,
       fill.color = colors, xlab = "PC #3", ylab = "Kappa mean",
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

  expl_pca_rq <- rq(expl_pca[, "Kappa_mean"] ~ expl_pca[, "PC2"], 
                    tau = seq(0.05, 0.95, 0.05))
  plot(expl_pca_rq)
  plot(summary(expl_pca_rq), ylim = c(-0.5, 0.5))

  expl_pca_rq <- rq(expl_pca[, "Kappa_mean"] ~ expl_pca[, "PC2"], 
                    tau = seq(0.05, 0.95, 0.05))
  summary(expl_pca_rq)

  
  