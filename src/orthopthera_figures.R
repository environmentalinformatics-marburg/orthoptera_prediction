# Libraries --------------------------------------------------------------------
# devtools::install_github("environmentalinformatics-marburg/gpm")
library(gpm)
library(ggplot2)
library(ggmap)
library(mapview)
library(rgdal)
library(raster)
library(sp)
library(reshape2)
library(quantreg)
library(latticeExtra)


# Set path ---------------------------------------------------------------------
if(Sys.info()["sysname"] == "Windows"){
  filepath_base <- "D:/orthoptera/data/"
} else {
  filepath_base <- "/media/tnauss/myWork/analysis/orthoptera/data/"
}

filepath_results <- paste0(filepath_base, "rdata/")
filepath_obsv <- paste0(filepath_base, "orthoptera/")
filepath_traits <- paste0(filepath_base, "traits/")
filepath_figures <- paste0(filepath_base, "figures/")


# Read data --------------------------------------------------------------------
load(paste0(filepath_results, "orthoptera_prediction_prevalence.RData"))
load(paste0(filepath_results, "orthoptera_figures_trait_subsets_01.RData"))
load(paste0(filepath_results, "orthoptera_figures_trait_subsets_02.RData"))

models_sub <- append(models_sub_01, models_sub_02)
tstat_sub <- append(tstat_sub_01, tstat_sub_02)
tstat_mean_sub <- append(tstat_mean_sub_01, tstat_mean_sub_02)
tstat_merged_sub <- rbind(tstat_merged_sub_01, tstat_merged_sub_02)

rm(models_sub_01, models_sub_02, tstat_sub_01, tstat_sub_02,
   tstat_mean_sub_01, tstat_mean_sub_02, 
   tstat_merged_sub_01, tstat_merged_sub_02)

dem <- raster(paste0(filepath_obsv, "DEM_UTM37S_WGS84_30m_Hemp.tif"))
traits <- read.table(paste0(filepath_traits, "Orthoptera_10_11_2015.csv"),
                     header = TRUE, sep = ";", dec = ".")
band_names <- read.table(paste0(filepath_obsv, "band_names.csv"), header = TRUE,
                         sep = ";")

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
                  data = obsv_samples)

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
                  data = obsv_samples)

png(paste0(filepath_figures, "fig_a01_maps.png"), 
    width = 1024 * 6, 
    height = 748 * 6, 
    units = "px", 
    res = 600)
map + obs + colourscale + labels_map_01 + style_map_01 + theme_bw() + facet_wrap(~ year)
graphics.off()


# Elevation statistics ---------------------------------------------------------
obsv_samples_utm <- obsv_samples
coordinates(obsv_samples_utm) <- ~lon+lat
projection(obsv_samples_utm) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs ")
obsv_samples_utm <- spTransform(obsv_samples_utm, projection(dem))

obsv_samples_utm@data$masl <- extract(dem, obsv_samples_utm)

summary(obsv_samples_utm@data$masl)
quantile(obsv_samples_utm@data$masl, probs = seq(0, 1, 0.1), na.rm = TRUE)


# NDVI statistics --------------------------------------------------------------
obsv_samples_utm@data$NDVI <- 
  (obsv_samples_utm@data$MYD09GA_sur_refl_b02_1 - obsv_samples_utm@data$MYD09GA_sur_refl_b01_1) / 
  (obsv_samples_utm@data$MYD09GA_sur_refl_b02_1 + obsv_samples_utm@data$MYD09GA_sur_refl_b01_1)

summary(obsv_samples_utm@data$NDVI)
quantile(obsv_samples_utm@data$NDVI, probs = seq(0, 1, 0.01), na.rm = TRUE)

quantile(obsv_samples_utm@data$diff_days_nocloud, seq(0, 1, 0.05))


# Prediction performance -------------------------------------------------------
tstat_merged_sub_boxplot <- tstat_merged_sub[grepl("Kappa_mean", tstat_merged_sub$variable),]
tstat_merged_sub_boxplot$variable <- factor(tstat_merged_sub_boxplot$variable,
                                            levels = as.character(unique(tstat_merged_sub_boxplot$variable))[order(as.character(unique(tstat_merged_sub_boxplot$variable)), decreasing = FALSE)])

ggplot(data = tstat_merged_sub_boxplot[grepl("Kappa_mean", tstat_merged_sub_boxplot$variable),], 
       aes(x = variable, y = value, fill = variable)) + 
  geom_boxplot(notch = TRUE) + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0), legend.position = "none")

perf <- do.call("rbind", lapply(unique(as.character(tstat_merged_sub_boxplot$variable)), function(x){
  data.frame(var = x,
             med = median(tstat_merged_sub_boxplot$value[as.character(tstat_merged_sub_boxplot$variable) == x]))
}))
perf[order(perf$med, decreasing = TRUE),]  


levels <- tstat_merged_sub$Response[
  tstat_merged_sub$variable == "Kappa_mean_mspc"][
    order(tstat_merged_sub$value[tstat_merged_sub$variable == "Kappa_mean_mspc"],
          decreasing = FALSE)]
tstat_merged_sub$Response <- factor(tstat_merged_sub$Response, levels = levels)

# ggplot(data = tstat_merged_sub[grepl("Kappa_mean", tstat_merged_sub$variable),], 
#        aes(x = Response, y = value, color = variable, group = variable)) + 
#   geom_line() + 
#   coord_flip()

testvar <- c("Kappa_mean_mspc", "Kappa_mean_lspc_lspt_mspc", "Kappa_mean_mspc_lspt_asl")
plot_top3 <- ggplot(data = tstat_merged_sub[tstat_merged_sub$variable %in% testvar, ], 
                    aes(x = Response, y = value, color = variable, group = variable)) + 
  geom_line() + 
  coord_flip()


plot_scors <- lapply(tstat_sub, function(x){
  plotClassPerformance(x, scores = c("Kappa", "ETS"))  
})


plot_scors <- lapply(tstat_sub, function(x){
  plotClassPerformance(x, scores = c("Kappa"))  
})




modstat <- tstat_sub$tstat_mspc
scores_mean <- paste0(scores, "_mean")
df_mean <- modstat[[1]]
df_mean$Response <- 
  factor(df_mean$Response, levels = df_mean$Response[order(df_mean$Kappa_mean, decreasing = FALSE)])
df_mean_melt <- melt(df_mean, id.vars = "Response")

ggplot(data = df_mean_melt[df_mean_melt$variable %in% scores_mean, ], 
       aes(x = Response, y = value, color = variable, group = variable)) + 
  geom_point() + 
  geom_line(colour = "#1f78b4", size = 1.5) + 
  geom_hline(yintercept = c(0.20, 0.40), colour="#b2df8a", linetype = "longdash") +
  theme_bw() +
  theme(plot.title = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15, face = "italic")) +
  labs(x = "Species", y = "Kappa") +
  coord_flip() +
  theme(legend.position="none")

df <- do.call("rbind", modstat[[2]])
df$Response <- factor(df$Response, levels(df$Response)[order(df_mean$Kappa_mean, 
                                                             decreasing = FALSE)])
df_melt <- melt(df, id.vars = "Response")
indv <- nrow(df[df_melt$Response == unique(df_melt$Response)[[1]] & 
                  df_melt$variable == unique(df_melt$variable)[[1]], 
                ])
df_melt$variable_id <- paste0(df_melt$variable, "_", 
                              seq(indv))

ggplot(data = df_melt[df_melt$variable %in% scores, ], 
             aes(x = Response, y = value, colour = variable, group = variable_id)) + 
  geom_line(linetype = "twodash", alpha = 0.25, colour = "#a6cee3") + 
  geom_point(data = df_mean_melt[df_mean_melt$variable %in% scores_mean, ], 
             aes(x = Response, y = value, colour = variable, group = variable)) + 
  geom_line(data = df_mean_melt[df_mean_melt$variable %in% scores_mean, ], 
            aes(x = Response, y = value, colour = variable, group = variable), 
            colour = "#1f78b4", size = 1.5) + 
  geom_hline(yintercept = c(0.20, 0.40), colour="#b2df8a", linetype = "longdash") +
  theme_bw() +
  theme(plot.title = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  labs(x = "Species", y = "Kappa") +
  coord_flip() +
  theme(legend.position="none")





plot_scors$tstat_mspc
plot_scors$tstat_lspc_lspt_mspc
plot_scors$tstat_mspc_lspt_asl


png(paste0(filepath_figures, "fig_02_scores.png"), 
    width = 1024 * 6, 
    height = 748 * 6, 
    units = "px", 
    res = 600)
plot_top3
graphics.off()


png(paste0(filepath_figures, "fig_03_scores.png"), 
    width = 1024 * 6, 
    height = 748 * 6, 
    units = "px", 
    res = 600)
plot_scors$tstat_mspc
graphics.off()


# Variable importance ----------------------------------------------------------
plotVarImpHeatmap <- function(var_imp, xlab = "Variable", ylab = "Method",
                              vis_range = "minmax"){
  temp <- do.call("rbind", var_imp)
  temp$VARIABLE <- factor(temp$VARIABLE, 
                          levels = sort(as.character(unique(temp$VARIABLE))))
  if(vis_range == "minmax"){
    vis_range <- c(min(temp$mean), max(temp$mean))
  }
  clr <- colorRampPalette(brewer.pal(9, "YlOrRd"))
  lattice::levelplot(mean ~ RESPONSE * VARIABLE, data = temp,
                     col.regions = clr(101), at = seq(vis_range[1], vis_range[2], length.out = 101),
                     asp = 1, as.table = TRUE,
                     ylab = ylab, xlab = xlab,
                     scales = list(x = list(rot = 45, font = "italic")),
                     main = "Variable importance",
                     cex.title = 1,
                     colorkey = list(space = "top",
                                     width = 1, height = 0.75),
                     panel=function(...) {
                       grid::grid.rect(gp=grid::gpar(col=NA, fill="grey60"))
                       panel.levelplot(...)
                     })
  
}

var_imp <- lapply(models_sub, function(x){
  compVarImp(x, scale = FALSE)
})

var_imp_scale <- lapply(models_sub, function(x){
  compVarImp(x, scale = TRUE)
})

var_imp_plot <- lapply(var_imp, function(x){
  plotVarImp(x)
})


var_imp_scale_new <- lapply(var_imp_scale, function(x){
  lapply(x, function(y){
    y$RESPONSE <- as.factor(gsub("\\.", " ", y$RESPONSE))
    y <- merge(y, band_names, by.x = "VARIABLE", by.y = "Original")
    y$VARIABLE <- y$New
    y$New <- NULL
    y$VARIABLE <- factor(y$VARIABLE,levels(y$VARIABLE[order(y$VARIABLE, decreasing = TRUE)]))
    return(y)
  })
})

var_imp_heat <- lapply(var_imp_scale_new, function(x){
  plotVarImpHeatmap(x, xlab = NULL, ylab = "Band")
})



png(paste0(filepath_figures, "fig_05_variable_importance_01.png"), 
    width = 1024 * 6, 
    height = 748 * 6, 
    units = "px", 
    res = 600)
var_imp_heat$mspc
graphics.off()

png(paste0(filepath_figures, "fig_05_variable_importance_02.png"), 
    width = 1024 * 6, 
    height = 748 * 6, 
    units = "px", 
    res = 600)
var_imp_heat$lspc_lspt_mspc
graphics.off()

png(paste0(filepath_figures, "fig_05_variable_importance_03.png"), 
    width = 1024 * 6, 
    height = 748 * 6, 
    units = "px", 
    res = 600)
var_imp_heat$mspc_lspt_asl
graphics.off()


# Test dependency of prediction accuracy ---------------------------------------
obsv_samples_utm@data$pid <- paste0(obsv_samples_utm@data$plot, "_", obsv_samples_utm@data$date_observation)

test <- lapply(models_sub$mspc, function(x){
  act_test <- lapply(x, function(y){
    if(inherits(y$model, "try-error")){
      NULL
    } else {
      pid <- data.frame(pid = paste0(y$testing$META$plot, "_", as.Date(y$testing$META$date_observation, format="%Y-%j")))
      ndvi = merge(pid, obsv_samples_utm@data, by = "pid")$NDVI
      df <- data.frame(TARGET = y$response,
                       RESPONSE = y$testing$RESPONSE,
                       PREDICTION = y$testing$PREDICTED$pred,
                       NDVI = ndvi,
                       y$testing$META)
      return(df)
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

test_mean <- merge(tstat_sub$tstat_mspc[[1]], prevalence, by.x = "Response", by.y = "RESPONSE")
test_mean <- merge(test_mean, ndvi_mean, by = "Response")
test_mean <- merge(test_mean, asl_mean, by = "Response")
test_mean <- merge(test_mean, rainfall_mean, by = "Response")
test_mean <- merge(test_mean, species_mean, by = "Response")

traits[, 5:15] <- log(traits[, 5:15])
traits_mean <- aggregate(traits[, 5:15], by = list(traits$Name), FUN = mean,
                         na.rm = TRUE)
test_mean$Response <- gsub("\\.", " ", test_mean$Response) 
test_mean <- merge(test_mean, traits_mean, by.x = "Response", by.y = "Group.1")

nv <- function(x){
  (x - min(x)) * (1 - 0) / (max(x) - min(x)) + 0
}

test_mean$OCCURENCE_scaled <- nv(test_mean$OCCURENCE)
test_mean$NDVI_scaled <- nv(test_mean$NDVI)
test_mean$asl_scaled <- nv(test_mean$asl)
test_mean$rainfall_scaled <- nv(test_mean$rainfall)
test_mean$nr.of.species_scaled <- nv(test_mean$nr.of.species)

test_mean_unscaled_melt <- melt(test_mean[, colnames(test_mean) %in%
                                          c("Kappa_mean", "OCCURENCE", 
                                            "nr.of.species_scaled",
                                            "Gesamtlaenge")], 
                              id.vars = "Kappa_mean")

test_mean_scaled_melt <- melt(test_mean[, colnames(test_mean) %in%
                                   c("Kappa_mean", "OCCURENCE_scaled", 
                                     "asl_scaled", 
                                     "NDVI_scaled",
                                     "rainfall_scaled", 
                                     "nr.of.species_scaled")], 
                       id.vars = "Kappa_mean")

ggplot(data = test_mean_scaled_melt, aes(x = value, y = Kappa_mean)) + 
  geom_point() + 
  geom_smooth() + 
  facet_wrap(~variable,  scales = "free_y")


lm_eqn <- function(df){
  sm <- summary(lm(Kappa_mean ~ value, df));
  eq <- substitute(italic(r)^2~"="~r2*","~~italic(p)~"="~pv, 
                   list(r2 = format(sm$r.squared, digits = 3),
                        pv = format(pf(sm$fstatistic[1L], sm$fstatistic[2L], sm$fstatistic[3L],lower.tail = FALSE), digits = 2)))
  as.character(as.expression(eq));                 
}

lm_eqn_occurence <- lm_eqn(test_mean_unscaled_melt[as.character(test_mean_unscaled_melt$variable) == "OCCURENCE",])

plot_occurence <- ggplot(data = test_mean_unscaled_melt[as.character(test_mean_unscaled_melt$variable) == "OCCURENCE",], 
                         aes(x = value, y = Kappa_mean)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) + 
  theme_bw() +
  # theme(axis.title = element_text(size = 30), axis.text = element_text(size = 22)) + 
  # facet_wrap(~variable,  scales = "free") + 
  geom_text(x = 30, y = 0.42, label = lm_eqn_occurence, parse = TRUE) + 
  xlab("Occurence of the species") + ylab("Kappa (mean)")





lm_eqn_body_length <- lm_eqn(test_mean_unscaled_melt[as.character(test_mean_unscaled_melt$variable) == "Gesamtlaenge",])

plot_body_length <- ggplot(data = test_mean_unscaled_melt[as.character(test_mean_unscaled_melt$variable) == "Gesamtlaenge",], 
                         aes(x = value, y = Kappa_mean)) + 
  geom_point(size = 4) + 
  geom_smooth(method = "lm", se = FALSE) + 
  theme_bw() +
  theme(axis.title = element_text(size = 30), axis.text = element_text(size = 22)) + 
  # facet_wrap(~variable,  scales = "free") + 
  geom_text(x = 1.15, y = 0.30, label = lm_eqn_body_length, parse = TRUE, size = 6) + 
  xlab("Body length of the species") + ylab("Kappa (mean)")



png(paste0(filepath_figures, "fig_04a_scores_dependency.png"), 
    width = 748 * 6, 
    height = 748 * 6, 
    units = "px", 
    res = 600)
plot_occurence
graphics.off()

png(paste0(filepath_figures, "fig_04b_scores_dependency.png"), 
    width = 748 * 6, 
    height = 748 * 6, 
    units = "px", 
    res = 600)
plot_body_length
graphics.off()

# ggplot(data = test_mean_unscaled_melt[as.character(test_mean_unscaled_melt$variable) == "Gesamtlaenge",], 
#        aes(x = value, y = Kappa_mean)) + 
#   geom_point(size = 10) + 
#   geom_smooth(method = lm) + 
#   # facet_wrap(~variable,  scales = "free") + 
#   geom_text(x = 1.3, y = 0.40, label = lm_eqn_body_length, parse = TRUE, size = 12) + 
#   xlab("Body length of the species") + ylab("Kappa (mean)") + 
#   theme(axis.title = element_text(size = 60), axis.text = element_text(size = 44))

# ggplot(data = test_mean, aes(x = OCCURENCE, y = Kappa_mean)) + geom_point() + geom_smooth()
# ggplot(data = test_mean, aes(x = NDVI, y = Kappa_mean)) + geom_point() + geom_smooth()
# ggplot(data = test_mean, aes(x = asl, y = Kappa_mean)) + geom_point() + geom_smooth()
# ggplot(data = test_mean, aes(x = rainfall, y = Kappa_mean)) + geom_point() + geom_smooth()
# ggplot(data = test_mean, aes(x = nr.of.species,  y= Kappa_mean)) + geom_point() + geom_smooth()

test_utm_melt <- 
  melt(test_utm@data[, colnames(test_utm@data) %in% 
                       c("TARGET", "performance", "asl", "NDVI", "rainfall", "nr.of.species")],
       id.var = c("TARGET", "performance"))

png(paste0(filepath_figures, "fig_06_scores_dependency_boxplots.png"), 
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

tstat_best <- tstat_mean_sub$tstat_mspc_mean
# Adjust test statistics
tstat_best$Response <- gsub("\\.", " ", tstat_best$Response) 
tstat_best$Response[tstat_best$Response == "Cyrtacanthacris tatarica"] <-
  "Cyrtacanthacris tatarica tatarica"
tstat_best$Response[tstat_best$Response == "Gymnobothrus flexuosus"] <-
  "Gymnobothrus temporalis flexuosus"
tstat_best$Response_SHORT <- substr(tstat_best$Response, 1, 4)

# Adjust prevalence information
prevalence$RESPONSE <- gsub("\\.", " ", prevalence$RESPONSE) 
prevalence$RESPONSE[prevalence$RESPONSE == "Cyrtacanthacris tatarica"] <-
  "Cyrtacanthacris tatarica tatarica"
prevalence$RESPONSE[prevalence$RESPONSE == "Gymnobothrus flexuosus"] <-
  "Gymnobothrus temporalis flexuosus"

traits <- traits[traits$Name %in% tstat_best$Response, ]
traits[, 5:15] <- log(traits[, 5:15])
traits_mean <- aggregate(traits[, 5:15], by = list(traits$Name), FUN = mean,
                         na.rm = TRUE)
# Combine model statistics and traits and analyse information
expl <- cbind(tstat_best[tstat_best$Response %in% unique(traits_mean$Group.1),],
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


png(paste0(filepath_figures, "fig_07_pca_biplot.png"), 
    width = 1024 * 6, 
    height = 748 * 6, 
    units = "px", 
    res = 600)
biplot(pca,  choices = c(1,2))
graphics.off()


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

rq <- rq(expl_pca[, "Kappa_mean"] ~ expl_pca[, "PC1"], tau = 0.75)

png(paste0(filepath_figures, "fig_08_variable_importance_01.png"), 
    width = 1024 * 6, 
    height = 748 * 6, 
    units = "px", 
    res = 600)
xyplot(expl_pca[, "Kappa_mean"] ~ expl_pca[, "Gesamtlaenge"], data = expl_pca,
       fill.color = colors, xlab = list("Overall length",  cex=2), 
       ylab = list("Kappa mean",  cex=2), 
       panel = function(x, y,fill.color,...) {
         fill = fill.color
         panel.xyplot(x, y, pch=19, cex = 2, col=fill, grid = TRUE)
         panel.text(x, y-0.02, cex = 0.75, labels = expl_pca$Response)
         panel.lmline(x, y, ...)
         # panel.abline(rq$coefficients, col = "red")
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

summary(lm(expl_pca[, "Kappa_mean"] ~ expl_pca[, "Gesamtlaenge"]))


expl_pca_rq <- rq(expl_pca[, "Kappa_mean"] ~ expl_pca[, "PC1"], 
                  tau = seq(0.05, 0.95, 0.05))
plot(expl_pca_rq)
plot(summary(expl_pca_rq), ylim = c(-0.5, 0.5))

expl_pca_rq <- rq(expl_pca[, "Kappa_mean"] ~ expl_pca[, "PC1"], 
                  tau = seq(0.05, 0.95, 0.05))
summary(expl_pca_rq)


