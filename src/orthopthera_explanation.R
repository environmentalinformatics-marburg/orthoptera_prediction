setwd("D:/active/orthoptera/data")

# Libraries --------------------------------------------------------------------
# library(gpm)
library(ggplot2)
# library(ggvis)
library(latticeExtra)
library(corrplot)
library(quantreg)

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


# Read model test statistics and prevalence information ------------------------
load("processed/prevalence.rda")
load("processed/prevalent_species.rda")
load("processed/tstat_mean_rf.rda")

# Adjust test statistics
tstat_mean$RESPONSE <- gsub("\\.", " ", tstat_mean$RESPONSE) 
tstat_mean$RESPONSE[tstat_mean$RESPONSE == "Cyrtacanthacris tatarica"] <-
  "Cyrtacanthacris tatarica tatarica"
tstat_mean$RESPONSE[tstat_mean$RESPONSE == "Gymnobothrus flexuosus"] <-
  "Gymnobothrus temporalis flexuosus"
tstat_mean$RESPONSE_SHORT <- substr(tstat_mean$RESPONSE, 1, 4)

# Adjust prevalence information
prevalence$RESPONSE <- gsub("\\.", " ", prevalence$RESPONSE) 
prevalence$RESPONSE[prevalence$RESPONSE == "Cyrtacanthacris tatarica"] <-
  "Cyrtacanthacris tatarica tatarica"
prevalence$RESPONSE[prevalence$RESPONSE == "Gymnobothrus flexuosus"] <-
  "Gymnobothrus temporalis flexuosus"


# Read traits information ------------------------------------------------------
traits <- read.table("traits/Orthoptera_10_11_2015.csv",
                     header = TRUE, sep = ";", dec = ".")
traits <- traits[traits$Name %in% tstat_mean$RESPONSE, ]
traits[, 5:15] <- log(traits[, 5:15])
traits_mean <- aggregate(traits[, 5:15], by = list(traits$Name), FUN = mean,
                           na.rm = TRUE)


# Combine model statistics and traits and analyse information ------------------
expl <- cbind(tstat_mean, traits_mean, by.x = "RESPONSE", by.y = "Group.1")
expl <- cbind(expl, prevalence, by = "RESPONSE")

# PCA over traits
trait_names <- c("Gesamtlaenge", "PronotLaenge", "PronotBreite",
                 "AugenDurchm", "AugenAbs", "Ovipositor", "Femurhinten",
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

expl_cor <- cor(expl_pca[, values])

# corrplot(expl_cor, type ="lower")



# Visual analysis
expl_pca_lm <- lm(ETS_MEAN ~ Gesamtlaenge, data = expl_pca)
summary(expl_pca_lm)

expl_pca_rq <- rq(ETS_MEAN ~ PC1, data = expl_pca, tau = seq(0.01, 0.90, 0.10))
expl_pca_rq <- rq(ETS_MEAN ~ Gesamtlaenge, data = expl_pca, tau = seq(0.05, 0.95, 0.05))
summary(expl_pca_rq)
plot(expl_pca_rq)
plot(summary(expl_pca_rq), parm = 2)
plot(expl_pca_rq, parm=2)

png("expl_pca_rq_gesamtlaenge_slope.png")
plot(summary(expl_pca_rq), parm = 2, main = "Total length",
     cex = 1, pch = 19, xlim = c(-1, 1), ylim = c(-0.5, 0.5))
dev.off()

colors <- clrs_hcl(nrow(expl_pca))

ggplot(expl_pca, aes(x = FAR_MEAN, y = POD_MEAN, group = RESPONSE)) +
  geom_point(aes(colour = OCCURENCE), size = 5) +   
  geom_text(aes(label = substr(RESPONSE, 1, 5
                               )), size = 3, hjust=0.5, vjust=2) + 
  ggtitle("POD vs FAR") + 
  labs(x = "FAR (mean)", y = "POD (mean)", colour = "Occurence") + 
  scale_colour_gradientn(colours=colors)


png("expl_pca_ets_mean-Gesamtlaenge.png")
ggplot(expl_pca, aes(x = Gesamtlaenge, y = ETS_MEAN, group = RESPONSE)) +
  geom_point(aes(colour = OCCURENCE), size = 5) +   
  geom_text(aes(label = substr(RESPONSE, 1, 5
  )), size = 3, hjust=0.5, vjust=2) + 
  ggtitle("Orthoptera") + 
  labs(x = "Total length", y = "ETS mean", colour = "Occurence") + 
  scale_colour_gradientn(colours=colors)
dev.off()

png("expl_pca_pod-far.png")
ggplot(expl_pca, aes(x = FAR_MEAN, y = POD_MEAN, group = RESPONSE)) +
  geom_point(aes(colour = OCCURENCE), size = 5) +   
  geom_text(aes(label = substr(RESPONSE, 1, 5
  )), size = 3, hjust=0.5, vjust=2) + 
  ggtitle("Orthoptera") + 
  labs(x = "FAR mean", y = "POD mean", colour = "Occurence") + 
  scale_colour_gradientn(colours=colors)
dev.off()





# xyplot(ETS_MEAN ~ Gesamtlaenge, data = expl_pca,
#        fill.color = colors, xlab = "Total length (cm)", ylab = "ETS",
#        panel = function(x, y,fill.color,...) {
#          fill = fill.color
#          panel.xyplot(x, y, pch=19, cex = 2, col=fill, grid = TRUE)
#          #          panel.abline(expl_pca_lm)
#          #          panel.abline(expl_pca_rq$coefficients[,1], col = "green")
#          #          panel.abline(expl_pca_rq$coefficients[,2], col = "red")
#          #          panel.lmline(x, y, ...)
#          panel.text(x, y-0.02, cex = 0.75, labels = substr(expl_pca$RESPONSE, 1, 5))},
#        par.settings = custom.theme(region = colors, alpha = 1.0),
#        legend = list(top = list(fun = "draw.colorkey", 
#                                 args = list(list(space = "top",
#                                                  at = seq(min(expl_pca$OCCURENCE),max(expl_pca$OCCURENCE), length.out = nrow(expl_pca)),
#                                                  width = 1,
#                                                  col = colors)))))
# 
# 

png("expl_pca_pod-far.png")
xyplot(POD_MEAN ~ FAR_MEAN, z = expl_pca$OCCURENCE, data = expl_pca,
       fill.color = colors, xlab = "FAR mean", ylab = "POD mean",
       xlim = c(0,1), ylim = c(0,1),
       panel = function(x, y, z, fill.color,...) {
         fill = fill.color
         panel.levelplot.points(x = x, y = y, z = z, cex = 2)
         panel.abline(a = 0, b = 1,  lty = 2, col = "grey")
         #          panel.abline(expl_pca_rq$coefficients[,1], col = "green")
         #          panel.abline(expl_pca_rq$coefficients[,2], col = "red")
         #          panel.lmline(x, y, ...)
         panel.text(x, y-0.02, cex = 0.75, labels = substr(expl_pca$RESPONSE, 1, 5))},
       aspect = "iso",
       par.settings = custom.theme(region = colors, alpha = 1.0),
       legend = list(top = list(fun = "draw.colorkey", 
                                args = list(list(space = "top",
                                                 at = seq(min(expl_pca$OCCURENCE),max(expl_pca$OCCURENCE), length.out = nrow(expl_pca)),
                                                 width = 1,
                                                 col = colors)))))
dev.off()


png("expl_pca_ets-gesamtlaenge.png")
xyplot(ETS_MEAN ~ Gesamtlaenge, z = expl_pca$OCCURENCE, data = expl_pca,
       fill.color = colors, xlab = "Total length", ylab = "ETS",
       panel = function(x, y, z, fill.color,...) {
         fill = fill.color
         panel.levelplot.points(x = x, y = y, z = z, cex = 2)
         #          panel.abline(expl_pca_lm)
         panel.abline(expl_pca_rq$coefficients[,18], col = "black")
         #          panel.abline(expl_pca_rq$coefficients[,2], col = "red")
         #          panel.lmline(x, y, ...)
         panel.text(x, y-0.02, cex = 0.75, labels = substr(expl_pca$RESPONSE, 1, 5))},
       par.settings = custom.theme(region = colors, alpha = 1.0),
       legend = list(top = list(fun = "draw.colorkey", 
                                args = list(list(space = "top",
                                                 at = seq(min(expl_pca$OCCURENCE),max(expl_pca$OCCURENCE), length.out = nrow(expl_pca)),
                                                 width = 1,
                                                 col = colors)))))
dev.off()

clr <- colorRampPalette(brewer.pal(9, "YlOrRd"))
levelplot(OCCURENCE ~ factor(RESPONSE) * ETS_MEAN, data = expl_pca,
          col.regions = clr(101), at = seq(0, 100, 1),
          asp = 1, as.table = TRUE,
          ylab = ylab, xlab = xlab,
          scales = list(x = list(rot = 45)),
          main = "Variable importance",
          cex.title = 1,
          colorkey = list(space = "top",
                          width = 1, height = 0.75),
          panel=function(...) {
            grid.rect(gp=gpar(col=NA, fill="grey60"))
            panel.levelplot(...)
          })


scatter3D(x = expl_pca$OCCURENCE, y = expl_pca$ETS_MEAN, z = expl_pca$OCCURENCE)

# GGVIS
# colors <- clrs_hcl(2)
# slider <- input_slider(10, 1000)
# yv <- input_select(values, label = "y axis", map = as.name)
# xv <- input_select(values, label = "x axis", map = as.name)
# expl_cor %>% ggvis(xv, yv) %>%
#   layer_points(fill = ~OCCURENCE, size := slider) %>%
#   layer_text(text := ~RESPONSE_SHORT, fontSize = 3) %>%
#   scale_numeric("fill", range = colors) %>%
#   add_axis("x", title = names(xv))
