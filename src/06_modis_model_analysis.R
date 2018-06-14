# Predict orthoptera based on satellite observations
# Thomas Nauss

if(Sys.info()["sysname"] == "Windows"){
  source("D:/orthoptera/orthoptera_prediction/src/00_set_environment.R")
} else {
  source("/media/tnauss/myWork/analysis/orthoptera/orthoptera_prediction/src/00_set_environment.R")
}

compute <- TRUE

# Read modis model -------------------------------------------------------------
obsv_myd <- readRDS(file = paste0(path_results, "obsv_myd_modis_trainmodel.rds"))
obsv_mod <- readRDS(file = paste0(path_results, "obsv_mod_modis_trainmodel.rds"))

obsv_modis <- obsv_mod

tstat_modis <- compContTests(obsv_modis@model$rf_rfe, mean = TRUE)
tstat_modis_mean <- tstat_modis[[1]]
tstat_modis_indv <- tstat_modis[[2]]

# Merge with trait information -------------------------------------------------
obsv_traits <- readRDS(paste0(path_results, "obsv_traits.rds"))
tstat_modis_mean_traits <- merge(tstat_modis_mean, obsv_traits, by.x = "Response", by.y = "Name")
tstat_modis_indv_traits <- merge(tstat_modis_indv, obsv_traits, by.x = "Response", by.y = "Name")

tstat_modis_mean$Response <- factor(tstat_modis_mean$Response, tstat_modis_mean$Response[order(tstat_modis_mean$Kappa_mean)])
tstat_modis_indv$Response <- factor(tstat_modis_indv$Response, levels(tstat_modis_mean$Response))

tstat_modis_mean_traits$Response <- as.character(tstat_modis_mean_traits$Response)
tstat_modis_mean_traits$Response <- factor(tstat_modis_mean_traits$Response, tstat_modis_mean_traits$Response[order(tstat_modis_mean_traits$Kappa_mean)])

tstat_modis_indv_traits$Response <- as.character(tstat_modis_indv_traits$Response)
tstat_modis_indv_traits$Response <- factor(tstat_modis_indv_traits$Response, levels(tstat_modis_mean_traits$Response))


ggplot(data = tstat_modis_indv, aes(x = Response, y = Kappa)) + 
  geom_boxplot() + 
  geom_line(data = tstat_modis_mean, aes(x = Response, y = Kappa_mean, group = 1, color = "red")) +
  # geom_line(colour = "#1f78b4", size = 1.5) + 
  geom_hline(yintercept = c(0.20, 0.40), colour="#b2df8a", linetype = "longdash") +
  theme_bw() +
  theme(plot.title = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15, face = "italic")) +
  labs(x = "Species", y = "Kappa") +
  coord_flip() +
  theme(legend.position="none")


ggplot(data = tstat_modis_indv_traits, aes(x = Response, y = Kappa)) + 
  geom_boxplot() + 
  geom_line(data = tstat_modis_mean_traits, aes(x = Response, y = Kappa_mean, group = 1, colour = "red")) +
  geom_line(data = tstat_modis_mean_traits, aes(x = Response, y = HK_mean, group = 1, colour = "blue")) +
  # geom_line(colour = "#1f78b4", size = 1.5) + 
  geom_hline(yintercept = c(0.20, 0.40), colour="#b2df8a", linetype = "longdash") +
  theme_bw() +
  theme(plot.title = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15, face = "italic")) +
  labs(x = "Species", y = "Kappa") +
  coord_flip() +
  theme(legend.position="none")


ggplot(data = tstat_modis_mean_traits, aes(x = PC1, y = HK_mean)) + 
  geom_point() + 
  geom_hline(yintercept = c(0.20, 0.40), colour="#b2df8a", linetype = "longdash") +
  geom_smooth(method = "lm") + 
  theme_bw() +
  theme(plot.title = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15, face = "italic")) +
  labs(x = "Trait component 1", y = "ETS") +
  theme(legend.position="none")

summary(lm(PC1 ~ Kappa_mean, data = tstat_modis_mean_traits))
summary(lm(PC1 ~ HK_mean, data = tstat_modis_mean_traits))



var_imp <- compVarImp(obsv_modis@model$rf_rfe, scale = FALSE)

var_imp_scale <- compVarImp(obsv_modis@model$rf_rfe, scale = TRUE)

var_imp_plot <- plotVarImp(var_imp)

var_imp_heat <- plotVarImpHeatmap(var_imp_scale, xlab = "Species", ylab = "Band")


tstat_modis_mean <- merge(tstat_modis[[1]], obsv_modis@meta$input$MIN_OCCURENCE,
                        by.x = "Response", by.y="names")

tstat_modis_mean[order(tstat_modis_mean$Kappa_mean, decreasing = TRUE),]

ggplot(data = tstat_modis_mean, aes(x = mo_mean, y = Kappa_mean)) + geom_point() + geom_smooth()





modstat <- tstat_modis_sub$tstat_modis_mspc
scores_mean <- paste0(scores, "_mean")
df_mean <- modstat[[1]]
df_mean$Response <- 
  factor(df_mean$Response, levels = df_mean$Response[order(df_mean$Kappa_mean, decreasing = FALSE)])
df_mean_melt <- melt(df_mean, id.vars = "Response")


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

