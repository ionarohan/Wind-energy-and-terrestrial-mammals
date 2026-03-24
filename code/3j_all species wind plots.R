###########################################################
## Code for: Turbine visibility is a strong predictor of ##
## altered habitat selection by terrestrial mammals at a ##
######## wind energy facility in central New Mexico #######
###########################################################
#### Script for the creation of the plots showing the #####
##### effects on the the wind energy variables on the #####
#  habitat selection of terrestrial mammals in central NM #
###########################################################
################ Script by: Iona Rohan ####################
############ Contact: ionarohan12@gmail.com ###############
###########################################################
########## Date Last Modified: 23-March-2026 ##############
###########################################################

#Clear work environment
rm(list=ls())

#Note: If you opened this script through the .Rproj file, the only line you 
  #should need to change for the script to run (assuming packages are installed)   #is the homewd directory on line 24.

#Set home working directory
  # e.g. homewd = "C:/Users/ionar/Desktop/R Repository/Wind-energy-and-terrestrial-mammals/"
homewd = "<insert your folder here and end with a forward slash>"

#Set wd to data folder on your local computer 
setwd(paste0(homewd, "data/"))

# Load packages 
library(unmarked)
library(ggplot2)
library(tidyverse)
library(scales)
library(patchwork)
library(cowplot)

########################################################
################# MULE DEER PLOTS ######################
########################################################

#### Setup code for mule deer models ####

# Read in edited .csv
detHist <- read.csv(file = "mule deer detection hist.csv", row.names = 1)
# Read in site.covs.scaled
site.covs.scaled <- readRDS("site_covs_scaled.RData")
# Read in obsCovs.scaled
obsCovs.scaled <- readRDS("obsCovs_scaled.RData")

# Change from integer to numeric
detHist <- detHist %>%
  mutate(across(1:51, as.numeric))

# Create detection history 
occu.odhe <- unmarkedFrameOccu(y=detHist, 
                               siteCovs = site.covs.scaled, 
                               obsCovs = obsCovs.scaled)

# Read in unscaled site-level covariates
site.covs <- read.csv("site_covs.csv", nrows = 102, header = TRUE)

#### Turbine visibility * biotic community + shrub density model plot ####

Bio_Com_X_Turbine_Vis <- occu( ~ max_trig_dist + livestock.count ~ 
                                 as.factor(biotic_com_2) * 
                                 X150cm_turbine_vis + shrub_yucca_density, 
                                 occu.odhe, starts = c(0,5,0,2,5,-3,0,-2))

# Convert biotic_com_2 to factor
occu.odhe@siteCovs$biotic_com_2 <- as.factor(occu.odhe@siteCovs$biotic_com_2)
biotic_levels <- levels(occu.odhe@siteCovs$biotic_com_2)

# Fit model using scaled X150cm_turbine_vis
mod1 <- occu( ~ max_trig_dist + livestock.count ~ shrub_yucca_density +
                as.factor(biotic_com_2) * X150cm_turbine_vis, occu.odhe)

# Get mean and SD of X150cm_turbine_vis
wind_mean <- mean(site.covs$X150cm_turbine_vis, na.rm = TRUE)
wind_sd   <- sd(site.covs$X150cm_turbine_vis, na.rm = TRUE)

# Generate sequence of unscaled values, then scale them
turbine_seq_unscaled <- seq(min(site.covs$X150cm_turbine_vis, 
                                na.rm = TRUE),
                            max(site.covs$X150cm_turbine_vis, 
                                na.rm = TRUE),
                            length.out = 100)
turbine_seq_scaled <- (turbine_seq_unscaled - wind_mean) / wind_sd

# Create plot prediction dataset
plot_data <- expand.grid(
  X150cm_turbine_vis = turbine_seq_scaled,
  biotic_com_2 = factor(biotic_levels, levels = c("1", "2"), 
                        labels = c("Grassland", "Woodland")),
  shrub_yucca_density = 0
)

# Predict
preds <- predict(mod1, type = "state", newdata = plot_data)
plot_data$Predicted <- preds$Predicted
plot_data$lower <- preds$lower
plot_data$upper <- preds$upper
plot_data$X150cm_turbine_vis_unscaled <- turbine_seq_unscaled

# Add species label for combined plots 
plot_data$species <- "Mule deer"
deer_df_vis <- plot_data

# Create site_data for observed points
site_data <- site.covs %>%
  dplyr::select(X150cm_turbine_vis, biotic_com_2, shrub_yucca_density) %>%
  mutate(
    X150cm_turbine_vis_scaled = 
      (X150cm_turbine_vis - wind_mean) / wind_sd,
    shrub_yucca_density = 0
  )

# Recode biotic_com_2
site_data$biotic_com_2 <- factor(as.character(site_data$biotic_com_2),
                                 levels = c("1", "2"),
                                 labels = c("Grassland", "Woodland"))

# Predict for observed sites
site_preds <- predict(mod1, type = "state", newdata = data.frame(
  X150cm_turbine_vis = site_data$X150cm_turbine_vis_scaled,
  biotic_com_2 = site_data$biotic_com_2,
  shrub_yucca_density = site_data$shrub_yucca_density
))

site_data$Predicted <- site_preds$Predicted
site_data$lower <- site_preds$lower
site_data$upper <- site_preds$upper

# Plot
deer.vis <- ggplot(plot_data, 
                   aes(x = X150cm_turbine_vis_unscaled,
                       y = Predicted,
                       color = biotic_com_2,
                       fill = biotic_com_2,
                       linetype = biotic_com_2)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.45, color = NA) +
  geom_line(aes(y = lower), linewidth = 0.4) +
  geom_line(aes(y = upper), linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c(
    "Grassland" = "black",
    "Woodland"  = "black"
  )) +
  scale_fill_manual(values = c(
    "Grassland" = "#FDE725",
    "Woodland"  = "darkgray"
  )) +
  scale_linetype_manual(values = c(
    "Grassland" = "solid",
    "Woodland"  = "dashed"
  )) +
  scale_x_continuous(
    breaks = seq(0, 140, 5),
    labels = ifelse(seq(0, 140, 5) %% 20 == 0,
                    seq(0, 140, 5), "")
  ) +
  scale_y_continuous(
    breaks = seq(0, 1.00, 0.05),
    labels = ifelse(seq(0, 1.00, 0.05) %% 0.25 == 0,
                    sprintf("%.2f", seq(0, 1.00, 0.05)), "")
  ) +
  labs(
    x = "Turbine visibility",
    y = "Habitat use probability (ψ)",
    color = "Biotic community",
    fill = "Biotic community",
    linetype = "Biotic community"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14),
    axis.line  = element_line(size = 0.8),
    axis.ticks = element_line(size = 0.8),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.position = "right"
  )

deer.vis

#### Turbine distance + veg cover model plot ####

# Fit model using scaled wind variable
mod1 <- occu(~ max_trig_dist + livestock.count ~  
               veg_cover + turbine_dist, occu.odhe)

# Get mean and SD of wind variable
wind_mean <- mean(site.covs$turbine_dist, na.rm = TRUE)
wind_sd   <- sd(site.covs$turbine_dist, na.rm = TRUE)

# Generate sequence of unscaled values, then scale them
turbine_seq_unscaled <- seq(min(site.covs$turbine_dist, 
                                na.rm = TRUE),
                            max(site.covs$turbine_dist,
                                na.rm = TRUE),
                            length.out = 100)
turbine_seq_scaled <- (turbine_seq_unscaled - wind_mean) / wind_sd

# Create plot prediction dataset 
plot_data <- data.frame(
  turbine_dist = turbine_seq_scaled,
  veg_cover = 0,
  turbine_dist_unscaled = turbine_seq_unscaled
)

# Predict
preds <- predict(mod1, type = "state", newdata = plot_data)
plot_data$Predicted <- preds$Predicted
plot_data$lower <- preds$lower
plot_data$upper <- preds$upper

# Create site_data for observed points
site_data <- site.covs %>%
  dplyr::select(turbine_dist, veg_cover) %>%
  mutate(turbine_dist_scaled = 
           (turbine_dist - wind_mean) / wind_sd,
         veg_cover = 0
  )

# Predict for observed sites
site_preds <- predict(mod1, type = "state", newdata = data.frame(
  turbine_dist = site_data$turbine_dist_scaled,
  veg_cover = site_data$veg_cover
))

site_data$Predicted <- site_preds$Predicted
site_data$lower <- site_preds$lower
site_data$upper <- site_preds$upper

# Add species label for combined plots 
plot_data$species <- "Mule deer"
deer_df_dist <- plot_data

# Plot 
mule.dist <- ggplot(plot_data, aes(x = turbine_dist_unscaled, 
                                   y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, 
              fill = "darkgray") +
  geom_line(color = "black", linewidth = 0.7) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.05),
    labels = ifelse(
      seq(0, 1, by = 0.05) %% 0.25 == 0 |
        seq(0, 1, by = 0.05) == 0,   
      sprintf("%.2f", seq(0, 1, by = 0.05)),
      ""
    )
  ) +
  scale_x_continuous(    breaks = seq(0, 11000, by = 500),                     
                         labels = ifelse(seq(0, 11000, by = 500) %% 2000 == 0, 
                                         seq(0, 11000, by = 500), "")
  ) +
  labs(x = "Distance to nearest turbine (m)",
       y = "Habitat use probability (ψ)") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),   
    axis.text  = element_text(size = 14),  
    axis.line  = element_line(size = 0.8),  
    axis.ticks = element_line(size = 0.8)   
  )

mule.dist

########################################################
################### COYOTE PLOTS #######################
########################################################

#### Setup code for coyote models ####

# Read in edited .csv
detHist <- read.csv(file = "coyote detection hist.csv", 
                    row.names = 1)

# Change from integer to numeric
detHist %>%
  mutate(across(1:51, as.numeric))

# Create detection history 
occu.cala <- unmarkedFrameOccu(y=detHist, 
                               siteCovs = site.covs.scaled, 
                               obsCovs = obsCovs.scaled)

#### Turbine visibility model plot ####


# Fit model using scaled wind variable
mod1 <- occu( ~ people.active + cottontail.active ~ 
                X50cm_turbine_vis, occu.cala)

confint(mod1, type = "state", level = 0.85) 
confint(mod1, type = "det", level = 0.85)

# Get mean and SD of wind variable
wind_mean <- mean(site.covs$X50cm_turbine_vis, na.rm = TRUE)
wind_sd   <- sd(site.covs$X50cm_turbine_vis, na.rm = TRUE)

# Generate sequence of unscaled values, then scale them
turbine_seq_unscaled <- seq(min(site.covs$X50cm_turbine_vis, 
                                na.rm = TRUE),
                            max(site.covs$X50cm_turbine_vis,
                                na.rm = TRUE),
                            length.out = 100)
turbine_seq_scaled <- (turbine_seq_unscaled - wind_mean) / wind_sd

# Create plot prediction dataset 
plot_data <- data.frame(
  X50cm_turbine_vis = turbine_seq_scaled,
  X50cm_turbine_vis_unscaled = turbine_seq_unscaled
)

# Predict
preds <- predict(mod1, type = "state", newdata = plot_data)
plot_data$Predicted <- preds$Predicted
plot_data$lower <- preds$lower
plot_data$upper <- preds$upper

# Create site_data for observed points
site_data <- site.covs %>%
  dplyr::select(X50cm_turbine_vis) %>%
  mutate(X50cm_turbine_vis_scaled = 
           (X50cm_turbine_vis - wind_mean) / wind_sd
  )

# Predict for observed sites
site_preds <- predict(mod1, type = "state", newdata = data.frame(
  X50cm_turbine_vis = site_data$X50cm_turbine_vis_scaled
))

site_data$Predicted <- site_preds$Predicted
site_data$lower <- site_preds$lower
site_data$upper <- site_preds$upper

# Add species label for combined plots 
plot_data$species <- "Coyote"
coyote_df <- plot_data

# Plot 
coy.vis  <- ggplot(plot_data, aes(x = X50cm_turbine_vis_unscaled, 
                                  y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, 
              fill = "darkgray") +
  geom_line(color = "black", linewidth = 0.7) +
  scale_x_continuous(
    breaks = seq(0, 140, by = 5),                   # all ticks
    labels = ifelse(seq(0, 140, by = 5) %% 20 == 0, # label only multiples of 20
                    seq(0, 140, by = 5), ""), expand = expansion(mult = c(0, 0.02)) ) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Turbine visibility",
       y = "Habitat use probability (ψ)") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),   # axis labels bigger
    axis.text  = element_text(size = 14),   # tick labels bigger
    axis.line  = element_line(size = 0.8),  # slightly bolder axis lines
    axis.ticks = element_line(size = 0.8)   # bolder ticks
  )

coy.vis 

########################################################
################# GRAY FOX PLOTS #######################
########################################################
#### Setup code for gray fox models ####

# Read in edited .csv
detHist <- read.csv(file = "gray fox detection hist.csv", 
                    row.names = 1)

# Change from integer to numeric
detHist <- detHist %>%
  mutate(across(1:51, as.numeric))

# Create detection history 
occu.urci <- unmarkedFrameOccu(y = detHist, 
                               siteCovs = site.covs.scaled, 
                               obsCovs = obsCovs.scaled)

#### Turbine interior + slope + NDVI model plot ####

mod1 <- occu( ~ as.factor(cam_moved) 
              ~ as.factor(turbine_interior) +
                NDVI_1_5km + slope, occu.urci)

# New data for predictions
newdat <- data.frame(
  NDVI_1_5km = 0,
  slope = 0,
  turbine_interior = factor(c(0, 1))  # 0 = exterior, 1 = interior
)

# Predict detection probabilities
preds <- predict(mod1, type = "state", newdata = newdat)
plotdat <- cbind(newdat, preds)

# Add species label for combined plots 
grayfox_plotdat <- plotdat %>%
  dplyr::mutate(
    species = "Gray fox"
  )

# Rename factor levels 
plotdat$turbine_interior <- factor(plotdat$turbine_interior,
                                        levels = c(0, 1),
                                        labels = c("Exterior", "Interior"))

# Plot
gray.fox.int <- ggplot(plotdat, aes(x = turbine_interior, 
                                    y = Predicted, group = 1)) +
  geom_line(color = "black", linewidth = 0.8) +  
  geom_point(aes(color = turbine_interior), size = 5) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.2, linewidth = 0.9, color = "black") +  
  scale_color_manual(values = c("Exterior" = "#35B779", 
                                "Interior" = "steelblue")) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Wind farm proximity",
       y = "Probability of habitat selection (ψ)") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14),
    legend.position = "none"
  )

gray.fox.int

#### Turbine dist + NDVI + slope model plot ####

# Fit model using scaled wind variable
mod1 <- occu( ~ as.factor(cam_moved) ~
                turbine_dist + NDVI_1_5km + slope, occu.urci)

#Get mean and SD of wind variable
wind_mean <- mean(site.covs$turbine_dist, na.rm = TRUE)
wind_sd   <- sd(site.covs$turbine_dist, na.rm = TRUE)

# Generate sequence of unscaled values, then scale them
turbine_seq_unscaled <- seq(min(site.covs$turbine_dist, 
                                na.rm = TRUE),
                            max(site.covs$turbine_dist,
                                na.rm = TRUE),
                            length.out = 100)
turbine_seq_scaled <- (turbine_seq_unscaled - wind_mean) / wind_sd

# Create plot prediction dataset 
plot_data <- data.frame(
  turbine_dist = turbine_seq_scaled,
  slope = 0,
  NDVI_1_5km = 0,
  turbine_dist_unscaled = turbine_seq_unscaled
)

# Predict
preds <- predict(mod1, type = "state", newdata = plot_data)
plot_data$Predicted <- preds$Predicted
plot_data$lower <- preds$lower
plot_data$upper <- preds$upper

# Create site_data for observed points
site_data <- site.covs %>%
  dplyr::select( turbine_dist, slope, NDVI_1_5km) %>%
  mutate(turbine_dist_scaled = (turbine_dist - wind_mean) / wind_sd,
         slope = 0,
         NDVI_1_5km = 0
  )

# Predict for observed sites
site_preds <- predict(mod1, type = "state", newdata = data.frame(
  turbine_dist = site_data$turbine_dist_scaled,
  slope = site_data$slope,
  NDVI_1_5km = site_data$NDVI_1_5km
))

site_data$Predicted <- site_preds$Predicted
site_data$lower <- site_preds$lower
site_data$upper <- site_preds$upper

# Add species label for combined plots 
plot_data$species <- "Gray fox"
grayfox_df_dist <- plot_data

# Plot 
gray.fox.dist <- ggplot(plot_data, aes(x = turbine_dist_unscaled, 
                                       y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, 
              fill = "#00BA38") +
  geom_line(color = "#00BA38", linewidth = 0.7) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.05),
    labels = ifelse(
      seq(0, 1, by = 0.05) %% 0.25 == 0 |
        seq(0, 1, by = 0.05) == 0,   
      sprintf("%.2f", seq(0, 1, by = 0.05)),
      ""
    )
  ) +
  scale_x_continuous(    breaks = seq(0, 12000, by = 500),                     
                         labels = ifelse(seq(0, 12000, by = 500) %% 3000 == 0, 
                                         seq(0, 12000, by = 500), "")
  ) +
  labs(x = "Distance to nearest turbine (m)",
       y = "Probability of habitat selection (ψ)") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),  
    axis.text  = element_text(size = 14),  
    axis.line  = element_line(size = 0.8), 
    axis.ticks = element_line(size = 0.8)  
  )

gray.fox.dist 

#### Turbine density + slope + NDVI model plot ####

# Fit model using scaled wind variable
mod1 <- occu( ~ as.factor(cam_moved) ~
                turbine_density_1_5km + slope + NDVI_1_5km, occu.urci)

# Get mean and SD of wind variable
wind_mean <- mean(site.covs$turbine_density_1_5km, na.rm = TRUE)
wind_sd   <- sd(site.covs$turbine_density_1_5km, na.rm = TRUE)

# Generate sequence of unscaled values, then scale them
turbine_seq_unscaled <- seq(min(site.covs$turbine_density_1_5km, 
                                na.rm = TRUE),
                            max(site.covs$turbine_density_1_5km,
                                na.rm = TRUE),
                            length.out = 100)
turbine_seq_scaled <- (turbine_seq_unscaled - wind_mean) / wind_sd

# Create plot prediction dataset 
plot_data <- data.frame(
  turbine_density_1_5km = turbine_seq_scaled,
  X50cm_turbine_vis = 0,
  slope = 0,
  NDVI_1_5km = 0,
  turbine_density_1_5km_unscaled = turbine_seq_unscaled
)

# Predict
preds <- predict(mod1, type = "state", newdata = plot_data)
plot_data$Predicted <- preds$Predicted
plot_data$lower <- preds$lower
plot_data$upper <- preds$upper

# Create site_data for observed points
site_data <- site.covs %>%
  dplyr::select(X50cm_turbine_vis, turbine_density_1_5km, slope,
                NDVI_1_5km) %>%
  mutate(turbine_density_1_5km_scaled = 
           (turbine_density_1_5km - wind_mean) / wind_sd,
         X50cm_turbine_vis = 0,
         slope = 0,
         NDVI_1_5km = 0
  )

# Predict for observed sites
site_preds <- predict(mod1, type = "state", newdata = data.frame(
  turbine_density_1_5km = site_data$turbine_density_1_5km_scaled,
  X50cm_turbine_vis = site_data$X50cm_turbine_vis,
  slope = site_data$slope,
  NDVI_1_5km = site_data$NDVI_1_5km
))

site_data$Predicted <- site_preds$Predicted
site_data$lower <- site_preds$lower
site_data$upper <- site_preds$upper

breaks_x <- seq(0, 1.3, by = 0.05) # label every 4th break
labels_x <- ifelse(((seq_along(breaks_x) - 1) %% 4) == 0,
                   sprintf("%.1f", breaks_x),
                   "")

# Plot 
gray.fox.dense <- ggplot(plot_data, aes(x = turbine_density_1_5km_unscaled, 
                                        y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, 
              fill = "#35B779") + 
  geom_line(color = "#35B779", linewidth = 0.7) +
  scale_y_continuous(
    limits = c(0, 1)
  ) +
  scale_x_continuous(
    breaks = seq(0, 1.2, 0.2)
  ) +
  labs(x = "Turbine density",
       y = "Probability of habitat selection (ψ)") +
  theme_classic() +
  guides(fill = "none") +
  theme(
    legend.position = "none"
  )

gray.fox.dense

#### Turbine visibility + turbine density + slope additive model plot ####
# Hold turbine density at mean and plot turbine visibility #

# Fit model using scaled wind variable
mod1 <- occu( ~ as.factor(cam_moved) 
              ~ X50cm_turbine_vis + slope + turbine_density_1_5km, occu.urci)

# Get mean and SD of wind variable
wind_mean <- mean(site.covs$X50cm_turbine_vis, na.rm = TRUE)
wind_sd   <- sd(site.covs$X50cm_turbine_vis, na.rm = TRUE)

# Generate sequence of unscaled values, then scale them
turbine_seq_unscaled <- seq(min(site.covs$X50cm_turbine_vis, 
                                na.rm = TRUE),
                            max(site.covs$X50cm_turbine_vis,
                                na.rm = TRUE),
                            length.out = 100)
turbine_seq_scaled <- (turbine_seq_unscaled - wind_mean) / wind_sd

# Create plot prediction dataset 
plot_data <- data.frame(
  X50cm_turbine_vis = turbine_seq_scaled,
  slope = 0,
  turbine_density_1_5km = 0,
  X50cm_turbine_vis_unscaled = turbine_seq_unscaled
)

# Predict
preds <- predict(mod1, type = "state", newdata = plot_data)
plot_data$Predicted <- preds$Predicted
plot_data$lower <- preds$lower
plot_data$upper <- preds$upper

# Create site_data for observed points
site_data <- site.covs %>%
  dplyr::select(X50cm_turbine_vis, slope, turbine_density_1_5km) %>%
  mutate(X50cm_turbine_vis_scaled = 
           (X50cm_turbine_vis - wind_mean) / wind_sd,
         slope = 0, 
         turbine_density_1_5km = 0
  )

# Predict for observed sites
site_preds <- predict(mod1, type = "state", newdata = data.frame(
  X50cm_turbine_vis = site_data$X50cm_turbine_vis_scaled,
  turbine_density_1_5km = site_data$turbine_density_1_5km,
  slope = site_data$slope
))

site_data$Predicted <- site_preds$Predicted
site_data$lower <- site_preds$lower
site_data$upper <- site_preds$upper

# Add species label for combined plots 
plot_data$species <- "Gray fox"
grayfox_df_vis <- plot_data

# Plot 
gray.fox.vis <- ggplot(plot_data, aes(x = X50cm_turbine_vis_unscaled, 
                                      y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, 
              fill = "darkgray") +
  geom_line(color = "black", linewidth = 0.7) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(
    breaks = seq(0, 140, by = 5),                  
    labels = ifelse(seq(0, 140, by = 5) %% 20 == 0, 
                    seq(0, 140, by = 5), ""), 
    expand = expansion(mult = c(0, 0.02)) ) +
  labs(x = "Turbine visibility",
       y = "Habitat use probability (ψ)") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),  
    axis.text  = element_text(size = 14),   
    axis.line  = element_line(size = 0.8),  
    axis.ticks = element_line(size = 0.8)   
  )

gray.fox.vis 

########################################################
############## AMERICAN BADGER PLOTS ###################
########################################################

#### Setup code for American badger models ####

# Read in edited .csv
detHist <- read.csv(file = "badger detection hist.csv", 
                    row.names = 1)

# Change from integer to numeric

detHist <- detHist %>%
  mutate(across(1:51, as.numeric))

# Create detection history 
occu.tata <- unmarkedFrameOccu(y=detHist, 
                               siteCovs = site.covs.scaled, 
                               obsCovs = obsCovs.scaled)

#### Access rd density * bio comm + vertical cov + slope model plot ####

# Convert biotic_com_2 to factor
occu.tata@siteCovs$biotic_com_2 <- as.factor(occu.tata@siteCovs$biotic_com_2)
biotic_levels <- levels(occu.tata@siteCovs$biotic_com_2)

mod1 <- occu( ~ detection_angle + as.factor(precip.cat)
              ~ turbine_rd_density_1_6km *
                as.factor(biotic_com_2) + vertical_cover + slope, occu.tata)

# Get mean and SD of turbine_rd_density_1_6km
wind_mean <- mean(site.covs$turbine_rd_density_1_6km, na.rm = TRUE)
wind_sd   <- sd(site.covs$turbine_rd_density_1_6km, na.rm = TRUE)

# Generate sequence of unscaled values, then scale them
turbine_seq_unscaled <- seq(min(site.covs$turbine_rd_density_1_6km, 
                                na.rm = TRUE),
                            max(site.covs$turbine_rd_density_1_6km, 
                                na.rm = TRUE),
                            length.out = 100)
turbine_seq_scaled <- (turbine_seq_unscaled - wind_mean) / wind_sd

# Create plot prediction dataset
plot_data <- expand.grid(
  turbine_rd_density_1_6km  = turbine_seq_scaled,
  biotic_com_2 = factor(biotic_levels, levels = c("1", "2"), 
                        labels = c("Grassland", "Woodland")),
  slope = 0,
  vertical_cover = 0
)

plot_data <- plot_data %>% 
  dplyr::filter(biotic_com_2 == "Woodland")

# Predict
preds <- predict(mod1, type = "state", newdata = plot_data)
plot_data$Predicted <- preds$Predicted
plot_data$lower <- preds$lower
plot_data$upper <- preds$upper
plot_data$turbine_rd_density_1_6km_unscaled <- turbine_seq_unscaled

badger_plot <- plot_data %>%
  mutate(
    rd_density_unscaled = turbine_rd_density_1_6km_unscaled,
    species = "American badger"
  ) %>%
  dplyr::select(rd_density_unscaled, Predicted, lower, upper, species)

# Create site_data for observed points
site_data <- site.covs %>%
  dplyr::select(turbine_rd_density_1_6km , biotic_com_2, slope, vertical_cover) %>%
  mutate(
    turbine_rd_density_1_6km_scaled = 
      (turbine_rd_density_1_6km - wind_mean) / wind_sd,
    slope = 0,
    vertical_cover = 0
  )

# Recode biotic_com_2
site_data$biotic_com_2 <- factor(as.character(site_data$biotic_com_2),
                                 levels = c("1", "2"),
                                 labels = c("Grassland", "Woodland"))

site_data <- site_data %>% 
  dplyr::filter(biotic_com_2 == "Woodland")

# Predict for observed sites
site_preds <- predict(mod1, type = "state", newdata = data.frame(
  turbine_rd_density_1_6km = site_data$turbine_rd_density_1_6km_scaled,
  biotic_com_2 = site_data$biotic_com_2,
  slope = site_data$slope,
  vertical_cover = site_data$vertical_cover
))

site_data$Predicted <- site_preds$Predicted
site_data$lower <- site_preds$lower
site_data$upper <- site_preds$upper

# Plot
badger.rd.dense <- ggplot(plot_data, aes(x = turbine_rd_density_1_6km_unscaled, 
                                         y = Predicted,
                                         linetype = biotic_com_2, 
                                         fill = biotic_com_2)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.45, color = NA) +
  geom_line(aes(y = lower), linewidth = 0.4) +
  geom_line(aes(y = upper), linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c(
    "Grassland" = "black",
    "Woodland"  = "black"
  )) +
  scale_fill_manual(values = c(
    "Grassland" = "darkgray",
    "Woodland"  = "#F66B6A"
  )) +
  scale_linetype_manual(values = c(
    "Grassland" = "dashed",
    "Woodland"  = "solid"
  )) +
  scale_x_continuous(
    breaks = seq(0, 1.5, by = 0.1),                   
    labels = ifelse(seq(0, 1.5, by = 0.1) %% 0.5 == 0,
                    seq(0, 1.5, by = 0.1), "")) +
  scale_y_continuous(
    breaks = seq(0, 1.00, by = 0.05),
    labels = ifelse(
      seq(0, 1.00, by = 0.05) %% 0.25 == 0,
      sprintf("%.2f", seq(0, 1.00, by = 0.05)),
      ""
    )
  ) +
  labs(x = "Access road density within 1.6 km radius",
       y = "Probability of habitat selection (ψ)",
       linetype = "Biotic community",
       fill = "Biotic community") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14),
    axis.line  = element_line(size = 0.8),
    axis.ticks = element_line(size = 0.8),
    legend.text = element_text(size = 12), 
    legend.title = element_text(size = 13), 
    legend.position = "right"
  )

badger.rd.dense
  ## Warning message is okay, woodland only plotted as grassland uninformative

########################################################
############### STRIPED SKUNK PLOTS ####################
########################################################

#### Setup code for striped skunk models ####

# Read in edited .csv
detHist <- read.csv(file = "skunk detection hist.csv", 
                    row.names = 1)

# Change from integer to numeric
detHist %>%
  mutate(across(1:51, as.numeric))

# Create detection history 
occu.meme <- unmarkedFrameOccu(y=detHist, 
                               siteCovs = site.covs.scaled, 
                               obsCovs = obsCovs.scaled)

#### Access road density model plot ####

# Fit model using scaled wind variable
mod1 <- occu( ~ veg_cover_cam + woodland_percent_1_3km 
              ~  turbine_rd_density_1_3km, occu.meme)

# Get mean and SD of wind variable
wind_mean <- mean(site.covs$turbine_rd_density_1_3km, na.rm = TRUE)
wind_sd   <- sd(site.covs$turbine_rd_density_1_3km, na.rm = TRUE)

# Generate sequence of unscaled values, then scale them
turbine_seq_unscaled <- seq(min(site.covs$turbine_rd_density_1_3km, 
                                na.rm = TRUE),
                            max(site.covs$turbine_rd_density_1_3km,
                                na.rm = TRUE),
                            length.out = 100)
turbine_seq_scaled <- (turbine_seq_unscaled - wind_mean) / wind_sd

# Create plot prediction dataset 
plot_data <- data.frame(
  turbine_rd_density_1_3km = turbine_seq_scaled,
  turbine_rd_density_1_3km_unscaled = turbine_seq_unscaled
)

# Predict
preds <- predict(mod1, type = "state", newdata = plot_data)
plot_data$Predicted <- preds$Predicted
plot_data$lower <- preds$lower
plot_data$upper <- preds$upper

skunk_plot <- plot_data %>%
  mutate(
    rd_density_unscaled = turbine_rd_density_1_3km_unscaled,
    species = "Striped skunk"
  ) %>%
  dplyr::select(rd_density_unscaled, Predicted, lower, upper, species)


# Create site_data for observed points
site_data <- site.covs %>%
  dplyr::select(turbine_rd_density_1_3km) %>%
  mutate(turbine_rd_density_1_3km_scaled = 
           (turbine_rd_density_1_3km - wind_mean) / wind_sd
  )

# Predict for observed sites
site_preds <- predict(mod1, type = "state", newdata = data.frame(
  turbine_rd_density_1_3km = site_data$turbine_rd_density_1_3km_scaled
))

site_data$Predicted <- site_preds$Predicted
site_data$lower <- site_preds$lower
site_data$upper <- site_preds$upper

# Plot 

skunk.rd.dense <- ggplot(plot_data, aes(x = turbine_rd_density_1_3km_unscaled, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "#B4DD2C") +
  geom_line(color = "#B4DD2C", linewidth = 0.7) +
  scale_y_continuous(
    limits = c(0, 1) ) +
  scale_x_continuous(
    breaks = seq(0, 1.6, by = 0.4)) +
  labs(x = "Access road density within 1.3 km radius",
       y = "Probability of habitat selection (ψ)") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 12),
    axis.line  = element_line(size = 0.8),
    axis.ticks = element_line(size = 0.8),
    axis.ticks.length = unit(0.18, "cm")
  )

skunk.rd.dense

########################################################
########## BLACK-TAILED JACKRABBIT PLOTS ###############
########################################################
#### Setup code for black-tailed jackrabbit models ####

# Read in edited .csv
detHist <- read.csv(file = "jackrabbit detection hist.csv", 
                    row.names = 1)

# Change from integer to numeric
detHist %>%
  mutate(across(1:51, as.numeric))

# Create detection history 
occu.leca <- unmarkedFrameOccu(y=detHist, 
                               siteCovs = site.covs.scaled, 
                               obsCovs = obsCovs.scaled)

#### Access rd dist + vis + slope + woodland percent (keep vis constant) ####

# Fit model using scaled wind variable
mod1 <-  occu( ~ as.factor(precip.cat) + 
                 veg_cover_cam_under_1m
               ~ slope + woodland_percent_0_9km +
                 turbine_rd_dist + X50cm_turbine_vis, occu.leca)

# Get mean and SD of wind variable
wind_mean <- mean(site.covs$turbine_rd_dist, na.rm = TRUE)
wind_sd   <- sd(site.covs$turbine_rd_dist, na.rm = TRUE)

# Generate sequence of unscaled values, then scale them
turbine_seq_unscaled <- seq(min(site.covs$turbine_rd_dist, 
                                na.rm = TRUE),
                            max(site.covs$turbine_rd_dist,
                                na.rm = TRUE),
                            length.out = 100)
turbine_seq_scaled <- (turbine_seq_unscaled - wind_mean) / wind_sd

# Create plot prediction dataset 
plot_data <- data.frame(
  turbine_rd_dist = turbine_seq_scaled,
  woodland_percent_0_9km = 0,
  slope = 0,
  X50cm_turbine_vis = 0,
  turbine_rd_dist_unscaled = turbine_seq_unscaled
)

# Predict
preds <- predict(mod1, type = "state", newdata = plot_data)
plot_data$Predicted <- preds$Predicted
plot_data$lower <- preds$lower
plot_data$upper <- preds$upper

# Create site_data for observed points
site_data <- site.covs %>%
  dplyr::select(turbine_rd_dist, slope, woodland_percent_0_9km, 
                X50cm_turbine_vis) %>%
  mutate(turbine_rd_dist_scaled = 
           (turbine_rd_dist - wind_mean) / wind_sd,
         woodland_percent_0_9km = 0,
         slope = 0,
         X50cm_turbine_vis = 0
  )

# Predict for observed sites
site_preds <- predict(mod1, type = "state", newdata = data.frame(
  turbine_rd_dist = site_data$turbine_rd_dist_scaled,
  slope = site_data$slope,
  woodland_percent_0_9km = site_data$woodland_percent_0_9km,
  X50cm_turbine_vis = site_data$X50cm_turbine_vis 
))

site_data$Predicted <- site_preds$Predicted
site_data$lower <- site_preds$lower
site_data$upper <- site_preds$upper

# Add species label for combined plots 
plot_data$species <- "Black-tailed jackrabbit"
jack_df <- plot_data

# Plot 

jack.rd.dist <- ggplot(plot_data, aes(x = turbine_rd_dist_unscaled, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "#31688E") +
  geom_line(color = "#31688E", linewidth = 1) +
  scale_x_continuous(
    breaks = seq(0, 11000, by = 2000),                     
    labels = ifelse(seq(0, 11000, by = 2000) %% 2000 == 0, 
                    seq(0, 11000, by = 2000), "")
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Distance to nearest access road (m)",
       y = "Probability of habitat selection (ψ)") +
  theme_classic() 


jack.rd.dist

#### Access rd dist + vis + slope + woodland percent (keep rd dist constant) ####

# Fit model using scaled wind variable
mod1 <-  occu( ~ as.factor(precip.cat) + 
                 veg_cover_cam_under_1m
               ~ slope + woodland_percent_0_9km +
                 turbine_rd_dist + X50cm_turbine_vis, occu.leca)

# Get mean and SD of wind variable
wind_mean <- mean(site.covs$X50cm_turbine_vis, na.rm = TRUE)
wind_sd   <- sd(site.covs$X50cm_turbine_vis, na.rm = TRUE)

# Generate sequence of unscaled values, then scale them
turbine_seq_unscaled <- seq(min(site.covs$X50cm_turbine_vis, 
                                na.rm = TRUE),
                            max(site.covs$X50cm_turbine_vis,
                                na.rm = TRUE),
                            length.out = 100)
turbine_seq_scaled <- (turbine_seq_unscaled - wind_mean) / wind_sd

# Create plot prediction dataset 
plot_data <- data.frame(
  X50cm_turbine_vis = turbine_seq_scaled,
  woodland_percent_0_9km = 0,
  slope = 0,
  turbine_rd_dist = 0,
  X50cm_turbine_vis_unscaled = turbine_seq_unscaled
)

# Predict
preds <- predict(mod1, type = "state", newdata = plot_data)
plot_data$Predicted <- preds$Predicted
plot_data$lower <- preds$lower
plot_data$upper <- preds$upper

# Create site_data for observed points
site_data <- site.covs %>%
  dplyr::select(X50cm_turbine_vis, slope, woodland_percent_0_9km, 
                turbine_rd_dist) %>%
  mutate(X50cm_turbine_vis_scaled = 
           (X50cm_turbine_vis - wind_mean) / wind_sd,
         woodland_percent_0_9km = 0,
         slope = 0,
         turbine_rd_dist = 0
  )

# Predict for observed sites
site_preds <- predict(mod1, type = "state", newdata = data.frame(
  X50cm_turbine_vis = site_data$X50cm_turbine_vis_scaled,
  slope = site_data$slope,
  woodland_percent_0_9km = site_data$woodland_percent_0_9km,
  turbine_rd_dist = site_data$turbine_rd_dist
))

site_data$Predicted <- site_preds$Predicted
site_data$lower <- site_preds$lower
site_data$upper <- site_preds$upper

# Add species label for combined plots 
plot_data$species <- "Black-tailed jackrabbit"
jack_df_vis <- plot_data

# Plot 
jack.vis <- ggplot(plot_data, aes(x = X50cm_turbine_vis_unscaled, 
                                  y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, 
              fill = "darkgray") +
  geom_line(color = "black", linewidth = 0.7) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(
    breaks = seq(0, 140, by = 5),  
    expand = expansion(add = c(0, 0.05)),
    labels = ifelse(seq(0, 140, by = 5) %% 20 == 0, 
                    seq(0, 140, by = 5), "")) +
  labs(x = "Turbine visibility",
       y = "Occupancy probability (ψ)") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),  
    axis.text  = element_text(size = 14),   
    axis.line  = element_line(size = 0.8), 
    axis.ticks = element_line(size = 0.8)   
  )

jack.vis

#### Turbine vis + turbine interior + slope + woodland percent ####
# keep vis constant

# Fit model using scaled wind variable
mod1 <-  occu( ~ as.factor(precip.cat) + veg_cover_cam_under_1m
               ~ slope + woodland_percent_0_9km + as.factor(turbine_interior) +
                 X50cm_turbine_vis, occu.leca)

# New data for predictions
newdat <- data.frame(
  woodland_percent_0_9km  = 0,
  slope = 0,
  X50cm_turbine_vis = 0,
  turbine_interior = factor(c(0, 1))  # 0 = exterior, 1 = interior
)

# Predict detection probabilities
preds <- predict(mod1, type = "state", newdata = newdat)
plotdat <- cbind(newdat, preds)

# Add species label for combined plots 
btjr_plotdat <- plotdat %>%
  dplyr::mutate(
    species = "Black-tailed jackrabbit"
  )

# Rename factor levels 
plotdat$turbine_interior <- factor(plotdat$turbine_interior,
                                   levels = c(0, 1),
                                   labels = c("Exterior", "Interior"))

# Plot
btjr.int <- ggplot(plotdat, aes(x = turbine_interior, 
                                    y = Predicted, group = 1)) +
  geom_line(color = "black", linewidth = 0.8) +  
  geom_point(aes(color = turbine_interior), size = 5) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.2, linewidth = 0.9, color = "black") +  
  scale_color_manual(values = c("Exterior" = "#35B779", 
                                "Interior" = "steelblue")) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Wind farm proximity",
       y = "Probability of habitat selection (ψ)") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14),
    legend.position = "none"
  )

btjr.int

########################################################
############# COMBINED SPECIES PLOTS ###################
########################################################

species_cols <- c(
  "Mule deer" = "#FDE725",
  "Coyote" = "#440154",
  "Gray fox" = "#35B779",
  "Striped skunk" = "#B4DD2C",
  "Black-tailed jackrabbit" = "#31688E",
  "American badger" = "#1F9E89"
)

species_scales <- list(
  scale_color_manual(
    values = species_cols,
    limits = names(species_cols)
  ),
  scale_fill_manual(
    values = species_cols,
    limits = names(species_cols)
  )
)

# ================================
# TURBINE INTERIOR PLOT
# ================================

combined_plotdat <- dplyr::bind_rows(btjr_plotdat, grayfox_plotdat)

combined_plotdat$turbine_interior <- factor(
  combined_plotdat$turbine_interior,
  levels = c(0, 1),
  labels = c("Exterior", "Interior")
)

int_plot <-  ggplot(combined_plotdat,
                    aes(x = turbine_interior,
                        y = Predicted,
                        group = species,
                        color = species)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.15, linewidth = 0.9) +
  scale_color_manual(values = species_cols) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Wind farm proximity",
       y = "Probability of habitat selection (ψ)",
       color = "Species") +
  theme_classic() 

int_plot

# ================================
# TURBINE VISIBILITY PLOT
# ================================

deer_df_vis <- deer_df_vis |>
  rename(X50cm_turbine_vis_unscaled = X150cm_turbine_vis_unscaled)

deer_grass_df <- deer_df_vis |>
  filter(biotic_com_2 == "Grassland")

turbine_vis_all <- bind_rows(
  grayfox_df_vis,
  coyote_df,
  jack_df_vis,
  deer_grass_df
)

vis_plot <- ggplot(
  turbine_vis_all,
  aes(X50cm_turbine_vis_unscaled, Predicted,
      color = species, fill = species)
) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.25, color = NA) +
  geom_line(linewidth = 1) +
  species_scales +
  scale_x_continuous(
    breaks = seq(0, 140, 20)
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Turbine visibility",
    y = "Probability of habitat selection (ψ)"
  ) +
  theme_classic() +
  guides(fill = "none") +
  theme(
    legend.position = "none"
  )

vis_plot

# ================================
# DISTANCE TO TURBINE PLOT
# ================================

turbine_dist <- bind_rows(grayfox_df_dist, deer_df_dist)

dist_plot <- ggplot(
  turbine_dist,
  aes(turbine_dist_unscaled, Predicted,
      color = species, fill = species)
) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.25, color = NA) +
  geom_line(linewidth = 1) +
  species_scales +
  scale_x_continuous(
    breaks = seq(0, 11000, 2000)
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Distance to nearest turbine (m)",
    y = "Probability of habitat selection (ψ)"
  ) +
  theme_classic() +
  guides(fill = "none") +
  theme(
    legend.position = "none"
  )

dist_plot

# ================================
# ACCESS RD DENSITY PLOT
# ================================

combined_plot <- bind_rows(badger_plot, skunk_plot)

rd_dense_plot <-  ggplot(combined_plot,
       aes(x = rd_density_unscaled,
           y = Predicted,
           color = species,
           fill = species)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.25, color = NA) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = species_cols) +
  scale_fill_manual(values = species_cols) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 1.6, by = 0.4)) +
  labs(x = "Access road density",
       y = "Probability of habitat selection (ψ)",
       color = "Species",
       fill  = "Species") +
  
  theme_classic()

rd_dense_plot

# ================================
# COMBINE FIRST 4 PLOTS
# ================================

# Remove plot legends
int_plot       <- int_plot          +  
  theme(
  legend.position = "none",
  axis.text  = element_text(size = 10),
  axis.title = element_text(size = 12)
)
vis_plot       <- vis_plot          +   
  theme(
    legend.position = "none",
    axis.text  = element_text(size = 10),
    axis.title = element_text(size = 12)
  )
dist_plot      <- dist_plot         +  
  theme(
    legend.position = "none",
    axis.text  = element_text(size = 10),
    axis.title = element_text(size = 12)
  )
gray.fox.dense <- gray.fox.dense   +  
  theme(
    legend.position = "none",
    axis.text  = element_text(size = 10),
    axis.title = element_text(size = 12)
)

# Create a dummy plot for the legend 
legend_plot <- ggplot(
  data.frame(
    species = factor(names(species_cols), levels = names(species_cols)),
    x = 1,
    y = 1
  ),
  aes(x, y, color = species)
) +
  geom_line(linewidth = 1) +
  scale_color_manual(
    values = species_cols,
    guide = guide_legend(ncol = 1)  
  ) +
  labs(color = "Species") +
  theme_classic() +
  theme(
    legend.position = "right",   
    legend.title = element_text(hjust = 0.5),
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    panel.background  = element_rect(fill = "transparent"),
    plot.background   = element_rect(fill = "transparent")
  )

# Extract the legend
shared_legend <- get_legend(legend_plot)

# Combine the four panels 
panel_grid <- plot_grid(
  int_plot,
  vis_plot,
  dist_plot,
  gray.fox.dense,
  ncol = 2,
  align = "hv"
)

# Add legend to the right 
final_figure <- plot_grid(
  panel_grid,
  ncol = 2,                  
  rel_widths = c(1, 0.2)     
)

# Draw final figure
final_figure

# Save layout to .png
ggsave(paste0(homewd, "figures/turbine.plots.png"), final_figure, 
       width = 8, height = 7.5, dpi = 300)

# ================================
# COMBINE LAST 2 PLOTS
# ================================

# Remove plots legends 
jack.rd.dist <- jack.rd.dist +
  theme(
    legend.position = "none",
    axis.text  = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

rd_dense_plot <- rd_dense_plot +
  theme(
    legend.position = "none",
    axis.text  = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

# Create a dummy plot for the legend 
legend_plot <- ggplot(
  data.frame(
    species = factor(names(species_cols), levels = names(species_cols)),
    x = 1,
    y = 1
  ),
  aes(x, y, color = species)
) +
  geom_line(linewidth = 1) +
  scale_color_manual(
    values = species_cols,
    guide = guide_legend(ncol = 1)  
  ) +
  labs(color = "Species") +
  theme_classic() +
  theme(
    legend.position = "right",   
    legend.title = element_text(hjust = 0.5),
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    panel.background  = element_rect(fill = "transparent"),
    plot.background   = element_rect(fill = "transparent")
  )

# Extract the legend 
shared_legend <- get_legend(legend_plot)

# Combine the 2 panels 
panel_grid <- plot_grid(
  jack.rd.dist,
  rd_dense_plot,
  ncol = 2,
  align = "hv"
)

# Add legend to the right 
final_figure_2 <- plot_grid(
  panel_grid,
  ncol = 2,                  
  rel_widths = c(1, 0.2)     
)

# Draw final figure 
final_figure_2

# Save layout to .png
ggsave(paste0(homewd, "figures/road.plots.png"), final_figure_2,  
       width = 8, height = 3.5, dpi = 300)

# Note that the two figures were combined to produce Figure 2 in the manuscript

########################################################
######################### END ##########################
########################################################