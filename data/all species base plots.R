###########################################################
###### Code for the creation of the plots showing the #####
### effects on the the non-wind energy variables on the ###
#  habitat selection of terrestrial mammals in central NM #
###########################################################
################ Script by: Iona Rohan ####################
############ Contact: ionarohan12@gmail.com ###############
###########################################################
########## Date Last Modified: 09-March-2026 ##############
###########################################################

###########################################################
####### NOTE:RUN THIS CODE AFTER "pre_model_code.R" #######
###########################################################

#Clear work environment
rm(list=ls())

#Note: If you opened this script through the .Rproj file, the only line you 
#should need to change for the script to run (assuming packages are installed) 
#is the homewd directory on line 24.

#Set home working directory
#homewd = "C:/Users/ionar/Desktop/R Repository/Wind Energy/"
homewd = "<insert your folder here and end with a forward slash>"

#Set wd to data folder on your local computer 
setwd(paste0(homewd, "data/"))

###########################################################
# MULE DEER BASE MODELS #
###########################################################

#### Setup code for mule deer occupancy models ####

# Read in edited .csv
detHist <- read.csv(file = "mule deer detection hist.csv", row.names = 1)

# Change from integer to numeric
detHist %>%
  mutate(across(1:51, as.numeric))

# Create detection history 
occu.odhe <- unmarkedFrameOccu(y=detHist, 
                               siteCovs = site.covs.scaled, 
                               obsCovs = obsCovs.scaled)

# ================================
# Mule deer detection variables 
# ================================

# Top detection model

mod1 <- occu( ~ max_trig_dist + NDVI_1_9km ~ 1, occu.odhe)

#### Plot max trigger distance ####

base_mean <- mean(site.covs$max_trig_dist, na.rm = TRUE)
base_sd   <- sd(site.covs$max_trig_dist, na.rm = TRUE)

# Generate sequence of unscaled values then scale them
seq_unscaled <- seq(min(site.covs$max_trig_dist, 
                        na.rm = TRUE),
                    max(site.covs$max_trig_dist,
                        na.rm = TRUE),
                    length.out = 100)
seq_scaled <- (seq_unscaled - base_mean) / base_sd

# Create plot prediction dataset 
plot_data <- data.frame(
  max_trig_dist = seq_scaled,
  max_trig_dist_unscaled = seq_unscaled,
  NDVI_1_9km = 0
)

# Predict
preds <- predict(mod1, type = "det", newdata = plot_data)
plot_data$Predicted <- preds$Predicted
plot_data$lower <- preds$lower
plot_data$upper <- preds$upper

# Create site_data for observed points
site_data <- site.covs %>%
  dplyr::select(NDVI_1_9km, max_trig_dist) %>%
  mutate(max_trig_dist_scaled = 
           (max_trig_dist - base_mean) / base_sd,
         NDVI_1_9km = 0
  )

# Predict for observed sites
site_preds <- predict(mod1, type = "det", newdata = data.frame(
  max_trig_dist = site_data$max_trig_dist_scaled,
  NDVI_1_9km = site_data$NDVI_1_9km
))

site_data$Predicted <- site_preds$Predicted
site_data$lower <- site_preds$lower
site_data$upper <- site_preds$upper

# Plot 
deer.trig  <- ggplot(plot_data, aes(x = max_trig_dist_unscaled, 
                                     y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, 
              fill = "darkgray") +
  geom_line(color = "black", linewidth = 0.7) +
  scale_x_continuous(breaks = pretty_breaks(n = 6)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Maximum trigger distance (m)",
       y = "Probability of detection (p)") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),   # axis labels bigger
    axis.text  = element_text(size = 14),   # tick labels bigger
    axis.line  = element_line(size = 0.8),  # slightly bolder axis lines
    axis.ticks = element_line(size = 0.8)   # bolder ticks
  )

deer.trig

#### Plot NDVI ####

base_mean <- mean(site.covs$NDVI_1_9km, na.rm = TRUE)
base_sd   <- sd(site.covs$NDVI_1_9km, na.rm = TRUE)

# Generate sequence of unscaled values, then scale them
seq_unscaled <- seq(min(site.covs$NDVI_1_9km, 
                        na.rm = TRUE),
                    max(site.covs$NDVI_1_9km,
                        na.rm = TRUE),
                    length.out = 100)
seq_scaled <- (seq_unscaled - base_mean) / base_sd

# Create plot prediction dataset 
plot_data <- data.frame(
  NDVI_1_9km = seq_scaled,
  NDVI_1_9km_unscaled = seq_unscaled,
  max_trig_dist = 0
)

# Predict
preds <- predict(mod1, type = "det", newdata = plot_data)
plot_data$Predicted <- preds$Predicted
plot_data$lower <- preds$lower
plot_data$upper <- preds$upper

# Create site_data for observed points
site_data <- site.covs %>%
  dplyr::select(NDVI_1_9km, max_trig_dist) %>%
  mutate(NDVI_1_9km_scaled = 
           (NDVI_1_9km - base_mean) / base_sd,
         max_trig_dist = 0
  )

# Predict for observed sites
site_preds <- predict(mod1, type = "det", newdata = data.frame(
  NDVI_1_9km = site_data$NDVI_1_9km_scaled,
  max_trig_dist = site_data$max_trig_dist
))

site_data$Predicted <- site_preds$Predicted
site_data$lower <- site_preds$lower
site_data$upper <- site_preds$upper

# Plot 
deer.ndvi  <- ggplot(plot_data, aes(x = NDVI_1_9km_unscaled, 
                                    y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, 
              fill = "darkgray") +
  geom_line(color = "black", linewidth = 0.7) +
  scale_x_continuous(breaks = pretty_breaks(n = 6)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "NDVI within 1.9 km radius",
       y = "Probability of detection (p)") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),   # axis labels bigger
    axis.text  = element_text(size = 14),   # tick labels bigger
    axis.line  = element_line(size = 0.8),  # slightly bolder axis lines
    axis.ticks = element_line(size = 0.8)   # bolder ticks
  )

deer.ndvi

# ================================
# Mule deer habitat selection variables 
# ================================

# Full mule deer base model 
mod1 <- occu( ~ max_trig_dist + NDVI_1_9km  
              ~ shrub_yucca_density + veg_cover, occu.odhe)

#### Plot shrub density ####

base_mean <- mean(site.covs$shrub_yucca_density, na.rm = TRUE)
base_sd   <- sd(site.covs$shrub_yucca_density, na.rm = TRUE)

# Generate sequence of unscaled values, then scale them
seq_unscaled <- seq(min(site.covs$shrub_yucca_density, 
                        na.rm = TRUE),
                    max(site.covs$shrub_yucca_density,
                        na.rm = TRUE),
                    length.out = 100)
seq_scaled <- (seq_unscaled - base_mean) / base_sd

# Create plot prediction dataset 
plot_data <- data.frame(
  shrub_yucca_density = seq_scaled,
  shrub_yucca_density_unscaled = seq_unscaled,
  veg_cover = 0
)

# Predict
preds <- predict(mod1, type = "state", newdata = plot_data)
plot_data$Predicted <- preds$Predicted
plot_data$lower <- preds$lower
plot_data$upper <- preds$upper

# Create site_data for observed points
site_data <- site.covs %>%
  dplyr::select(shrub_yucca_density, veg_cover) %>%
  mutate(shrub_yucca_density_scaled = 
           (shrub_yucca_density - base_mean) / base_sd,
         veg_cover = 0
  )

# Predict for observed sites
site_preds <- predict(mod1, type = "state", newdata = data.frame(
  shrub_yucca_density = site_data$shrub_yucca_density_scaled,
  veg_cover = site_data$veg_cover
))

site_data$Predicted <- site_preds$Predicted
site_data$lower <- site_preds$lower
site_data$upper <- site_preds$upper

# Plot 
deer.shrub  <- ggplot(plot_data, aes(x = shrub_yucca_density_unscaled, 
                                    y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, 
              fill = "darkgray") +
  geom_line(color = "black", linewidth = 0.7) +
  scale_x_continuous(breaks = pretty_breaks(n = 6)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Shrub density",
       y = "Probability of habitat selection (ψ)") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),   # axis labels bigger
    axis.text  = element_text(size = 14),   # tick labels bigger
    axis.line  = element_line(size = 0.8),  # slightly bolder axis lines
    axis.ticks = element_line(size = 0.8)   # bolder ticks
  )

deer.shrub

#### Plot visual obstruction ####

base_mean <- mean(site.covs$veg_cover, na.rm = TRUE)
base_sd   <- sd(site.covs$veg_cover, na.rm = TRUE)

# Generate sequence of unscaled values, then scale them
seq_unscaled <- seq(min(site.covs$veg_cover, 
                        na.rm = TRUE),
                    max(site.covs$veg_cover,
                        na.rm = TRUE),
                    length.out = 100)
seq_scaled <- (seq_unscaled - base_mean) / base_sd

# Create plot prediction dataset 
plot_data <- data.frame(
  veg_cover = seq_scaled,
  veg_cover_unscaled = seq_unscaled,
  shrub_yucca_density = 0
)

# Predict
preds <- predict(mod1, type = "state", newdata = plot_data)
plot_data$Predicted <- preds$Predicted
plot_data$lower <- preds$lower
plot_data$upper <- preds$upper

# Create site_data for observed points
site_data <- site.covs %>%
  dplyr::select(shrub_yucca_density, veg_cover) %>%
  mutate(veg_cover_scaled = 
           (veg_cover - base_mean) / base_sd,
         shrub_yucca_density = 0
  )

# Predict for observed sites
site_preds <- predict(mod1, type = "state", newdata = data.frame(
  veg_cover = site_data$veg_cover_scaled,
  shrub_yucca_density = site_data$shrub_yucca_density
))

site_data$Predicted <- site_preds$Predicted
site_data$lower <- site_preds$lower
site_data$upper <- site_preds$upper

# Plot 
deer.vis.obs  <- ggplot(plot_data, aes(x = veg_cover_unscaled, 
                                     y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, 
              fill = "darkgray") +
  geom_line(color = "black", linewidth = 0.7) +
  scale_x_continuous(breaks = pretty_breaks(n = 6)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Vegetation cover",
       y = "Probability of habitat selection (ψ)") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),   # axis labels bigger
    axis.text  = element_text(size = 14),   # tick labels bigger
    axis.line  = element_line(size = 0.8),  # slightly bolder axis lines
    axis.ticks = element_line(size = 0.8)   # bolder ticks
  )

deer.vis.obs

# Create combined plot layout 

final_plot <-
  (deer.trig  | deer.ndvi ) /
  (deer.shrub | deer.vis.obs)   +
  plot_layout(widths = c(1,1), heights = c(1,1)) +
  plot_annotation(
    title = "Mule deer",
    tag_levels = "a",          
    tag_prefix = "(",           
    tag_suffix = ")",           
    theme = theme(
      plot.tag = element_text(size = 40, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 15),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
      plot.tag.position = c(0.02, 0.95),  
      axis.title  = element_text(size = 20),
      axis.text   = element_text(size = 10)
    ) 
  )

final_plot

# Save as .png
#ggsave("mule_deer_base_plots.png", final_plot, width = 12, 
#height = 12, dpi = 300)


###########################################################
# PRONGHORN BASE MODELS #
###########################################################

#### Setup code for pronghorn occupancy models ####

# Read in edited .csv
detHist <- read.csv(file = "pronghorn detection hist.csv", row.names = 1)

# Change from integer to numeric
detHist %>%
  mutate(across(1:51, as.numeric))

# Create detection history 
occu.anam <- unmarkedFrameOccu(y=detHist, 
                               siteCovs = site.covs.scaled, 
                               obsCovs = obsCovs.scaled)

# ================================
# Pronghorn detection variables 
# ================================

# Top detection model
mod1 <- occu( ~ veg_cover_cam ~ 1, occu.anam)

#### Visual obstruction transect 1 plot ####

base_mean <- mean(site.covs$veg_cover_cam, na.rm = TRUE)
base_sd   <- sd(site.covs$veg_cover_cam, na.rm = TRUE)

# Generate sequence of unscaled values, then scale them
seq_unscaled <- seq(min(site.covs$veg_cover_cam, 
                        na.rm = TRUE),
                    max(site.covs$veg_cover_cam,
                        na.rm = TRUE),
                    length.out = 100)
seq_scaled <- (seq_unscaled - base_mean) / base_sd

# Create plot prediction dataset 
plot_data <- data.frame(
  veg_cover_cam = seq_scaled,
  veg_cover_cam_unscaled = seq_unscaled
)

# Predict
preds <- predict(mod1, type = "det", newdata = plot_data)
plot_data$Predicted <- preds$Predicted
plot_data$lower <- preds$lower
plot_data$upper <- preds$upper

# Create site_data for observed points
site_data <- site.covs %>%
  dplyr::select(veg_cover_cam) %>%
  mutate(veg_cover_cam_scaled = 
           (veg_cover_cam - base_mean) / base_sd
  )

# Predict for observed sites
site_preds <- predict(mod1, type = "det", newdata = data.frame(
  veg_cover_cam = site_data$veg_cover_cam_scaled
))

site_data$Predicted <- site_preds$Predicted
site_data$lower <- site_preds$lower
site_data$upper <- site_preds$upper

# Plot 
prong.vis.obs  <- ggplot(plot_data, aes(x = veg_cover_cam_unscaled, 
                                    y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, 
              fill = "darkgray") +
  geom_line(color = "black", linewidth = 0.7) +
  scale_x_continuous(breaks = pretty_breaks(n = 6)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Veg cover at camera",
       y = "Probability of detection (p)") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),   # axis labels bigger
    axis.text  = element_text(size = 14),   # tick labels bigger
    axis.line  = element_line(size = 0.8),  # slightly bolder axis lines
    axis.ticks = element_line(size = 0.8)   # bolder ticks
  )

prong.vis.obs

# ================================
# Pronghorn habitat selection variables 
# ================================

# Full pronghorn base model 
mod1 <- occu( ~ veg_cover_cam ~ slope, occu.anam)

base_mean <- mean(site.covs$slope, na.rm = TRUE)
base_sd   <- sd(site.covs$slope, na.rm = TRUE)

# Generate sequence of unscaled values, then scale them
seq_unscaled <- seq(min(site.covs$slope, 
                        na.rm = TRUE),
                    max(site.covs$slope,
                        na.rm = TRUE),
                    length.out = 100)
seq_scaled <- (seq_unscaled - base_mean) / base_sd

# Create plot prediction dataset 
plot_data <- data.frame(
  slope = seq_scaled,
  slope_unscaled = seq_unscaled
)

# Predict
preds <- predict(mod1, type = "state", newdata = plot_data)
plot_data$Predicted <- preds$Predicted
plot_data$lower <- preds$lower
plot_data$upper <- preds$upper

# Create site_data for observed points
site_data <- site.covs %>%
  dplyr::select(slope) %>%
  mutate(slope_scaled = 
           (slope - base_mean) / base_sd
  )

# Predict for observed sites
site_preds <- predict(mod1, type = "state", newdata = data.frame(
  slope = site_data$slope_scaled
))

site_data$Predicted <- site_preds$Predicted
site_data$lower <- site_preds$lower
site_data$upper <- site_preds$upper

# Plot 
prong.slope  <- ggplot(plot_data, aes(x = slope_unscaled, 
                                       y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, 
              fill = "darkgray") +
  geom_line(color = "black", linewidth = 0.7) +
  scale_x_continuous(breaks = pretty_breaks(n = 6)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Slope (degrees)",
       y = "Probability of habitat selection (ψ)") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),   # axis labels bigger
    axis.text  = element_text(size = 14),   # tick labels bigger
    axis.line  = element_line(size = 0.8),  # slightly bolder axis lines
    axis.ticks = element_line(size = 0.8)   # bolder ticks
  )

prong.slope

# Create combined plot layout 

final_plot <-
  (prong.vis.obs ) /
  (prong.slope )  +
  plot_layout(widths = c(1,1), heights = c(1,1)) +
  plot_annotation(
    title = "Pronghorn",
    tag_levels = "a",          
    tag_prefix = "(",           
    tag_suffix = ")",           
    theme = theme(
      plot.tag = element_text(size = 40, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 15),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
      plot.tag.position = c(0.02, 0.95),  
      axis.title  = element_text(size = 20),
      axis.text   = element_text(size = 10)
    ) 
  )

final_plot 

# Save layout as .png
#ggsave("prong_base_plots.png", final_plot, width = 12, height = 12, dpi = 300)


###########################################################
# COYOTE BASE MODELS #
###########################################################

#### Setup code for coyote occupancy models ####

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

# ================================
# Coyote detection variables 
# ================================

# Top detection model

mod1 <- occu( ~ human.active + cow.active 
              ~ 1, occu.cala, starts = c(2, 0, -3, 0))

#### Human activity and cow activity plots ####

#  Compute mean & SD from original unscaled obsCovs
human_vec <- as.numeric(as.matrix(obsCovs$human.active))
cow_vec   <- as.numeric(as.matrix(obsCovs$cow.active))

human_mean <- mean(human_vec, na.rm = TRUE)
human_sd   <- sd(human_vec, na.rm = TRUE)

cow_mean <- mean(cow_vec, na.rm = TRUE)
cow_sd   <- sd(cow_vec, na.rm = TRUE)


# Generate sequences of unscaled values
seq_human_unscaled <- seq(min(human_vec, na.rm = TRUE),
                          max(human_vec, na.rm = TRUE),
                          length.out = 100)
seq_cow_unscaled   <- seq(min(cow_vec, na.rm = TRUE),
                          max(cow_vec, na.rm = TRUE),
                          length.out = 100)

# Scale them for prediction
seq_human_scaled <- (seq_human_unscaled - human_mean) / human_sd
seq_cow_scaled   <- (seq_cow_unscaled - cow_mean) / cow_sd

# Build prediction datasets
pred_human_df <- data.frame(
  human.active = seq_human_scaled,
  cow.active   = 0  # hold cow at mean
)

pred_cow_df <- data.frame(
  human.active = 0,  # hold human at mean
  cow.active   = seq_cow_scaled
)

# Predict detection probability
pred_human <- predict(mod1, type = "det", newdata = pred_human_df)
pred_human$human_unscaled <- seq_human_unscaled

pred_cow <- predict(mod1, type = "det", newdata = pred_cow_df)
pred_cow$cow_unscaled <- seq_cow_unscaled

# Plot human activity effect
coy.human <- ggplot(pred_human, aes(x = human_unscaled, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, 
              fill = "darkgray") +
  geom_line(color = "black", linewidth = 0.7) +
  scale_x_continuous(breaks = pretty_breaks(n = 6)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Hours humans detected", y = "Probability of detection (p)") +
  theme_classic() +
  theme(axis.title = element_text(size = 16),
        axis.text  = element_text(size = 14))

coy.human

# Plot cow activity effect

coy.cow <- ggplot(pred_cow, aes(x = cow_unscaled, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, 
              fill = "darkgray") +
  geom_line(color = "black", linewidth = 0.7) +
  scale_x_continuous(breaks = pretty_breaks(n = 6)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Hours cattle detected", y = "Probability of detection (p)") +
  theme_classic() +
  theme(axis.title = element_text(size = 16),
        axis.text  = element_text(size = 14))

coy.cow

# Create combined plot layout 

final_plot <-
  (coy.human | coy.cow) +
  plot_layout(widths = c(1,1), heights = c(1,1)) +
  plot_annotation(
    title = "Coyote",
    tag_levels = "a",          
    tag_prefix = "(",           
    tag_suffix = ")",           
    theme = theme(
      plot.tag = element_text(size = 40, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 15),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
      plot.tag.position = c(0.02, 0.95),  
      axis.title  = element_text(size = 20),
      axis.text   = element_text(size = 10)
    ) 
  )

final_plot

# Save layout to .png
#ggsave("coy_base_plots.png", final_plot, width = 12, height = 11, dpi = 300)

###########################################################
# KIT FOX BASE MODELS #
###########################################################

#### Setup code for kit fox occupancy models ####

# Read in edited .csv
detHist <- read.csv(file = "kit fox detection hist.csv", 
                    row.names = 1)

# Change from integer to numeric
detHist %>%
  mutate(across(1:51, as.numeric))

# Create detection history 
occu.vuma <- unmarkedFrameOccu(y=detHist, 
                               siteCovs = site.covs.scaled, 
                               obsCovs = obsCovs.scaled)

# ================================
# Kit fox detection variables 
# ================================

# Top detection model
mod1 <- occu( ~ as.factor(cam_moved) + NDVI_1_9km ~ 1, occu.vuma)

#### Camera moved plot ####

newdat <- data.frame(
  cam_moved = factor(c(1, 2)),  # 1=not moved, 2=moved
  NDVI_1_9km = 0
)

preds <- predict(mod1, type = "det", newdata = newdat)
plotdat <- cbind(newdat, preds)

kit.fox.cam <- ggplot(plotdat, aes(x = cam_moved, y = Predicted, group = 1)) +
  geom_line(color = "black", linewidth = 1) +
  geom_point(aes(color = cam_moved), size = 4) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15) +
  scale_color_manual(values = c("#D55E00", "steelblue"),
                     labels = c("Not moved", "Moved")) +
  scale_x_discrete(labels = c("Not moved", "Moved")) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Camera moved?",
       y = "Probability of detection (p)",
       color = "Camera status") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14),
    legend.position = "none"
  )

kit.fox.cam

#### NDVI plot ####

mod1 <- occu(~ cam_moved + NDVI_1_9km ~ 1, occu.vuma)

base_mean <- mean(site.covs$NDVI_1_9km, na.rm = TRUE)
base_sd   <- sd(site.covs$NDVI_1_9km, na.rm = TRUE)

# Generate sequence of unscaled values, then scale them
seq_unscaled <- seq(min(site.covs$NDVI_1_9km, 
                        na.rm = TRUE),
                    max(site.covs$NDVI_1_9km,
                        na.rm = TRUE),
                    length.out = 100)
seq_scaled <- (seq_unscaled - base_mean) / base_sd

# Create plot prediction dataset 
plot_data <- data.frame(
  NDVI_1_9km = seq_scaled,
  NDVI_1_9km_unscaled = seq_unscaled,
  cam_moved = 1
)

# Predict
preds <- predict(mod1, type = "det", newdata = plot_data)
plot_data$Predicted <- preds$Predicted
plot_data$lower <- preds$lower
plot_data$upper <- preds$upper

# Create site_data for observed points
site_data <- site.covs %>%
  dplyr::select(NDVI_1_9km, cam_moved) %>%
  mutate(NDVI_1_9km_scaled = 
           (NDVI_1_9km - base_mean) / base_sd,
         cam_moved = 1
  )

# Predict for observed sites
site_preds <- predict(mod1, type = "det", newdata = data.frame(
  NDVI_1_9km = site_data$NDVI_1_9km_scaled,
  cam_moved = site_data$cam_moved
))

site_data$Predicted <- site_preds$Predicted
site_data$lower <- site_preds$lower
site_data$upper <- site_preds$upper

# Plot 
kit.fox.ndvi  <- ggplot(plot_data, aes(x = NDVI_1_9km_unscaled, 
                                       y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, 
              fill = "darkgray") +
  geom_line(color = "black", linewidth = 0.7) +
  scale_x_continuous(breaks = pretty_breaks(n = 4)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "NDVI within 1.9 km radius",
       y = "Probability of detection (p)") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),   # axis labels bigger
    axis.text  = element_text(size = 14),   # tick labels bigger
    axis.line  = element_line(size = 0.8),  # slightly bolder axis lines
    axis.ticks = element_line(size = 0.8)   # bolder ticks
  )

kit.fox.ndvi

# ================================
# Kit fox habitat selection variables 
# ================================

# Full kit fox base model 
mod1 <- occu( ~  as.factor(cam_moved) + NDVI_1_9km 
              ~ coy_count_avg + jackrabbit_count_avg, occu.vuma, 
                starts = c(-1, -1, -1, -1, -1, -1))

#### Coyote count plot ####

base_mean <- mean(site.covs$coy_count_avg, na.rm = TRUE)
base_sd   <- sd(site.covs$coy_count_avg, na.rm = TRUE)

# Generate sequence of unscaled values, then scale them
seq_unscaled <- seq(min(site.covs$coy_count_avg, 
                        na.rm = TRUE),
                    max(site.covs$coy_count_avg,
                        na.rm = TRUE),
                    length.out = 100)
seq_scaled <- (seq_unscaled - base_mean) / base_sd

# Create plot prediction dataset 
plot_data <- data.frame(
  coy_count_avg = seq_scaled,
  coy_count_avg_unscaled = seq_unscaled,
  jackrabbit_count_avg = 0
)

# Predict
preds <- predict(mod1, type = "state", newdata = plot_data)
plot_data$Predicted <- preds$Predicted
plot_data$lower <- preds$lower
plot_data$upper <- preds$upper

# Create site_data for observed points
site_data <- site.covs %>%
  dplyr::select(coy_count_avg, jackrabbit_count_avg) %>%
  mutate(coy_count_avg_scaled = 
           (coy_count_avg - base_mean) / base_sd,
         jackrabbit_count_avg = 0
  )

# Predict for observed sites
site_preds <- predict(mod1, type = "state", newdata = data.frame(
  coy_count_avg = site_data$coy_count_avg_scaled,
  jackrabbit_count_avg = site_data$jackrabbit_count_avg 
))

site_data$Predicted <- site_preds$Predicted
site_data$lower <- site_preds$lower
site_data$upper <- site_preds$upper

# Plot 
kit.fox.coy  <- ggplot(plot_data, aes(x = coy_count_avg_unscaled, 
                                        y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, 
              fill = "darkgray") +
  geom_line(color = "black", linewidth = 0.7) +
  scale_x_continuous(breaks = pretty_breaks(n = 4)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Index of coyote activity",
       y = "Probability of habitat selection (ψ)") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),   # axis labels bigger
    axis.text  = element_text(size = 14),   # tick labels bigger
    axis.line  = element_line(size = 0.8),  # slightly bolder axis lines
    axis.ticks = element_line(size = 0.8)   # bolder ticks
  )

kit.fox.coy

#### Jackrabbit count plot ####

base_mean <- mean(site.covs$jackrabbit_count_avg, na.rm = TRUE)
base_sd   <- sd(site.covs$jackrabbit_count_avg, na.rm = TRUE)

# Generate sequence of unscaled values, then scale them
seq_unscaled <- seq(min(site.covs$jackrabbit_count_avg, 
                        na.rm = TRUE),
                    max(site.covs$jackrabbit_count_avg,
                        na.rm = TRUE),
                    length.out = 100)
seq_scaled <- (seq_unscaled - base_mean) / base_sd

# Create plot prediction dataset 
plot_data <- data.frame(
  jackrabbit_count_avg = seq_scaled,
  jackrabbit_count_avg_unscaled = seq_unscaled,
  coy_count_avg = 0
)

# Predict
preds <- predict(mod1, type = "state", newdata = plot_data)
plot_data$Predicted <- preds$Predicted
plot_data$lower <- preds$lower
plot_data$upper <- preds$upper

# Create site_data for observed points
site_data <- site.covs %>%
  dplyr::select(coy_count_avg, jackrabbit_count_avg) %>%
  mutate(jackrabbit_count_avg_scaled = 
           (jackrabbit_count_avg - base_mean) / base_sd,
         coy_count_avg = 0
  )

# Predict for observed sites
site_preds <- predict(mod1, type = "state", newdata = data.frame(
  jackrabbit_count_avg = site_data$jackrabbit_count_avg_scaled,
  coy_count_avg = site_data$coy_count_avg 
))

site_data$Predicted <- site_preds$Predicted
site_data$lower <- site_preds$lower
site_data$upper <- site_preds$upper

# Plot 
kit.fox.jack  <- ggplot(plot_data, aes(x = jackrabbit_count_avg_unscaled, 
                                      y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, 
              fill = "darkgray") +
  geom_line(color = "black", linewidth = 0.7) +
  scale_x_continuous(breaks = pretty_breaks(n = 4)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Index of jackrabbit activity",
       y = "Probability of habitat selection (ψ)") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),   # axis labels bigger
    axis.text  = element_text(size = 14),   # tick labels bigger
    axis.line  = element_line(size = 0.8),  # slightly bolder axis lines
    axis.ticks = element_line(size = 0.8)   # bolder ticks
  )

kit.fox.jack

# Create combined plot layout
final_plot <-
  (kit.fox.cam | kit.fox.ndvi ) /
  (kit.fox.coy | kit.fox.jack) +
  plot_layout(widths = c(1,1), heights = c(1,1)) +
  plot_annotation(
    title = "Kit fox",
    tag_levels = "a",          
    tag_prefix = "(",           
    tag_suffix = ")",           
    theme = theme(
      plot.tag = element_text(size = 40, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 15),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
      plot.tag.position = c(0.02, 0.95),  
      axis.title  = element_text(size = 20),
      axis.text   = element_text(size = 10)
    ) 
  )

final_plot

# Save plot layout to .png
#ggsave("kit_fox_base_plots.png", final_plot, width = 12, 
#height = 12, dpi = 300)


###########################################################
# GRAY FOX BASE MODELS #
###########################################################

#### Setup code for gray fox occupancy models ####

# Read in edited .csv
detHist <- read.csv(file = "gray fox detection hist.csv", 
                    row.names = 1)

# Change from integer to numeric
detHist <- detHist %>%
  mutate(across(1:51, as.numeric))

sum(as.matrix(detHist), na.rm = TRUE)

# Create detection history 
occu.urci <- unmarkedFrameOccu(y = detHist, 
                               siteCovs = site.covs.scaled, 
                               obsCovs = obsCovs.scaled)

# ================================
# Gray fox detection variables 
# ================================

# Top detection model
mod1 <- occu( ~ as.factor(cam_moved) ~ 1, occu.urci,starts = c(-1, 0, 0))   

#### Camera moved plot ####

newdat <- data.frame(
  cam_moved = factor(c(1, 2))  # categories
)

preds <- predict(mod1, type = "det", newdata = newdat)
plotdat <- cbind(newdat, preds)

gray.fox.cam <- ggplot(plotdat, aes(x = cam_moved, y = Predicted, group = 1)) +
  geom_line(color = "black", linewidth = 1) +
  geom_point(aes(color = cam_moved), size = 4) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15) +
  scale_color_manual(values = c("#D55E00", "steelblue"),
                     labels = c("Not moved", "Moved")) +
  scale_x_discrete(labels = c("Not moved", "Moved")) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Camera moved?",
       y = "Probability of detection (p)",
       color = "Camera status") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14),
    legend.position = "none"
  )

gray.fox.cam

# ================================
# Gray fox habitat selection variables 
# ================================

# Full gray fox base model 

mod1 <- occu( ~ as.factor(cam_moved) 
              ~ slope + NDVI_1_5km, occu.urci)

#### NDVI plot ####

base_mean <- mean(site.covs$NDVI_1_5km, na.rm = TRUE)
base_sd   <- sd(site.covs$NDVI_1_5km, na.rm = TRUE)

# Generate sequence of unscaled values, then scale them
seq_unscaled <- seq(min(site.covs$NDVI_1_5km, 
                        na.rm = TRUE),
                    max(site.covs$NDVI_1_5km,
                        na.rm = TRUE),
                    length.out = 100)
seq_scaled <- (seq_unscaled - base_mean) / base_sd

# Create plot prediction dataset 
plot_data <- data.frame(
  NDVI_1_5km = seq_scaled,
  NDVI_1_5km_unscaled = seq_unscaled,
  slope = 0
)

# Predict
preds <- predict(mod1, type = "state", newdata = plot_data)
plot_data$Predicted <- preds$Predicted
plot_data$lower <- preds$lower
plot_data$upper <- preds$upper

# Create site_data for observed points
site_data <- site.covs %>%
  dplyr::select(NDVI_1_5km, slope) %>%
  mutate(NDVI_1_5km_scaled = 
           (NDVI_1_5km - base_mean) / base_sd,
         slope = 0
  )

# Predict for observed sites
site_preds <- predict(mod1, type = "state", newdata = data.frame(
  NDVI_1_5km = site_data$NDVI_1_5km_scaled,
  slope = site_data$slope 
))

site_data$Predicted <- site_preds$Predicted
site_data$lower <- site_preds$lower
site_data$upper <- site_preds$upper

# Plot 
gray.fox.ndvi  <- ggplot(plot_data, aes(x = NDVI_1_5km_unscaled, 
                                        y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, 
              fill = "darkgray") +
  geom_line(color = "black", linewidth = 0.7) +
  scale_x_continuous(breaks = pretty_breaks(n = 4)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "NDVI within 1.5 km radius",
       y = "Probability of habitat selection (ψ)") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),   # axis labels bigger
    axis.text  = element_text(size = 14),   # tick labels bigger
    axis.line  = element_line(size = 0.8),  # slightly bolder axis lines
    axis.ticks = element_line(size = 0.8),
    plot.margin = margin(t = 10, r = 10, b = 20, l = 10) 
  )

gray.fox.ndvi

#### Slope plot ####

base_mean <- mean(site.covs$slope, na.rm = TRUE)
base_sd   <- sd(site.covs$slope, na.rm = TRUE)

# Generate sequence of unscaled values, then scale them
seq_unscaled <- seq(min(site.covs$slope, 
                        na.rm = TRUE),
                    max(site.covs$slope,
                        na.rm = TRUE),
                    length.out = 100)
seq_scaled <- (seq_unscaled - base_mean) / base_sd

# Create plot prediction dataset 
plot_data <- data.frame(
  slope = seq_scaled,
  slope_unscaled = seq_unscaled,
  NDVI_1_5km = 0
)

# Predict
preds <- predict(mod1, type = "state", newdata = plot_data)
plot_data$Predicted <- preds$Predicted
plot_data$lower <- preds$lower
plot_data$upper <- preds$upper

# Create site_data for observed points
site_data <- site.covs %>%
  dplyr::select(slope, NDVI_1_5km) %>%
  mutate(slope_scaled = 
           (slope - base_mean) / base_sd,
         NDVI_1_5km = 0
  )

# Predict for observed sites
site_preds <- predict(mod1, type = "state", newdata = data.frame(
  slope = site_data$slope_scaled,
  NDVI_1_5km = site_data$NDVI_1_5km
))

site_data$Predicted <- site_preds$Predicted
site_data$lower <- site_preds$lower
site_data$upper <- site_preds$upper

# Plot 
gray.fox.slope  <- ggplot(plot_data, aes(x = slope_unscaled, 
                                      y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, 
              fill = "darkgray") +
  geom_line(color = "black", linewidth = 0.7) +
  scale_x_continuous(breaks = pretty_breaks(n = 6)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Slope (degrees)",
       y = "Probability of habitat selection (ψ)") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),   # axis labels bigger
    axis.text  = element_text(size = 14),   # tick labels bigger
    axis.line  = element_line(size = 0.8),  # slightly bolder axis lines
    axis.ticks = element_line(size = 0.8)   # bolder ticks
  )

gray.fox.slope

# Create combined plot layout

final_plot <-
  (gray.fox.cam | gray.fox.ndvi | gray.fox.slope ) +
  plot_layout(widths = c(1), heights = c(1,1)) +
  plot_annotation(
    title = "Gray fox",
    tag_levels = "a",    
    tag_prefix = "(",           
    tag_suffix = ")",           
    theme = theme(
      plot.tag = element_text(size = 40, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 15),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
      plot.tag.position = c(0.02, 0.95),  
      axis.title  = element_text(size = 20),
      axis.text   = element_text(size = 10)
    ) 
  )

final_plot

# Save layout to .png
#ggsave("gray_fox_base_plots.png", final_plot, width = 12, 
#height = 10, dpi = 300)

###########################################################
# AMERICAN BADGER BASE MODELS #
###########################################################

#### Setup code for American badger occupancy models ####

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

# ================================
# American badger detection variables 
# ================================

# Top detection model
mod1 <- occu( ~ detection_angle + as.factor(precip.cat) ~ 1, occu.tata)

#### Precipitation plot ####

newdat <- data.frame(
  precip.cat = factor(c(0, 1)),  # categories (0=no, 1=yes),
  detection_angle = 0
)

preds <- predict(mod1, type = "det", newdata = newdat)
plotdat <- cbind(newdat, preds)

badger.precip <- ggplot(plotdat, aes(x = precip.cat, y = Predicted, group = 1)) +
  geom_line(color = "black", linewidth = 1) +
  geom_point(aes(color = precip.cat), size = 4) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15) +
  scale_color_manual(values = c("#D55E00", "steelblue"),
                     labels = c("No", "Yes")) +
  scale_x_discrete(labels = c("No", "Yes")) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Sampling occasion precipitation",
       y = "Probability of detection (p)") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14),
    legend.position = "none"
  )

badger.precip

#### Detection angle #### 

mod1 <- occu( ~ detection_angle ~ 1, occu.tata)

base_mean <- mean(site.covs$detection_angle, na.rm = TRUE)
base_sd   <- sd(site.covs$detection_angle, na.rm = TRUE)

# Generate sequence of unscaled values, then scale them
seq_unscaled <- seq(min(site.covs$detection_angle, 
                        na.rm = TRUE),
                    max(site.covs$detection_angle,
                        na.rm = TRUE),
                    length.out = 100)
seq_scaled <- (seq_unscaled - base_mean) / base_sd

# Create plot prediction dataset 
plot_data <- data.frame(
  detection_angle = seq_scaled,
  detection_angle_unscaled = seq_unscaled
)

# Predict
preds <- predict(mod1, type = "det", newdata = plot_data)
plot_data$Predicted <- preds$Predicted
plot_data$lower <- preds$lower
plot_data$upper <- preds$upper

# Create site_data for observed points
site_data <- site.covs %>%
  dplyr::select(detection_angle) %>%
  mutate(detection_angle_scaled = 
           (detection_angle - base_mean) / base_sd
  )

# Predict for observed sites
site_preds <- predict(mod1, type = "det", newdata = data.frame(
  detection_angle = site_data$detection_angle
))

site_data$Predicted <- site_preds$Predicted
site_data$lower <- site_preds$lower
site_data$upper <- site_preds$upper

# Plot 
badger.angle  <- ggplot(plot_data, aes(x = detection_angle_unscaled, 
                                       y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, 
              fill = "darkgray") +
  geom_line(color = "black", linewidth = 0.7) +
  scale_x_continuous(breaks = pretty_breaks(n = 4)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Detection angle",
       y = "Probability of detection (p)") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),   # axis labels bigger
    axis.text  = element_text(size = 14),   # tick labels bigger
    axis.line  = element_line(size = 0.8),  # slightly bolder axis lines
    axis.ticks = element_line(size = 0.8)   # bolder ticks
  )

badger.angle

# ================================
# American badger habitat selection variables 
# ================================

# Full American badger base model 
mod1 <- occu( ~ detection_angle + as.factor(precip.cat)
              ~  vertical_cover + slope, occu.tata)

#### Slope plot ####

base_mean <- mean(site.covs$slope, na.rm = TRUE)
base_sd   <- sd(site.covs$slope, na.rm = TRUE)

# Generate sequence of unscaled values, then scale them
seq_unscaled <- seq(min(site.covs$slope, 
                        na.rm = TRUE),
                    max(site.covs$slope,
                        na.rm = TRUE),
                    length.out = 100)
seq_scaled <- (seq_unscaled - base_mean) / base_sd

# Create plot prediction dataset 
plot_data <- data.frame(
  slope = seq_scaled,
  slope_unscaled = seq_unscaled,
  vertical_cover = 0
)

# Predict
preds <- predict(mod1, type = "state", newdata = plot_data)
plot_data$Predicted <- preds$Predicted
plot_data$lower <- preds$lower
plot_data$upper <- preds$upper

# Create site_data for observed points
site_data <- site.covs %>%
  dplyr::select(slope, vertical_cover) %>%
  mutate(slope_scaled = 
           (slope - base_mean) / base_sd,
         vertical_cover = 0
  )

# Predict for observed sites
site_preds <- predict(mod1, type = "state", newdata = data.frame(
  slope = site_data$slope_scaled,
  vertical_cover = site_data$vertical_cover
))

site_data$Predicted <- site_preds$Predicted
site_data$lower <- site_preds$lower
site_data$upper <- site_preds$upper

# Plot 
badger.slope  <- ggplot(plot_data, aes(x = slope_unscaled, 
                                         y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, 
              fill = "darkgray") +
  geom_line(color = "black", linewidth = 0.7) +
  scale_x_continuous(breaks = pretty_breaks(n = 6)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Slope (degrees)",
       y = "Probability of habitat selection (ψ)") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),   # axis labels bigger
    axis.text  = element_text(size = 14),   # tick labels bigger
    axis.line  = element_line(size = 0.8),  # slightly bolder axis lines
    axis.ticks = element_line(size = 0.8)   # bolder ticks
  )

badger.slope

#### Vertical cover plot ####

base_mean <- mean(site.covs$vertical_cover, na.rm = TRUE)
base_sd   <- sd(site.covs$vertical_cover, na.rm = TRUE)

# Generate sequence of unscaled values, then scale them
seq_unscaled <- seq(min(site.covs$vertical_cover, 
                        na.rm = TRUE),
                    max(site.covs$vertical_cover,
                        na.rm = TRUE),
                    length.out = 100)
seq_scaled <- (seq_unscaled - base_mean) / base_sd

# Create plot prediction dataset 
plot_data <- data.frame(
  vertical_cover = seq_scaled,
  vertical_cover_unscaled = seq_unscaled,
  slope = 0
)

# Predict
preds <- predict(mod1, type = "state", newdata = plot_data)
plot_data$Predicted <- preds$Predicted
plot_data$lower <- preds$lower
plot_data$upper <- preds$upper

# Create site_data for observed points
site_data <- site.covs %>%
  dplyr::select(slope, vertical_cover) %>%
  mutate(vertical_cover_scaled = 
           (vertical_cover - base_mean) / base_sd,
         slope = 0
  )

# Predict for observed sites
site_preds <- predict(mod1, type = "state", newdata = data.frame(
  vertical_cover = site_data$vertical_cover_scaled,
  slope = site_data$slope
))

site_data$Predicted <- site_preds$Predicted
site_data$lower <- site_preds$lower
site_data$upper <- site_preds$upper

# Plot 
badger.vert.cov  <- ggplot(plot_data, aes(x = vertical_cover_unscaled, 
                                       y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, 
              fill = "darkgray") +
  geom_line(color = "black", linewidth = 0.7) +
  scale_x_continuous(breaks = pretty_breaks(n = 6)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Herbaceous cover height (in)",
       y = "Probability of habitat selection (ψ)") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),   # axis labels bigger
    axis.text  = element_text(size = 14),   # tick labels bigger
    axis.line  = element_line(size = 0.8),  # slightly bolder axis lines
    axis.ticks = element_line(size = 0.8)   # bolder ticks
  )

badger.vert.cov 

# Create combined plot layout
final_plot <-
  (badger.precip | badger.angle ) /
  (badger.vert.cov | badger.slope) +
  plot_layout(widths = c(1,1), heights = c(1,1)) +
  plot_annotation(
    title = "American badger",
    tag_levels = "a",          
    tag_prefix = "(",           
    tag_suffix = ")",           
    theme = theme(
      plot.tag = element_text(size = 40, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 15),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
      plot.tag.position = c(0.02, 0.95),  
      axis.title  = element_text(size = 20),
      axis.text   = element_text(size = 10)
    ) 
  )

final_plot

# Save layout to .png
#ggsave("badger_base_plots.png", final_plot, width = 12, height = 12, dpi = 300)


###########################################################
# STIPED SKUNK BASE MODELS #
###########################################################

#### Setup code for striped skunk occupancy models ####

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

# ================================
# Striped skunk detection variables 
# ================================

# Top detection model
mod1 <-  occu( ~ veg_cover_cam + woodland_percent_1_3km ~ 1, 
                 occu.meme)

#### Visual obstruction transect 1 plot ####

base_mean <- mean(site.covs$veg_cover_cam, na.rm = TRUE)
base_sd   <- sd(site.covs$veg_cover_cam, na.rm = TRUE)

# Generate sequence of unscaled values, then scale them
seq_unscaled <- seq(min(site.covs$veg_cover_cam, 
                        na.rm = TRUE),
                    max(site.covs$veg_cover_cam,
                        na.rm = TRUE),
                    length.out = 100)
seq_scaled <- (seq_unscaled - base_mean) / base_sd

# Create plot prediction dataset 
plot_data <- data.frame(
  veg_cover_cam = seq_scaled,
  veg_cover_cam_unscaled = seq_unscaled,
  woodland_percent_1_3km = 0
)

# Predict
preds <- predict(mod1, type = "det", newdata = plot_data)
plot_data$Predicted <- preds$Predicted
plot_data$lower <- preds$lower
plot_data$upper <- preds$upper

# Create site_data for observed points
site_data <- site.covs %>%
  dplyr::select(veg_cover_cam, woodland_percent_1_3km) %>%
  mutate(veg_cover_cam_scaled = 
           (veg_cover_cam - base_mean) / base_sd,
         woodland_percent_1_3km = 0
  )

# Predict for observed sites
site_preds <- predict(mod1, type = "det", newdata = data.frame(
  veg_cover_cam = site_data$veg_cover_cam_scaled,
  woodland_percent_1_3km = site_data$woodland_percent_1_3km
))

site_data$Predicted <- site_preds$Predicted
site_data$lower <- site_preds$lower
site_data$upper <- site_preds$upper

# Plot 
skunk.vis.obs  <- ggplot(plot_data, aes(x = veg_cover_cam_unscaled, 
                                        y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, 
              fill = "darkgray") +
  geom_line(color = "black", linewidth = 0.7) +
  scale_x_continuous(breaks = pretty_breaks(n = 4)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Veg cover at camera",
       y = "Probability of detection (p)") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),   # axis labels bigger
    axis.text  = element_text(size = 14),   # tick labels bigger
    axis.line  = element_line(size = 0.8),  # slightly bolder axis lines
    axis.ticks = element_line(size = 0.8)   # bolder ticks
  )

skunk.vis.obs

#### Percent woodland plots ####

base_mean <- mean(site.covs$woodland_percent_1_3km, na.rm = TRUE)
base_sd   <- sd(site.covs$woodland_percent_1_3km, na.rm = TRUE)

# Generate sequence of unscaled values, then scale them
seq_unscaled <- seq(min(site.covs$woodland_percent_1_3km, 
                        na.rm = TRUE),
                    max(site.covs$woodland_percent_1_3km,
                        na.rm = TRUE),
                    length.out = 100)
seq_scaled <- (seq_unscaled - base_mean) / base_sd

# Create plot prediction dataset 
plot_data <- data.frame(
  woodland_percent_1_3km = seq_scaled,
  woodland_percent_1_3km_unscaled = seq_unscaled,
  veg_cover_cam = 0
)

# Predict
preds <- predict(mod1, type = "det", newdata = plot_data)
plot_data$Predicted <- preds$Predicted
plot_data$lower <- preds$lower
plot_data$upper <- preds$upper

# Create site_data for observed points
site_data <- site.covs %>%
  dplyr::select(veg_cover_cam, woodland_percent_1_3km) %>%
  mutate(woodland_percent_1_3km_scaled = 
           (woodland_percent_1_3km - base_mean) / base_sd,
         veg_cover_cam = 0
  )

# Predict for observed sites
site_preds <- predict(mod1, type = "det", newdata = data.frame(
  woodland_percent_1_3km = site_data$woodland_percent_1_3km_scaled,
  veg_cover_cam = site_data$veg_cover_cam
))

site_data$Predicted <- site_preds$Predicted
site_data$lower <- site_preds$lower
site_data$upper <- site_preds$upper

# Plot 
skunk.wood  <- ggplot(plot_data, aes(x = woodland_percent_1_3km_unscaled, 
                                        y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, 
              fill = "darkgray") +
  geom_line(color = "black", linewidth = 0.7) +
  scale_x_continuous(breaks = pretty_breaks(n = 4)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Woodland within 1.3 km radius",
       y = "Probability of detection (p)") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),   # axis labels bigger
    axis.text  = element_text(size = 14),   # tick labels bigger
    axis.line  = element_line(size = 0.8),  # slightly bolder axis lines
    axis.ticks = element_line(size = 0.8)   # bolder ticks
  )

skunk.wood

# ================================
# Striped skunk habitat selection variables 
# ================================

# Full striped skunk base model 
mod1 <-  occu( ~ veg_cover_cam + woodland_percent_1_3km 
               ~ herbaceous_cov, occu.meme)

#### Percent grass plot ####

base_mean <- mean(site.covs$herbaceous_cov, na.rm = TRUE)
base_sd   <- sd(site.covs$herbaceous_cov, na.rm = TRUE)

# Generate sequence of unscaled values, then scale them
seq_unscaled <- seq(min(site.covs$herbaceous_cov, 
                        na.rm = TRUE),
                    max(site.covs$herbaceous_cov,
                        na.rm = TRUE),
                    length.out = 100)
seq_scaled <- (seq_unscaled - base_mean) / base_sd

# Create plot prediction dataset 
plot_data <- data.frame(
  herbaceous_cov = seq_scaled,
  herbaceous_cov_unscaled = seq_unscaled
)

# Predict
preds <- predict(mod1, type = "state", newdata = plot_data)
plot_data$Predicted <- preds$Predicted
plot_data$lower <- preds$lower
plot_data$upper <- preds$upper

# Create site_data for observed points
site_data <- site.covs %>%
  dplyr::select(herbaceous_cov) %>%
  mutate(herbaceous_cov_scaled = 
           (herbaceous_cov - base_mean) / base_sd
  )

# Predict for observed sites
site_preds <- predict(mod1, type = "state", newdata = data.frame(
  herbaceous_cov = site_data$herbaceous_cov_scaled
))

site_data$Predicted <- site_preds$Predicted
site_data$lower <- site_preds$lower
site_data$upper <- site_preds$upper

# Plot 
skunk.grass  <- ggplot(plot_data, aes(x = herbaceous_cov_unscaled, 
                                     y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, 
              fill = "darkgray") +
  geom_line(color = "black", linewidth = 0.7) +
  scale_x_continuous(breaks = pretty_breaks(n = 4)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Herbaceous ground cover",
       y = "Probability of habitat selection (ψ)") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),   # axis labels bigger
    axis.text  = element_text(size = 14),   # tick labels bigger
    axis.line  = element_line(size = 0.8),  # slightly bolder axis lines
    axis.ticks = element_line(size = 0.8)   # bolder ticks
  )

skunk.grass

# Create combined layout plot

final_plot <-
  (skunk.vis.obs | skunk.wood | skunk.grass ) +
  plot_layout(widths = c(1), heights = c(1,1)) +
  plot_annotation(
    title = "Striped skunk",
    tag_levels = "a", 
    tag_prefix = "(",           
    tag_suffix = ")",           
    theme = theme(
      plot.tag = element_text(size = 40, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 15),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
      plot.tag.position = c(0.02, 0.95),  
      axis.title  = element_text(size = 20),
      axis.text   = element_text(size = 10)
    ) 
  )

final_plot

# Save layout as .png
#ggsave("skunk_base_plots.png", final_plot, width = 12, height = 10, dpi = 300)


###########################################################
# BLACK-TAILED JACKRABBITT BASE MODELS #
###########################################################

#### Setup code for black-tailed jackrabbit occupancy models ####

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


# ================================
# Black-tailed jackrabbit detection variables 
# ================================

# Top detection model
mod1 <- occu( ~ veg_cover_cam_under_1m ~ 1, occu.leca)

#### Visual obstruction under 1 m plot ####

# Get mean and SD of wind variable
base_mean <- mean(site.covs$veg_cover_cam_under_1m, na.rm = TRUE)
base_sd   <- sd(site.covs$veg_cover_cam_under_1m, na.rm = TRUE)

# Generate sequence of unscaled values, then scale them
seq_unscaled <- seq(min(site.covs$veg_cover_cam_under_1m, 
                        na.rm = TRUE),
                    max(site.covs$veg_cover_cam_under_1m,
                        na.rm = TRUE),
                    length.out = 100)
seq_scaled <- (seq_unscaled - base_mean) / base_sd

# Create plot prediction dataset 
plot_data <- data.frame(
  veg_cover_cam_under_1m = seq_scaled,
  veg_cover_cam_under_1m_unscaled = seq_unscaled
)

# Predict
preds <- predict(mod1, type = "det", newdata = plot_data)
plot_data$Predicted <- preds$Predicted
plot_data$lower <- preds$lower
plot_data$upper <- preds$upper

# Create site_data for observed points
site_data <- site.covs %>%
  dplyr::select(veg_cover_cam_under_1m) %>%
  mutate(veg_cover_cam_under_1m_scaled = 
           (veg_cover_cam_under_1m - base_mean) / base_sd
  )

# Predict for observed sites
site_preds <- predict(mod1, type = "det", newdata = data.frame(
  veg_cover_cam_under_1m = 
    site_data$veg_cover_cam_under_1m_scaled
))

site_data$Predicted <- site_preds$Predicted
site_data$lower <- site_preds$lower
site_data$upper <- site_preds$upper

# Plot 
jack.vis.obs  <- ggplot(plot_data, aes(x = 
                                   veg_cover_cam_under_1m_unscaled, 
                                       y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, 
              fill = "darkgray") +
  geom_line(color = "black", linewidth = 0.7) +
  scale_x_continuous(breaks = pretty_breaks(n = 6)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Veg cover under 1 m",
       y = "Probability of detection (p)") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),   # axis labels bigger
    axis.text  = element_text(size = 14),   # tick labels bigger
    axis.line  = element_line(size = 0.8),  # slightly bolder axis lines
    axis.ticks = element_line(size = 0.8)   # bolder ticks
  )

jack.vis.obs

#### Precipitation ####

mod1 <- occu( ~ veg_cover_cam_under_1m + as.factor(precip.cat)
              ~ 1, occu.leca)

newdat <- data.frame(
  precip.cat = factor(c(0, 1)),  # categories (0=no, 1=yes),
  veg_cover_cam_under_1m  = 0
)

preds <- predict(mod1, type = "det", newdata = newdat)
plotdat <- cbind(newdat, preds)

jack.precip <- ggplot(plotdat, aes(x = precip.cat, y = Predicted, group = 1)) +
  geom_line(color = "black", linewidth = 1) +
  geom_point(aes(color = precip.cat), size = 4) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15) +
  scale_color_manual(values = c("#D55E00", "steelblue"),
                     labels = c("No", "Yes")) +
  scale_x_discrete(labels = c("No", "Yes")) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Sampling occasion precipitation",
       y = "Probability of detection (p)") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14),
    legend.position = "none"
  )

jack.precip

# ================================
# Black-tailed jackrabbit habitat selection variables 
# ================================

# Full black-tailed jackrabbit base model 
mod1 <- occu( ~  as.factor(precip.cat) + veg_cover_cam_under_1m
              ~  woodland_percent_0_9km + slope, occu.leca)

#### Slope plot ####

# Get mean and SD of wind variable
base_mean <- mean(site.covs$slope, na.rm = TRUE)
base_sd   <- sd(site.covs$slope, na.rm = TRUE)

# Generate sequence of unscaled values, then scale them
seq_unscaled <- seq(min(site.covs$slope, 
                        na.rm = TRUE),
                    max(site.covs$slope,
                        na.rm = TRUE),
                    length.out = 100)
seq_scaled <- (seq_unscaled - base_mean) / base_sd

# Create plot prediction dataset 
plot_data <- data.frame(
  slope = seq_scaled,
  slope_unscaled = seq_unscaled,
  woodland_percent_0_9km = 0
)

# Predict
preds <- predict(mod1, type = "state", newdata = plot_data)
plot_data$Predicted <- preds$Predicted
plot_data$lower <- preds$lower
plot_data$upper <- preds$upper

# Create site_data for observed points
site_data <- site.covs %>%
  dplyr::select(slope, woodland_percent_0_9km) %>%
  mutate(slope_scaled = 
           (slope - base_mean) / base_sd,
         woodland_percent_0_9km = 0
  )

# Predict for observed sites
site_preds <- predict(mod1, type = "state", newdata = data.frame(
  slope = site_data$slope_scaled,
  woodland_percent_0_9km = site_data$woodland_percent_0_9km
))

site_data$Predicted <- site_preds$Predicted
site_data$lower <- site_preds$lower
site_data$upper <- site_preds$upper

# Plot 
jack.slope  <- ggplot(plot_data, aes(x = slope_unscaled, 
                                     y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, 
              fill = "darkgray") +
  geom_line(color = "black", linewidth = 0.7) +
  scale_x_continuous(breaks = pretty_breaks(n = 6)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Slope",
       y = "Probability of site use (ψ)") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),   # axis labels bigger
    axis.text  = element_text(size = 14),   # tick labels bigger
    axis.line  = element_line(size = 0.8),  # slightly bolder axis lines
    axis.ticks = element_line(size = 0.8)   # bolder ticks
  )

jack.slope

#### Woodland proportion plot ####

# Get mean and SD of wind variable
base_mean <- mean(site.covs$woodland_percent_0_9km, na.rm = TRUE)
base_sd   <- sd(site.covs$woodland_percent_0_9km, na.rm = TRUE)

# Generate sequence of unscaled values, then scale them
seq_unscaled <- seq(min(site.covs$woodland_percent_0_9km, 
                        na.rm = TRUE),
                    max(site.covs$woodland_percent_0_9km,
                        na.rm = TRUE),
                    length.out = 100)
seq_scaled <- (seq_unscaled - base_mean) / base_sd

# Create plot prediction dataset 
plot_data <- data.frame(
  woodland_percent_0_9km = seq_scaled,
  woodland_percent_0_9km_unscaled = seq_unscaled,
  slope = 0
)

# Predict
preds <- predict(mod1, type = "state", newdata = plot_data)
plot_data$Predicted <- preds$Predicted
plot_data$lower <- preds$lower
plot_data$upper <- preds$upper

# Create site_data for observed points
site_data <- site.covs %>%
  dplyr::select(slope, woodland_percent_0_9km) %>%
  mutate(woodland_percent_0_9km_scaled = 
           (woodland_percent_0_9km - base_mean) / base_sd,
         slope = 0
  )

# Predict for observed sites
site_preds <- predict(mod1, type = "state", newdata = data.frame(
  woodland_percent_0_9km = site_data$woodland_percent_0_9km_scaled,
  slope = site_data$slope
))

site_data$Predicted <- site_preds$Predicted
site_data$lower <- site_preds$lower
site_data$upper <- site_preds$upper

# Plot 
jack.wood  <- ggplot(plot_data, aes(x = woodland_percent_0_9km_unscaled, 
                                    y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, 
              fill = "darkgray") +
  geom_line(color = "black", linewidth = 0.7) +
  scale_x_continuous(breaks = pretty_breaks(n = 6)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Woodland within 0.9 km",
       y = "Probability of habitat selection (ψ)") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),   # axis labels bigger
    axis.text  = element_text(size = 14),   # tick labels bigger
    axis.line  = element_line(size = 0.8),  # slightly bolder axis lines
    axis.ticks = element_line(size = 0.8),   # bolder ticks
    plot.margin = margin(t = 10, r = 10, b = 20, l = 10)
  )

jack.wood

# Create combined layout plot
final_plot <-
  (jack.precip | jack.vis.obs ) /
  (jack.wood   | jack.slope ) +
  plot_layout(widths = c(1,1), heights = c(1,1)) +
  plot_annotation(
    title = "Black-tailed jackrabbit",
    tag_levels = "a",          
    tag_prefix = "(",           
    tag_suffix = ")",           
    theme = theme(
      plot.tag = element_text(size = 40, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 15),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
      plot.tag.position = c(0.02, 0.95),  
      axis.title  = element_text(size = 20),
      axis.text   = element_text(size = 10)
    ) 
  )

final_plot

# Save layout to .png
#ggsave("jack_base_plots.png", final_plot, width = 12, height = 12, dpi = 300)


###########################################################
# DESERT COTTONTAIL BASE MODELS #
###########################################################

#### Setup code for desert cottontail occupancy models ####

# Read in edited .csv
detHist <- read.csv(file = "cottontail_detection_hist.csv", 
                    row.names = 1)

# Change from integer to numeric
detHist %>%
  mutate(across(1:51, as.numeric))

# Create detection history 
occu.syau <- unmarkedFrameOccu(y=detHist, 
                               siteCovs = site.covs.scaled, 
                               obsCovs = obsCovs.scaled)

# ================================
# Desert cottontail detection variables 
# ================================

# Top detection model
mod1 <- occu( ~  veg_cover_cam_under_1m + biotic_com_2
              ~  1, occu.syau)

#### Visual obstruction under 1m plot ####

# Fit model using scaled wind variable
mod1 <- occu( ~ veg_cover_cam_under_1m ~ 1, occu.syau)

# Get mean and SD of wind variable
base_mean <- mean(site.covs$veg_cover_cam_under_1m, na.rm = TRUE)
base_sd   <- sd(site.covs$veg_cover_cam_under_1m, na.rm = TRUE)

# Generate sequence of unscaled values, then scale them
seq_unscaled <- seq(min(site.covs$veg_cover_cam_under_1m, 
                                na.rm = TRUE),
                            max(site.covs$veg_cover_cam_under_1m,
                                na.rm = TRUE),
                            length.out = 100)
seq_scaled <- (seq_unscaled - base_mean) / base_sd

# Create plot prediction dataset 
plot_data <- data.frame(
  veg_cover_cam_under_1m = seq_scaled,
  veg_cover_cam_under_1m_unscaled = seq_unscaled
)

# Predict
preds <- predict(mod1, type = "det", newdata = plot_data)
plot_data$Predicted <- preds$Predicted
plot_data$lower <- preds$lower
plot_data$upper <- preds$upper

# Create site_data for observed points
site_data <- site.covs %>%
  dplyr::select(veg_cover_cam_under_1m) %>%
  mutate(veg_cover_cam_under_1m_scaled = 
           (veg_cover_cam_under_1m - base_mean) / base_sd
  )

# Predict for observed sites
site_preds <- predict(mod1, type = "det", newdata = data.frame(
  veg_cover_cam_under_1m = site_data$veg_cover_cam_under_1m_scaled
))

site_data$Predicted <- site_preds$Predicted
site_data$lower <- site_preds$lower
site_data$upper <- site_preds$upper

# Plot 
cotton.vis.obs  <- ggplot(plot_data, aes(x = veg_cover_cam_under_1m_unscaled, 
                                  y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, 
              fill = "darkgray") +
  geom_line(color = "black", linewidth = 0.7) +
  scale_x_continuous(breaks = pretty_breaks(n = 4)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Vegetation cover under 1 m",
       y = "Probability of detection (p)") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),   # axis labels bigger
    axis.text  = element_text(size = 14),   # tick labels bigger
    axis.line  = element_line(size = 0.8),  # slightly bolder axis lines
    axis.ticks = element_line(size = 0.8)   # bolder ticks
  )

cotton.vis.obs 

#### biotic community plot ####

mod1 <- occu( ~  veg_cover_cam_under_1m + biotic_com_2
              ~  1, occu.syau)

# New data for predictions
newdat <- data.frame(
  veg_cover_cam_under_1m = 0,
  biotic_com_2 = factor(c(1, 2))  # 1=Grassland, 2=Woodland
)

# Predict detection probabilities
preds <- predict(mod1, type = "det", newdata = newdat)
plotdat <- cbind(newdat, preds)

# Rename factor levels for clarity
plotdat$biotic_com_2 <- factor(plotdat$biotic_com_2,
                               levels = c(1, 2),
                               labels = c("Grassland", "Woodland"))

cotton.biotic <- ggplot(plotdat, aes(x = biotic_com_2, y = Predicted, group = 1)) +
  geom_line(color = "black", linewidth = 0.8) +  # line connecting the two points
  geom_point(aes(color = biotic_com_2), size = 5) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.2, linewidth = 0.9, color = "black") +  # error bars black
  scale_color_manual(values = c("Grassland" = "darkorange", 
                                "Woodland" = "steelblue")) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Biotic community",
       y = "Probability of detection (p)") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14),
    legend.position = "none"
  )

cotton.biotic

# ================================
# Desert cottontail habitat selection variables 
# ================================

# Full desert cottontail base model 

mod1 <- occu( ~  veg_cover_cam_under_1m + 
                 as.factor(biotic_com_2) 
              ~  NDVI_0_1km, occu.syau)

#### NDVI plot ####

# Get mean and SD of wind variable
base_mean <- mean(site.covs$NDVI_0_1km, na.rm = TRUE)
base_sd   <- sd(site.covs$NDVI_0_1km, na.rm = TRUE)

# Generate sequence of unscaled values, then scale them
seq_unscaled <- seq(min(site.covs$NDVI_0_1km, 
                        na.rm = TRUE),
                    max(site.covs$NDVI_0_1km,
                        na.rm = TRUE),
                    length.out = 100)
seq_scaled <- (seq_unscaled - base_mean) / base_sd

# Create plot prediction dataset 
plot_data <- data.frame(
  NDVI_0_1km = seq_scaled,
  NDVI_0_1km_unscaled = seq_unscaled
)

# Predict
preds <- predict(mod1, type = "state", newdata = plot_data)
plot_data$Predicted <- preds$Predicted
plot_data$lower <- preds$lower
plot_data$upper <- preds$upper

# Create site_data for observed points
site_data <- site.covs %>%
  dplyr::select(NDVI_0_1km) %>%
  mutate(NDVI_0_1km_scaled = 
           (NDVI_0_1km - base_mean) / base_sd
  )

# Predict for observed sites
site_preds <- predict(mod1, type = "state", newdata = data.frame(
  NDVI_0_1km = site_data$NDVI_0_1km_scaled
))

site_data$Predicted <- site_preds$Predicted
site_data$lower <- site_preds$lower
site_data$upper <- site_preds$upper

# Plot 
cotton.ndvi  <- ggplot(plot_data, aes(x = NDVI_0_1km_unscaled, 
                                         y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, 
              fill = "darkgray") +
  geom_line(color = "black", linewidth = 0.7) +
  scale_x_continuous(breaks = pretty_breaks(n = 4)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "NDVI 0.1 km radius",
       y = "Probability of habitat selection (ψ)") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),   # axis labels bigger
    axis.text  = element_text(size = 14),   # tick labels bigger
    axis.line  = element_line(size = 0.8),  # slightly bolder axis lines
    axis.ticks = element_line(size = 0.8)   # bolder ticks
  )

cotton.ndvi

# Create combined layout plot

final_plot <-
  (cotton.biotic | cotton.vis.obs | cotton.ndvi) +
  plot_layout(widths = c(1), heights = c(1)) +
  plot_annotation(
    title = "Desert cottontail",
    tag_levels = "a",          
    tag_prefix = "(",           
    tag_suffix = ")",           
    theme = theme(
      plot.tag = element_text(size = 40, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 15),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
      plot.tag.position = c(0.02, 0.95),  
      axis.title  = element_text(size = 20),
      axis.text   = element_text(size = 10)
    ) 
  )

final_plot

# Save layout to .png
#ggsave("cottontail_base_plots.png", final_plot, width = 12, 
       #height = 10, dpi = 300)


########################################################
######################### END ##########################
########################################################