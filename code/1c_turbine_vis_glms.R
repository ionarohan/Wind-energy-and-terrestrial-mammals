###########################################################
## Code for: Terrestrial mammals alter habitat selection ##
###### in response to turbine visibility at a central #####
############## New Mexico wind energy facility ############
###########################################################
### Script for the creation of generalized linear mixed ###
##### models to determine which variables best describe ###
#################### turbine visibility ###################
###########################################################
################ Script by: Iona Rohan ####################
############ Contact: ionarohan12@gmail.com ###############
###########################################################
########## Date Last Modified: 23-March-2026 ##############
###########################################################

#Clear work environment
rm(list=ls())

#Note: If you opened this script through the .Rproj file, the only line you 
  #should need to change for the script to run (assuming packages are installed)   #is     the homewd directory on line 24.

#Set home working directory
# e.g. homewd = "C:/Users/ionar/Desktop/R Repository/Wind-energy-and-terrestrial-mammals/"
homewd = "<insert your folder here and end with a forward slash>"

#Set wd to data folder on your local computer 
setwd(paste0(homewd, "data/"))

# Load packages 
library(ggplot2)

#### Model Setup ####

# Read in unscaled site-level covariates
site.covs <- read.csv("site_covs.csv", nrows = 102, header = TRUE)

## Results from these models are found in supporting information S3, 
## Tables S3.1 - S3.12

########################################################
############ RUN MODELS FOR 150 CM TURBINE VIS #########
########################################################

# Veg concealment cover 

model <- glm(X150cm_turbine_vis ~ veg_cover,
             data = site.covs,
             family = gaussian)

summary(model)

# Plot effect
ggplot(site.covs, aes(x = veg_cover, y = X150cm_turbine_vis)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "glm", method.args = list(family = "gaussian"), 
              se = TRUE) +
  theme_classic()

# Tree density 

model <- glm(X150cm_turbine_vis ~ tree_density_5_70,
             data = site.covs,
             family = gaussian)

summary(model)

# Plot effect
ggplot(site.covs, aes(x = tree_density_5_70, y = X150cm_turbine_vis)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "glm", method.args = list(family = "gaussian"), 
              se = TRUE) +
  theme_classic()

# Vertical herbaceous cover 

model <- glm(X150cm_turbine_vis ~ vertical_cover,
             data = site.covs,
             family = gaussian)

summary(model)

# Plot effect
ggplot(site.covs, aes(x = vertical_cover, y = X150cm_turbine_vis)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "glm", method.args = list(family = "gaussian"), 
              se = TRUE) +
  theme_classic()

# Shrub density 

model <- glm(X150cm_turbine_vis ~ shrub_yucca_density,
             data = site.covs,
             family = gaussian)

summary(model)

# Plot effect
ggplot(site.covs, aes(x = shrub_yucca_density, y = X150cm_turbine_vis)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "glm", method.args = list(family = "gaussian"), 
              se = TRUE) +
  theme_classic()

# Topographic position 

model <- glm(X150cm_turbine_vis ~ topo_pos,
             data = site.covs,
             family = gaussian)

summary(model)

# Plot effect
ggplot(site.covs, aes(x = topo_pos, y = X150cm_turbine_vis)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "glm", method.args = list(family = "gaussian"), 
              se = TRUE) +
  theme_classic()

# Slope

model <- glm(X150cm_turbine_vis ~ slope,
             data = site.covs,
             family = gaussian)

summary(model)

# Plot effect
ggplot(site.covs, aes(x = slope, y = X150cm_turbine_vis)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "glm", method.args = list(family = "gaussian"), 
              se = TRUE) +
  theme_classic()


########################################################
############ RUN MODELS FOR 50 CM TURBINE VIS ##########
########################################################

# Veg concealment cover 

model <- glm(X50cm_turbine_vis ~ veg_cover,
             data = site.covs,
             family = gaussian)

summary(model)

# Plot effect
ggplot(site.covs, aes(x = veg_cover, y = X50cm_turbine_vis)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "glm", method.args = list(family = "gaussian"), 
              se = TRUE) +
  theme_classic()

# Tree density 

model <- glm(X50cm_turbine_vis ~ tree_density_5_70,
             data = site.covs,
             family = gaussian)

summary(model)

# Plot effect
ggplot(site.covs, aes(x = tree_density_5_70, y = X50cm_turbine_vis)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "glm", method.args = list(family = "gaussian"), 
              se = TRUE) +
  theme_classic()

# Vertical herbaceous cover 

model <- glm(X50cm_turbine_vis ~ vertical_cover,
             data = site.covs,
             family = gaussian)

summary(model)

# Plot effect
ggplot(site.covs, aes(x = vertical_cover, y = X50cm_turbine_vis)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "glm", method.args = list(family = "gaussian"), 
              se = TRUE) +
  theme_classic()

# Shrub density 

model <- glm(X50cm_turbine_vis ~ shrub_yucca_density,
             data = site.covs,
             family = gaussian)

summary(model)

# Plot effect
ggplot(site.covs, aes(x = shrub_yucca_density, y = X50cm_turbine_vis)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "glm", method.args = list(family = "gaussian"), 
              se = TRUE) +
  theme_classic()

# Topographic position 

model <- glm(X50cm_turbine_vis ~ topo_pos,
             data = site.covs,
             family = gaussian)

summary(model)

# Plot effect
ggplot(site.covs, aes(x = topo_pos, y = X50cm_turbine_vis)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "glm", method.args = list(family = "gaussian"), 
              se = TRUE) +
  theme_classic()

# Slope

model <- glm(X150cm_turbine_vis ~ slope,
             data = site.covs,
             family = gaussian)

summary(model)

# Plot effect
ggplot(site.covs, aes(x = slope, y = X50cm_turbine_vis)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "glm", method.args = list(family = "gaussian"), 
              se = TRUE) +
  theme_classic()
