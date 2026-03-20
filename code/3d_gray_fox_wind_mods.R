###########################################################
########### Code for: wind models for gray foxes ##########
###########################################################
################ Script by: Iona Rohan ####################
############ Contact: ionarohan12@gmail.com ###############
###########################################################
########## Date Last Modified: 19-March-2026 ##############
###########################################################

#Clear work environment
rm(list=ls())

#Note: If you opened this script through the .Rproj file, the only line you 
  #should need to change for the script to run (assuming packages are installed) 
  #is the homewd directory on line 19.

#Set home working directory
  #e.g. homewd = "C:/Users/ionar/Desktop/R Repository/Wind-energy-and-terrestrial-mammals/"
homewd = "<insert your folder here and end with a forward slash>"

#Set wd to data folder on your local computer 
setwd(paste0(homewd, "data/"))

# Load packages 
library(unmarked)
library(tidyverse)
library(MuMIn)

###########################################################
# SETUP CODE FOR GRAY FOX OCCUPANCY MODELS #
###########################################################

# Read in edited .csv
detHist <- read.csv(file = "gray fox detection hist.csv", row.names = 1)
# Read in site.covs.scaled
site.covs.scaled <- readRDS("site_covs_scaled.RData")
# Read in obsCovs.scaled
obsCovs.scaled <- readRDS("obsCovs_scaled.RData")

# Change from integer to numeric
detHist <- detHist %>%
  mutate(across(1:51, as.numeric))

# Create detection history 
occu.urci <- unmarkedFrameOccu(y = detHist, 
                               siteCovs = site.covs.scaled, 
                               obsCovs = obsCovs.scaled)

#################################################
# CREATE THE WIND MODELS #
#################################################

#### Models with no change over all turbine variables ####

Habitat <- occu( ~ as.factor(cam_moved) ~ slope + NDVI_1_5km, occu.urci)

Null <- occu( ~ as.factor(cam_moved) ~ 1, occu.urci)

Canopy_Cov <- occu( ~ as.factor(cam_moved) ~ canopy_cov, occu.urci)

Canopy_Cov_Slope <- occu( ~ as.factor(cam_moved) 
                          ~ canopy_cov + slope, occu.urci)

NDVI_Occu <- occu( ~ as.factor(cam_moved) ~ NDVI_1_5km, occu.urci)

Slope_Occu <- occu( ~ as.factor(cam_moved) ~ slope, occu.urci)

Bio_Com_Slope <- occu( ~ as.factor(cam_moved) 
                       ~ slope + as.factor(biotic_com_2), occu.urci)

Bio_Com <- occu( ~ as.factor(cam_moved) ~ as.factor(biotic_com_2), occu.urci)

## The table with the comparison of the relative weight of evidence between
  #occupancy models separately ranked for the effect of each wind energy 
  #variable on the probability of habitat selection (ψ) of gray foxes is found 
  #in Table S3.16.

## Check correlations between wind variables and habitat variables 

# Read in site-level covariates
site.covs <- read.csv("site_covs.csv", nrows = 102, header = TRUE)
# site-level variable correlations
wind.hab.cor <- site.covs %>% 
  select(NDVI_1_5km, biotic_com_2, slope, canopy_cov,
         turbine_interior, X50cm_turbine_vis, turbine_dist,
         turbine_density_1_5km, turbine_rd_dist, turbine_rd_density_1_5km)

cors <- cor(wind.hab.cor, method='spearman')  
# none correlated above |0.7|

#### Turbine Interior Models ####

Turbine_Int_Occu <- occu( ~ as.factor(cam_moved) 
                          ~ as.factor(turbine_interior), occu.urci)

Habitat_Add_Turbine_Int <- occu( ~ as.factor(cam_moved) 
                                 ~ as.factor(turbine_interior) +
                                   NDVI_1_5km + slope, occu.urci)

NDVI_Add_Turbine_Int <- occu( ~ as.factor(cam_moved) 
                              ~ as.factor(turbine_interior) +
                                NDVI_1_5km, occu.urci)

Slope_Add_Turbine_Int <- occu( ~ as.factor(cam_moved) 
                               ~ as.factor(turbine_interior) +
                                 slope, occu.urci)

Bio_Com_X_Turbine_Int_Slope <- occu( ~ as.factor(cam_moved) 
                                     ~ as.factor(turbine_interior) *
                                       as.factor(biotic_com_2) + slope, 
                                     occu.urci,
                                     starts = c(0, -3, -1, -1, 1, 0.5, 2))

NDVI_X_Turbine_Int_Slope <- occu( ~ as.factor(cam_moved) 
                                  ~ as.factor(turbine_interior) * 
                                    NDVI_1_5km + slope, occu.urci)

# 85% CI for bio community interaction 
# Extract coefficients and VCOV matrix
coef_est <- coef(Bio_Com_X_Turbine_Int_Slope, type = "state")
vcov_mat <- vcov(Bio_Com_X_Turbine_Int_Slope, type = "state")

# Exact coefficient names from the model
intercept_name   <- "psi(Int)"
turbine_name     <- "psi(as.factor(turbine_interior)1)"
habitat_name     <- "psi(as.factor(biotic_com_2)2)"
interaction_name <- 
  "psi(as.factor(turbine_interior)1:as.factor(biotic_com_2)2)"
slope_name       <- "psi(slope)"

# Grassland, outside wind farm (reference) 
beta_grassland <- coef_est[intercept_name]
se_grassland <- sqrt(vcov_mat[intercept_name, intercept_name])

# Grassland, inside wind farm
beta_grassland_inside <- coef_est[intercept_name] + coef_est[turbine_name]
var_grassland_inside <- vcov_mat[intercept_name, intercept_name] +
  vcov_mat[turbine_name, turbine_name] +
  2 * vcov_mat[intercept_name, turbine_name]
se_grassland_inside <- sqrt(var_grassland_inside)

# Woodland, outside wind farm 
beta_woodland <- coef_est[intercept_name] + coef_est[habitat_name]
var_woodland <- vcov_mat[intercept_name, intercept_name] +
  vcov_mat[habitat_name, habitat_name] +
  2 * vcov_mat[intercept_name, habitat_name]
se_woodland <- sqrt(var_woodland)

# Woodland, inside wind farm 
beta_woodland_inside <- coef_est[intercept_name] + coef_est[turbine_name] +
  coef_est[habitat_name] + coef_est[interaction_name]

var_woodland_inside <- vcov_mat[intercept_name, intercept_name] +
  vcov_mat[turbine_name, turbine_name] +
  vcov_mat[habitat_name, habitat_name] +
  vcov_mat[interaction_name, interaction_name] +
  2 * vcov_mat[intercept_name, turbine_name] +
  2 * vcov_mat[intercept_name, habitat_name] +
  2 * vcov_mat[intercept_name, interaction_name] +
  2 * vcov_mat[turbine_name, habitat_name] +
  2 * vcov_mat[turbine_name, interaction_name] +
  2 * vcov_mat[habitat_name, interaction_name]

se_woodland_inside <- sqrt(var_woodland_inside)

# z-scores for 85% CI
z <- qnorm(c(0.075, 0.925))

# Confidence intervals 
ci_grassland_out <- beta_grassland + z * se_grassland
ci_grassland_in  <- beta_grassland_inside + z * se_grassland_inside
ci_woodland_out  <- beta_woodland + z * se_woodland
ci_woodland_in   <- beta_woodland_inside + z * se_woodland_inside

# Tidy output 
ci_table <- data.frame(
  Habitat = c("Grassland", "Grassland", "Woodland", "Woodland"),
  Turbine_Location = c("Outside", "Inside", "Outside", "Inside"),
  Estimate = c(beta_grassland, beta_grassland_inside,
               beta_woodland, beta_woodland_inside),
  CI_85_lower = c(ci_grassland_out[1], ci_grassland_in[1],
                  ci_woodland_out[1], ci_woodland_in[1]),
  CI_85_upper = c(ci_grassland_out[2], ci_grassland_in[2],
                  ci_woodland_out[2], ci_woodland_in[2])
)

ci_table

# Calculate the 85% CI for the NDVI interaction 
# Extract coefficients and VCOV matrix
coef_est <- coef(NDVI_X_Turbine_Int_Slope, type = "state")
vcov_mat <- vcov(NDVI_X_Turbine_Int_Slope, type = "state")

#  NDVI effect OUTSIDE wind farm (reference = 0) 
beta_outside <- coef_est["psi(NDVI_1_5km)"]

se_outside <- sqrt(
  vcov_mat["psi(NDVI_1_5km)", "psi(NDVI_1_5km)"]
)

# NDVI effect INSIDE wind farm (reference + interaction)
beta_inside <- coef_est["psi(NDVI_1_5km)"] +
  coef_est["psi(as.factor(turbine_interior)1:NDVI_1_5km)"]

var_inside <-
  vcov_mat["psi(NDVI_1_5km)", "psi(NDVI_1_5km)"] +
  vcov_mat["psi(as.factor(turbine_interior)1:NDVI_1_5km)",
           "psi(as.factor(turbine_interior)1:NDVI_1_5km)"] +
  2 * vcov_mat["psi(NDVI_1_5km)",
               "psi(as.factor(turbine_interior)1:NDVI_1_5km)"]

se_inside <- sqrt(var_inside)

# z-scores for 85% CI 
z <- qnorm(c(0.075, 0.925))

# Confidence intervals
ci_outside <- beta_outside + z * se_outside
ci_inside  <- beta_inside  + z * se_inside

# Tidy output 
ci_table <- data.frame(
  Turbine_Location = c("Outside wind farm", "Inside wind farm"),
  Estimate = c(beta_outside, beta_inside),
  CI_85_lower = c(ci_outside[1], ci_inside[1]),
  CI_85_upper = c(ci_outside[2], ci_inside[2])
)

ci_table

# AICc Table with all possible turbine interior models 

model_list <- list(
  "NDVI + slope" = Habitat,
  "Null" = Null,
  "Slope" = Slope_Occu,
  "NDVI" = NDVI_Occu,
  "Biotic community + slope" = Bio_Com_Slope,
  "Biotic community" = Bio_Com,
  "Turbine interior" = Turbine_Int_Occu,
  "Turbine interior + NDVI + slope" = Habitat_Add_Turbine_Int,
  "Turbine interior + NDVI" = NDVI_Add_Turbine_Int,
  "Turbine interior + slope" = Slope_Add_Turbine_Int,
  "Turbine interior x NDVI + slope" = NDVI_X_Turbine_Int_Slope,
  "Turbine interior x biotic community + slope" = Bio_Com_X_Turbine_Int_Slope
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine interior" = "Null",
  "Turbine interior + NDVI + slope" = "NDVI + slope",
  "Turbine interior + NDVI" = "NDVI",
  "Turbine interior + slope" = "Slope",
  "Turbine interior x NDVI + slope" = "NDVI + slope",
  "Turbine interior x biotic community + slope" = "Biotic community + slope"
)

# Ensure model_list matches wind_pairs
stopifnot(all(names(wind_pairs) %in% names(model_list)))
stopifnot(all(unlist(wind_pairs) %in% names(model_list)))

# The beta values (ß), standard errors, and 85% confidence intervals for 
  #parameter estimates within the best-supported model describing the effect of
  #the interior of the wind farm on the probability of habitat selection (ψ) for 
  #gray foxes are listed in Table S3.25.

#### Turbine Visibility Models ####

Turbine_Vis_Occu <- occu( ~ as.factor(cam_moved) ~ X50cm_turbine_vis, occu.urci)

Habitat_Add_Turbine_Vis <- occu( ~ as.factor(cam_moved) 
                                 ~ X50cm_turbine_vis + NDVI_1_5km + slope, 
                                 occu.urci)

NDVI_Add_Turbine_Vis <- occu( ~ as.factor(cam_moved) 
                              ~ X50cm_turbine_vis + NDVI_1_5km, occu.urci)

Slope_Add_Turbine_Vis <- occu( ~ as.factor(cam_moved)
                               ~ X50cm_turbine_vis + slope, occu.urci)

Bio_Com_X_Turbine_Vis_Slope <- occu( ~ as.factor(cam_moved) 
                                     ~ X50cm_turbine_vis * 
                                       as.factor(biotic_com_2) + slope, 
                                     occu.urci,
                                     starts = c(0, -3, 0, 1, -1, 0, -3))

Canopy_X_Turbine_Vis_Slope <- occu( ~ as.factor(cam_moved) ~ X50cm_turbine_vis *
                                      canopy_cov + slope, occu.urci,
                                    starts = c(-2, -3, -5, -1, -1, 0, -3))

NDVI_X_Turbine_Vis_Slope <- occu( ~ as.factor(cam_moved) ~ X50cm_turbine_vis *
                                    NDVI_1_5km + slope, occu.urci)

# 85% CI for bio community interaction                             
# Extract coefficients and VCOV matrix
coef_est <- coef(Bio_Com_X_Turbine_Vis_Slope, type = "state")
vcov_mat <- vcov(Bio_Com_X_Turbine_Vis_Slope, type = "state")

beta_grassland <- coef_est["psi(X50cm_turbine_vis)"]
se_grassland <- sqrt(vcov_mat["psi(X50cm_turbine_vis)", 
                              "psi(X50cm_turbine_vis)"])

beta_woodland <- coef_est["psi(X50cm_turbine_vis)"] + 
  coef_est["psi(X50cm_turbine_vis:as.factor(biotic_com_2)2)"]

var_woodland <- 
  vcov_mat["psi(X50cm_turbine_vis)", "psi(X50cm_turbine_vis)"] +
  vcov_mat["psi(X50cm_turbine_vis:as.factor(biotic_com_2)2)",
           "psi(X50cm_turbine_vis:as.factor(biotic_com_2)2)"] +
  2 * vcov_mat["psi(X50cm_turbine_vis)", 
               "psi(X50cm_turbine_vis:as.factor(biotic_com_2)2)"]
se_woodland <- sqrt(var_woodland)

# z-scores for 85% CI
z <- qnorm(c(0.075, 0.925))

# Confidence intervals
ci_grassland <- beta_grassland + z * se_grassland
ci_woodland  <- beta_woodland  + z * se_woodland

# Tidy output
ci_table <- data.frame(
  Biotic_Community = c("Grassland", "Woodland"),
  Estimate = c(beta_grassland, beta_woodland),
  CI_85_lower = c(ci_grassland[1], ci_woodland[1]),
  CI_85_upper = c(ci_grassland[2], ci_woodland[2])
)

ci_table

# AICc Table with all possible turbine visibility models 

model_list <- list(
  "NDVI + slope" = Habitat,
  "Null" = Null,
  "Biotic community" = Bio_Com,
  "Biotic community + slope" = Bio_Com_Slope,
  "Slope" = Slope_Occu,
  "NDVI" = NDVI_Occu,
  "Canopy cover" = Canopy_Cov,
  "Canopy cover + slope" = Canopy_Cov_Slope,
  "Turbine visibility" = Turbine_Vis_Occu,
  "Turbine visibility + NDVI + slope" = Habitat_Add_Turbine_Vis,
  "Turbine visibility + NDVI" = NDVI_Add_Turbine_Vis,
  "Turbine visibility + slope" = Slope_Add_Turbine_Vis,
  "Turbine visibility x NDVI + slope" = NDVI_X_Turbine_Vis_Slope,
  "Turbine visibility x biotic community + slope" = Bio_Com_X_Turbine_Vis_Slope,
  "Turbine visibility x canopy cover + slope" = Canopy_X_Turbine_Vis_Slope
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine visibility" = "Null",
  "Turbine visibility + NDVI + slope" = "NDVI + slope",
  "Turbine visibility + NDVI" = "NDVI",
  "Turbine visibility + slope" = "Slope",
  "Turbine visibility x NDVI + slope" = "NDVI + slope",
  "Turbine visibility x biotic community + slope" = "Biotic community + slope",
  "Turbine visibility x canopy cover + slope" = "Canopy cover + slope"
)

# Ensure model_list matches wind_pairs
stopifnot(all(names(wind_pairs) %in% names(model_list)))
stopifnot(all(unlist(wind_pairs) %in% names(model_list)))

# The beta values (ß), standard errors, and 85% confidence intervals for 
  #parameter estimates within the best-supported model describing the effect of
  #the visibility of wind turbines from 50 cm off the ground on the probability 
  #of habitat selection (ψ) for gray foxes are listed in Table S3.26.

#### Turbine Distance Models ####

Turbine_Dist_Occu <- occu( ~ as.factor(cam_moved) ~ turbine_dist, occu.urci)

Habitat_Add_Turbine_Dist <- occu( ~ as.factor(cam_moved) 
                                  ~ turbine_dist + NDVI_1_5km + slope, 
                                  occu.urci)

NDVI_Add_Turbine_Dist <- occu( ~ as.factor(cam_moved) 
                               ~ turbine_dist + NDVI_1_5km, occu.urci)

Slope_Add_Turbine_Dist <- occu( ~ as.factor(cam_moved) 
                                ~ turbine_dist + slope, occu.urci)

Bio_Com_X_Turbine_Dist <- occu( ~ as.factor(cam_moved) 
                                ~ turbine_dist * as.factor(biotic_com_2) + 
                                  slope, occu.urci)

NDVI_X_Turbine_Dist <- occu( ~ as.factor(cam_moved) 
                             ~ turbine_dist * NDVI_1_5km + slope, occu.urci)

# 85% CI for bio community interaction 
# Extract coefficients and VCOV matrix
coef_est <- coef(Bio_Com_X_Turbine_Dist, type = "state")
vcov_mat <- vcov(Bio_Com_X_Turbine_Dist, type = "state")

beta_grassland <- coef_est["psi(turbine_dist)"]
se_grassland <- sqrt(vcov_mat["psi(turbine_dist)", 
                              "psi(turbine_dist)"])

beta_woodland <- coef_est["psi(turbine_dist)"] + 
  coef_est["psi(turbine_dist:as.factor(biotic_com_2)2)"]

var_woodland <- 
  vcov_mat["psi(turbine_dist)", "psi(turbine_dist)"] +
  vcov_mat["psi(turbine_dist:as.factor(biotic_com_2)2)",
           "psi(turbine_dist:as.factor(biotic_com_2)2)"] +
  2 * vcov_mat["psi(turbine_dist)", 
               "psi(turbine_dist:as.factor(biotic_com_2)2)"]
se_woodland <- sqrt(var_woodland)

# z-scores for 85% CI
z <- qnorm(c(0.075, 0.925))

# Confidence intervals
ci_grassland <- beta_grassland + z * se_grassland
ci_woodland  <- beta_woodland  + z * se_woodland

# Tidy output
ci_table <- data.frame(
  Biotic_Community = c("Grassland", "Woodland"),
  Estimate = c(beta_grassland, beta_woodland),
  CI_85_lower = c(ci_grassland[1], ci_woodland[1]),
  CI_85_upper = c(ci_grassland[2], ci_woodland[2])
)

ci_table 

# AICc Table with all possible turbine distance models 

model_list <- list(
  "NDVI + slope" = Habitat,
  "Null" = Null,
  "Biotic community + slope" = Bio_Com_Slope,
  "Slope" = Slope_Occu,
  "NDVI" = NDVI_Occu,
  "Turbine distance" = Turbine_Dist_Occu,
  "Turbine distance + NDVI + slope" = Habitat_Add_Turbine_Dist,
  "Turbine distance + NDVI" = NDVI_Add_Turbine_Dist,
  "Turbine distance + slope" = Slope_Add_Turbine_Dist,
  "Turbine distance x NDVI + slope" = NDVI_X_Turbine_Dist,
  "Turbine distance x biotic community + slope" = Bio_Com_X_Turbine_Dist
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine distance" = "Null",
  "Turbine distance + NDVI + slope" = "NDVI + slope",
  "Turbine distance + NDVI" = "NDVI",
  "Turbine distance + slope" = "Slope",
  "Turbine distance x NDVI + slope" = "NDVI + slope",
  "Turbine distance x biotic community + slope" = "Biotic community + slope" 
)

# Ensure model_list matches wind_pairs
stopifnot(all(names(wind_pairs) %in% names(model_list)))
stopifnot(all(unlist(wind_pairs) %in% names(model_list)))

# The beta values (ß), standard errors, and 85% confidence intervals for 
  #parameter estimates within the best-supported model describing the effect of
  #the distance to nearest wind turbine on the probability of habitat selection 
  #(ψ) for gray foxes are listed in Table S3.27.

#### Turbine Density Models ####

Turbine_Dense_Occu <- occu( ~ as.factor(cam_moved) 
                            ~ turbine_density_1_5km, occu.urci)

Habitat_Add_Turbine_Dense <- occu( ~ as.factor(cam_moved) 
                                   ~ turbine_density_1_5km + NDVI_1_5km + slope, 
                                     occu.urci)

NDVI_Add_Turbine_Dense <- occu( ~ as.factor(cam_moved) 
                                ~ turbine_density_1_5km + NDVI_1_5km, occu.urci)

Slope_Add_Turbine_Dense <- occu( ~ as.factor(cam_moved) 
                                 ~ turbine_density_1_5km + slope, occu.urci)

Bio_Com_X_Turbine_Dense <- occu( ~ as.factor(cam_moved) 
                                 ~ turbine_density_1_5km * 
                                   as.factor(biotic_com_2) + slope, occu.urci,
                                   starts = c(0, 1, 0, 0.5, 1, 0, -3))

NDVI_X_Turbine_Dense <- occu( ~ as.factor(cam_moved) 
                              ~ turbine_density_1_5km * NDVI_1_5km + slope, 
                                occu.urci)

# 85% CI for bio community interaction 

# Extract coefficients and VCOV matrix
coef_est <- coef(Bio_Com_X_Turbine_Dense, type = "state")
vcov_mat <- vcov(Bio_Com_X_Turbine_Dense, type = "state")

beta_grassland <- coef_est["psi(turbine_density_1_5km)"]
se_grassland <- sqrt(vcov_mat["psi(turbine_density_1_5km)", 
                              "psi(turbine_density_1_5km)"])

beta_woodland <- coef_est["psi(turbine_density_1_5km)"] + 
  coef_est["psi(turbine_density_1_5km:as.factor(biotic_com_2)2)"]

var_woodland <- 
  vcov_mat["psi(turbine_density_1_5km)", "psi(turbine_density_1_5km)"] +
  vcov_mat["psi(turbine_density_1_5km:as.factor(biotic_com_2)2)",
           "psi(turbine_density_1_5km:as.factor(biotic_com_2)2)"] +
  2 * vcov_mat["psi(turbine_density_1_5km)", 
               "psi(turbine_density_1_5km:as.factor(biotic_com_2)2)"]
se_woodland <- sqrt(var_woodland)

# z-scores for 85% CI
z <- qnorm(c(0.075, 0.925))

# Confidence intervals
ci_grassland <- beta_grassland + z * se_grassland
ci_woodland  <- beta_woodland  + z * se_woodland

# Tidy output
ci_table <- data.frame(
  Biotic_Community = c("Grassland", "Woodland"),
  Estimate = c(beta_grassland, beta_woodland),
  CI_85_lower = c(ci_grassland[1], ci_woodland[1]),
  CI_85_upper = c(ci_grassland[2], ci_woodland[2])
)

ci_table

# AICc Table with all possible turbine density models 

model_list <- list(
  "NDVI + slope" = Habitat,
  "Null" = Null,
  "Biotic community + slope" = Bio_Com_Slope,
  "Slope" = Slope_Occu,
  "NDVI" = NDVI_Occu,
  "Turbine density" = Turbine_Dense_Occu,
  "Turbine density + NDVI + slope" = Habitat_Add_Turbine_Dense,
  "Turbine density + NDVI" = NDVI_Add_Turbine_Dense,
  "Turbine density + slope" = Slope_Add_Turbine_Dense,
  "Turbine density x NDVI + slope" = NDVI_X_Turbine_Dense,
  "Turbine density x biotic community + slope" = Bio_Com_X_Turbine_Dense
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine density" = "Null",
  "Turbine density + NDVI + slope" = "NDVI + slope",
  "Turbine density + NDVI" = "NDVI",
  "Turbine density + slope" = "Slope",
  "Turbine density x NDVI + slope" = "NDVI + slope",
  "Turbine density x biotic community + slope" = 
    "Biotic community + slope"
)

# Ensure model_list matches wind_pairs
stopifnot(all(names(wind_pairs) %in% names(model_list)))
stopifnot(all(unlist(wind_pairs) %in% names(model_list)))

# The beta values (ß), standard errors, and 85% confidence intervals for 
  #parameter estimates within the best-supported model describing the effect of
  #the density of wind turbines within a 1.5 km squared radius surrounding the 
  #site on the probability of habitat selection (ψ) for gray foxes are listed 
  #in Table S3.28.

#### Access Road Distance Models ####

Turbine_Rd_Dist_Occu <- occu( ~ as.factor(cam_moved) 
                              ~ turbine_rd_dist, occu.urci)

Habitat_Add_Turbine_Rd_Dist <- occu( ~ as.factor(cam_moved) 
                                     ~ turbine_rd_dist + NDVI_1_5km + slope, 
                                     occu.urci)

NDVI_Add_Turbine_Rd_Dist <- occu( ~ as.factor(cam_moved) 
                                  ~ turbine_rd_dist + NDVI_1_5km, occu.urci)

Slope_Add_Turbine_Rd_Dist <- occu( ~ as.factor(cam_moved) 
                                   ~ turbine_rd_dist + slope, occu.urci)

Bio_Com_X_Turbine_Rd_Dist <- occu( ~ as.factor(cam_moved) 
                                   ~ turbine_rd_dist * as.factor(biotic_com_2) +
                                     slope, occu.urci)

NDVI_X_Turbine_Rd_Dist <- occu( ~ as.factor(cam_moved) 
                                ~ turbine_rd_dist * NDVI_1_5km + slope, 
                                occu.urci)

# 85% CI for bio community interaction 
# Extract coefficients and VCOV matrix
coef_est <- coef(Bio_Com_X_Turbine_Rd_Dist, type = "state")
vcov_mat <- vcov(Bio_Com_X_Turbine_Rd_Dist, type = "state")

beta_grassland <- coef_est["psi(turbine_rd_dist)"]
se_grassland <- sqrt(vcov_mat["psi(turbine_rd_dist)", 
                              "psi(turbine_rd_dist)"])

beta_woodland <- coef_est["psi(turbine_rd_dist)"] + 
  coef_est["psi(turbine_rd_dist:as.factor(biotic_com_2)2)"]

var_woodland <- 
  vcov_mat["psi(turbine_rd_dist)", "psi(turbine_rd_dist)"] +
  vcov_mat["psi(turbine_rd_dist:as.factor(biotic_com_2)2)",
           "psi(turbine_rd_dist:as.factor(biotic_com_2)2)"] +
  2 * vcov_mat["psi(turbine_rd_dist)", 
               "psi(turbine_rd_dist:as.factor(biotic_com_2)2)"]
se_woodland <- sqrt(var_woodland)

# z-scores for 85% CI
z <- qnorm(c(0.075, 0.925))

# Confidence intervals
ci_grassland <- beta_grassland + z * se_grassland
ci_woodland  <- beta_woodland  + z * se_woodland

# Tidy output
ci_table <- data.frame(
  Biotic_Community = c("Grassland", "Woodland"),
  Estimate = c(beta_grassland, beta_woodland),
  CI_85_lower = c(ci_grassland[1], ci_woodland[1]),
  CI_85_upper = c(ci_grassland[2], ci_woodland[2])
)

ci_table

# AICc Table with all possible access rd distance models 

model_list <- list(
  "NDVI + slope" = Habitat,
  "Null" = Null,
  "Biotic community + slope" = Bio_Com_Slope,
  "Slope" = Slope_Occu,
  "NDVI" = NDVI_Occu,
  "Turbine road distance" = Turbine_Rd_Dist_Occu,
  "Turbine road distance + NDVI + slope" = Habitat_Add_Turbine_Rd_Dist,
  "Turbine road distance + NDVI" = NDVI_Add_Turbine_Rd_Dist,
  "Turbine road distance + slope" = Slope_Add_Turbine_Rd_Dist,
  "Turbine road distance x NDVI + slope" = NDVI_X_Turbine_Rd_Dist,
  "Turbine road distance x biotic community + slope" = Bio_Com_X_Turbine_Rd_Dist
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine road distance" = "Null",
  "Turbine road distance + NDVI + slope" = "NDVI + slope",
  "Turbine road distance + NDVI" = "NDVI",
  "Turbine road distance + slope" = "Slope",
  "Turbine road distance x NDVI + slope" = "NDVI + slope",
  "Turbine road distance x biotic community + slope" = "Biotic community + slope" 
)

# Ensure model_list matches wind_pairs
stopifnot(all(names(wind_pairs) %in% names(model_list)))
stopifnot(all(unlist(wind_pairs) %in% names(model_list)))

#### Access Road Density Models ####

Turbine_Rd_Dense_Occu <- occu( ~ as.factor(cam_moved) 
                               ~ turbine_rd_density_1_5km, occu.urci)

Habitat_Add_Turbine_Rd_Dense <- occu( ~ as.factor(cam_moved) 
                                      ~ turbine_rd_density_1_5km + NDVI_1_5km +
                                        slope, occu.urci)

NDVI_Add_Turbine_Rd_Dense <- occu( ~ as.factor(cam_moved) 
                                   ~ turbine_rd_density_1_5km + NDVI_1_5km, 
                                     occu.urci)

Slope_Add_Turbine_Rd_Dense <- occu( ~ as.factor(cam_moved) 
                                    ~ turbine_rd_density_1_5km + slope, 
                                      occu.urci)

Bio_Com_X_Turbine_Rd_Dense <- occu( ~ as.factor(cam_moved) 
                                    ~ turbine_rd_density_1_5km *
                                      as.factor(biotic_com_2) + slope, 
                                      occu.urci)

NDVI_X_Turbine_Rd_Dense <- occu( ~ as.factor(cam_moved) 
                                 ~ turbine_rd_density_1_5km * NDVI_1_5km + 
                                   slope, occu.urci)

# 85% CI for bio community interaction 
# Extract coefficients and VCOV matrix
coef_est <- coef(Bio_Com_X_Turbine_Rd_Dense, type = "state")
vcov_mat <- vcov(Bio_Com_X_Turbine_Rd_Dense, type = "state")

beta_grassland <- coef_est["psi(turbine_rd_density_1_5km)"]
se_grassland <- sqrt(vcov_mat["psi(turbine_rd_density_1_5km)", 
                              "psi(turbine_rd_density_1_5km)"])

beta_woodland <- coef_est["psi(turbine_rd_density_1_5km)"] + 
  coef_est["psi(turbine_rd_density_1_5km:as.factor(biotic_com_2)2)"]

var_woodland <- 
  vcov_mat["psi(turbine_rd_density_1_5km)", "psi(turbine_rd_density_1_5km)"] +
  vcov_mat["psi(turbine_rd_density_1_5km:as.factor(biotic_com_2)2)",
           "psi(turbine_rd_density_1_5km:as.factor(biotic_com_2)2)"] +
  2 * vcov_mat["psi(turbine_rd_density_1_5km)", 
               "psi(turbine_rd_density_1_5km:as.factor(biotic_com_2)2)"]
se_woodland <- sqrt(var_woodland)

# z-scores for 85% CI
z <- qnorm(c(0.075, 0.925))

# Confidence intervals
ci_grassland <- beta_grassland + z * se_grassland
ci_woodland  <- beta_woodland  + z * se_woodland

# Tidy output
ci_table <- data.frame(
  Biotic_Community = c("Grassland", "Woodland"),
  Estimate = c(beta_grassland, beta_woodland),
  CI_85_lower = c(ci_grassland[1], ci_woodland[1]),
  CI_85_upper = c(ci_grassland[2], ci_woodland[2])
)

ci_table 

# AICc Table with all possible access rd density models 

model_list <- list(
  "NDVI + slope" = Habitat,
  "Null" = Null,
  "Biotic community + slope" = Bio_Com_Slope,
  "Slope" = Slope_Occu,
  "NDVI" = NDVI_Occu,
  "Turbine road density" = Turbine_Rd_Dense_Occu,
  "Turbine road density + NDVI + slope" = Habitat_Add_Turbine_Rd_Dense,
  "Turbine road density + NDVI" = NDVI_Add_Turbine_Rd_Dense,
  "Turbine road density + slope" = Slope_Add_Turbine_Rd_Dense,
  "Turbine road density x NDVI + slope" = NDVI_X_Turbine_Rd_Dense,
  "Turbine road density x biotic community + slope" = Bio_Com_X_Turbine_Rd_Dense
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine road density" = "Null",
  "Turbine road density + NDVI + slope" = "NDVI + slope",
  "Turbine road density + NDVI" = "NDVI",
  "Turbine road density + slope" = "Slope",
  "Turbine road density x NDVI + slope" = "NDVI + slope",
  "Turbine road density x biotic community + slope" = "Biotic community + slope"
)

# Ensure model_list matches wind_pairs
stopifnot(all(names(wind_pairs) %in% names(model_list)))
stopifnot(all(unlist(wind_pairs) %in% names(model_list)))

#### Additive Wind Variables ####

# 2 wind variable mods 
Add1 <- occu( ~ as.factor(cam_moved) ~ X50cm_turbine_vis + 
                turbine_density_1_5km + slope, occu.urci)

Add2 <- occu( ~ as.factor(cam_moved) ~ X50cm_turbine_vis + 
                turbine_dist + slope, occu.urci)

Add3 <- occu( ~ as.factor(cam_moved) ~ X50cm_turbine_vis + 
                as.factor(turbine_interior) + slope, occu.urci)

Add4 <- occu( ~ as.factor(cam_moved) ~ X50cm_turbine_vis *
                as.factor(biotic_com_2) + 
                turbine_density_1_5km + slope, occu.urci)

Add5 <- occu( ~ as.factor(cam_moved) ~ X50cm_turbine_vis *
                as.factor(biotic_com_2) + 
                as.factor(turbine_interior) + slope, occu.urci)

Add6 <- occu( ~ as.factor(cam_moved) ~ X50cm_turbine_vis *
                as.factor(biotic_com_2) + 
                turbine_dist + slope, occu.urci)

Add7 <- occu( ~ as.factor(cam_moved) ~ X50cm_turbine_vis + 
                turbine_density_1_5km + NDVI_1_5km, occu.urci)

Add8 <- occu( ~ as.factor(cam_moved) ~ X50cm_turbine_vis + 
                turbine_dist + NDVI_1_5km, occu.urci)

Add9 <- occu( ~ as.factor(cam_moved) ~ X50cm_turbine_vis + 
                as.factor(turbine_interior) + NDVI_1_5km, occu.urci)

Add10 <- occu( ~ as.factor(cam_moved) ~ X50cm_turbine_vis + 
                turbine_dist + NDVI_1_5km + slope, occu.urci)

Add11 <- occu( ~ as.factor(cam_moved) ~ turbine_dist + 
                as.factor(turbine_interior) + NDVI_1_5km + slope, occu.urci)

Add12 <- occu( ~ as.factor(cam_moved) ~ X50cm_turbine_vis + 
                 as.factor(turbine_interior) + NDVI_1_5km + slope, occu.urci)

Add13 <- occu( ~ as.factor(cam_moved) ~ turbine_density_1_5km + 
                 X50cm_turbine_vis + NDVI_1_5km + slope, occu.urci)

Add14 <- occu( ~ as.factor(cam_moved) ~ turbine_dist + 
                 as.factor(turbine_interior) * NDVI_1_5km + slope, occu.urci)

Add15 <- occu( ~ as.factor(cam_moved) ~ X50cm_turbine_vis + 
                 as.factor(turbine_interior) * NDVI_1_5km + slope, occu.urci)

Add16 <- occu( ~ as.factor(cam_moved) ~ turbine_dist + 
                 as.factor(turbine_interior) * NDVI_1_5km, occu.urci)

Add17 <- occu( ~ as.factor(cam_moved) ~ X50cm_turbine_vis + 
                 as.factor(turbine_interior) * NDVI_1_5km, occu.urci)

#3 wind variable mods

Add18 <- occu( ~ as.factor(cam_moved) ~ turbine_dist + X50cm_turbine_vis + 
                 as.factor(turbine_interior) * NDVI_1_5km + slope, occu.urci)

Add19 <- occu( ~ as.factor(cam_moved) ~ turbine_dist + X50cm_turbine_vis + 
                 as.factor(turbine_interior) * NDVI_1_5km, occu.urci)

Add20 <- occu( ~ as.factor(cam_moved) ~ turbine_dist + X50cm_turbine_vis + 
                 as.factor(turbine_interior) + NDVI_1_5km, occu.urci)

Add21 <- occu( ~ as.factor(cam_moved) ~ turbine_dist + X50cm_turbine_vis + 
                 as.factor(turbine_interior) + NDVI_1_5km + slope, occu.urci)

# AICc Table with all possible additive models 

model_list <- list(
  "NDVI + slope" = Habitat,
  "Null" = Null,
  "Slope" = Slope_Occu,
  "NDVI" = NDVI_Occu,
  "Biotic community + slope" = Bio_Com_Slope,
  "Turbine visibility + turbine density + slope" = Add1,
  "Turbine visibility + turbine distance + slope" = Add2,
  "Turbine visibility + turbine interior + slope" = Add3,
  "Turbine visibility * biotic community + turbine density + slope" = Add4,
  "Turbine visibility * biotic community + turbine interior + slope" = Add5,
  "Turbine visibility * biotic community + turbine distance + slope" = Add6,
  "Turbine visibility + turbine density + NDVI" = Add7,
  "Turbine visibility + turbine distance + NDVI" = Add8,
  "Turbine visibility + turbine interior + NDVI" = Add9,
  "Turbine visibility + turbine distance + NDVI + slope" = Add10,
  "Turbine interior + turbine distance + NDVI + slope" = Add11,
  "Turbine visibility + turbine interior + NDVI + slope" = Add12,
  "Turbine visibility + turbine density + NDVI + slope" = Add13,
  "Turbine distance + turbine interior * NDVI + slope" = Add14,
  "Turbine visibility + turbine interior * NDVI + slope" = Add15,
  "Turbine distance + turbine visibility + turbine interior * NDVI + slope" = 
    Add18,
  "Turbine distance + turbine visibility + turbine interior + NDVI" = Add20,
  "Turbine distance + turbine visibility + turbine interior + NDVI + slope" = 
    Add21
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine visibility + turbine density + slope" = "Slope",
  "Turbine visibility + turbine distance + slope" = "Slope",
  "Turbine visibility + turbine interior + slope" = "Slope",
  "Turbine visibility * biotic community + turbine density + slope" = 
    "Biotic community + slope",
  "Turbine visibility * biotic community + turbine interior + slope" = 
    "Biotic community + slope",
  "Turbine visibility * biotic community + turbine distance + slope" = 
    "Biotic community + slope",
  "Turbine visibility + turbine density + NDVI" = "NDVI",
  "Turbine visibility + turbine distance + NDVI" = "NDVI",
  "Turbine visibility + turbine interior + NDVI" = "NDVI",
  "Turbine visibility + turbine distance + NDVI + slope" = "NDVI + slope",
  "Turbine interior + turbine distance + NDVI + slope" = "NDVI + slope",
  "Turbine visibility + turbine interior + NDVI + slope" = "NDVI + slope",
  "Turbine visibility + turbine density + NDVI + slope" = "NDVI + slope",
  "Turbine distance + turbine interior * NDVI + slope" = "NDVI + slope",
  "Turbine visibility + turbine interior * NDVI + slope" = "NDVI + slope",
  "Turbine distance + turbine visibility + turbine interior * NDVI + slope" = 
    "NDVI + slope",
  "Turbine distance + turbine visibility + turbine interior + NDVI + slope" = 
    "NDVI + slope"
)

# Ensure model_list matches wind_pairs
stopifnot(all(names(wind_pairs) %in% names(model_list)))
stopifnot(all(unlist(wind_pairs) %in% names(model_list)))

# The beta values (ß), standard errors, and 85% confidence intervals for 
  #parameter estimates within the best-supported model describing the effect of
  #the density of wind turbines within a 1.5 km squared radius surrounding the 
  #site and the visibility of turbines from 50 cm off the ground on the 
  #probability of habitat selection (ψ) for gray foxes are listed in Table S3.29.


#####################################################################
# CREATE MODEL SELECTION TABLES #
# RUN AFTER THE AICc TABLE IS CREATED FOR EACH GROUP OF WIND MODELS #
#####################################################################


# 1. Function to compute 85% CI for occupancy (psi) or detection
get_85CI <- function(model, type = "state") {
  est <- coef(model, type = type)
  if(length(est) == 0) return(NULL)
  
  vc <- vcov(model, type = type)
  se <- sqrt(diag(vc))
  z <- qnorm(0.925)  # 85% CI
  
  df <- data.frame(
    Parameter = names(est),
    CI_string = sprintf("%.2f [%.2f, %.2f]", est, est - z*se, est + z*se)
  )
  return(df)
}

# 2. Model selection table
top_mods <- model.sel(model_list)
top_mods_df <- as.data.frame(top_mods)
top_mods_df$Model <- rownames(top_mods_df)

# 3. Compute occupancy and detection CIs for all models
ci_psi_list <- lapply(model_list, get_85CI, type = "state")
ci_det_list <- lapply(model_list, get_85CI, type = "det")

# Helper to pivot wide and add model name
pivot_ci_wide <- function(ci_list) {
  df <- bind_rows(
    lapply(names(ci_list), function(name) {
      x <- ci_list[[name]]
      if(!is.null(x)) x$Model <- name
      x
    })
  )
  if(nrow(df) == 0) return(data.frame(Model = character(0)))
  
  df %>% pivot_wider(names_from = Parameter, values_from = CI_string)
}

ci_psi_wide <- pivot_ci_wide(ci_psi_list)
ci_det_wide <- pivot_ci_wide(ci_det_list)

# 4. Log-likelihood ratio tests for wind vs non-wind models
llr_results <- lapply(names(wind_pairs), function(mod_with_wind) {
  mod_without <- wind_pairs[[mod_with_wind]]
  m1 <- model_list[[mod_without]]
  m2 <- model_list[[mod_with_wind]]
  
  lrt <- 2 * (logLik(m2) - logLik(m1))
  df_diff <- attr(logLik(m2), "df") - attr(logLik(m1), "df")
  pval <- pchisq(lrt, df = df_diff, lower.tail = FALSE)
  
  data.frame(
    Model = mod_with_wind,
    LRT_stat = as.numeric(lrt),
    LRT_p = as.numeric(pval)
  )
})
llr_df <- bind_rows(llr_results)

# 5. Pairwise evidence ratios
w <- top_mods_df$weight
names(w) <- top_mods_df$Model
pairwise_ER <- lapply(names(wind_pairs), function(mod_with_wind) {
  mod_without <- wind_pairs[[mod_with_wind]]
  ER <- if(all(c(mod_with_wind, mod_without) %in% names(w))) 
    w[mod_with_wind] / w[mod_without] else NA
  data.frame(Model = mod_with_wind, Evidence_Ratio_vs_NonWind = ER)
})
pairwise_ER_df <- bind_rows(pairwise_ER)

# 6. Merge everything
final_table <- top_mods_df %>%
  left_join(ci_psi_wide, by = "Model") %>%
  left_join(ci_det_wide, by = "Model", suffix = c("_psi", "_det")) %>%
  left_join(llr_df, by = "Model") %>%
  left_join(pairwise_ER_df, by = "Model")

# 7. Export tables to .csv files

#write.csv(final_table, "gray_fox_turbine_interior.csv", row.names = FALSE)

#write.csv(final_table, "gray_fox_turbine_vis.csv", row.names = FALSE)

#write.csv(final_table, "gray_fox_turbine_dist.csv", row.names = FALSE)

#write.csv(final_table, "gray_fox_turbine_density.csv", row.names = FALSE)

#write.csv(final_table, "gray_fox_turbine_rd_dist.csv", row.names = FALSE)

#write.csv(final_table, "gray_fox_turbine_rd_density.csv", row.names = FALSE)

#write.csv(final_table, "gray_fox_additive_mods.csv", row.names = FALSE)

########################################################
######################### END ##########################
########################################################