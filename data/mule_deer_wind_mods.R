###########################################################
########### Code for: wind models for mule deer ###########
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
#is the homewd directory on line 23

#Set home working directory
#homewd = "C:/Users/ionar/Desktop/R Repository/Wind Energy/"
homewd = "<insert your folder here and end with a forward slash>"

#Set wd to data folder on your local computer 
setwd(paste0(homewd, "data/"))

###########################################################
# SETUP CODE FOR MULE DEER OCCUPANCY MODELS #
###########################################################

# Read in edited .csv
detHist <- read.csv(file = "mule deer detection hist.csv", row.names = 1)

# Change from integer to numeric
detHist <- detHist %>%
  mutate(across(1:51, as.numeric))

sum(as.matrix(detHist), na.rm = TRUE)

# Create detection history 
occu.odhe <- unmarkedFrameOccu(y=detHist, 
                               siteCovs = site.covs.scaled, 
                               obsCovs = obsCovs.scaled)

#################################################
# CREATE THE WIND MODELS #
#################################################

#### Models with no change over all turbine variables ####

Habitat <- occu( ~ max_trig_dist + NDVI_1_9km 
                 ~ shrub_yucca_density + veg_cover, occu.odhe)

Vis_Obs_Occu <- occu( ~ max_trig_dist + NDVI_1_9km 
                      ~ veg_cover, occu.odhe)

Canopy_Shrub <- occu( ~ max_trig_dist + NDVI_1_9km 
                      ~ shrub_yucca_density + canopy_cov, occu.odhe)

Shrub_Dense_Occu <- occu( ~ max_trig_dist + NDVI_1_9km 
                          ~ shrub_yucca_density, occu.odhe)

Bio_Com_Occu <- occu( ~ max_trig_dist + NDVI_1_9km 
                      ~ as.factor(biotic_com_2), occu.odhe)

Bio_Com_Shrub <- occu( ~ max_trig_dist + NDVI_1_9km 
                       ~ as.factor(biotic_com_2) + shrub_yucca_density, 
                         occu.odhe)

Null <- occu( ~ max_trig_dist + NDVI_1_9km ~ 1, occu.odhe)

## The table with the comparison of the relative weight of evidence between
  #occupancy models separately ranked for the effect of each wind energy 
  #variable on the probability of habitat selection (ψ) of mule deer is found 
  #in Table S3.13.

#### Turbine Interior Models ####

Bio_Com_X_Turbine_Int_Add_Shrub <- occu( ~ max_trig_dist + NDVI_1_9km 
                                         ~ as.factor(turbine_interior) *
                                           as.factor(biotic_com_2) +
                                           shrub_yucca_density, occu.odhe)

Vis_Obs_X_Turbine_Int_Add_Shrub <- occu( ~ max_trig_dist + NDVI_1_9km 
                                         ~ as.factor(turbine_interior) *
                                           veg_cover +
                                           shrub_yucca_density, occu.odhe)

Habitat_Add_Turbine_Int <- occu( ~ max_trig_dist + NDVI_1_9km  
                                 ~ shrub_yucca_density + veg_cover +
                                   as.factor(turbine_interior), occu.odhe)

Turbine_Int_Occu <- occu( ~ max_trig_dist + NDVI_1_9km 
                          ~ as.factor(turbine_interior), occu.odhe)

Vis_Obs_Add_Turbine_Int <- occu( ~ max_trig_dist + NDVI_1_9km 
                                 ~ as.factor(turbine_interior) +
                                   veg_cover, occu.odhe)

Shrub_Dense_Add_Turbine_Int <- occu( ~ max_trig_dist + NDVI_1_9km 
                                     ~ as.factor(turbine_interior) +
                                       shrub_yucca_density, occu.odhe)   

# 85% CI for interior x vis obs
# Extract coefficients and VCOV matrix
coef_est <- coef(Vis_Obs_X_Turbine_Int_Add_Shrub, type = "state")
vcov_mat <- vcov(Vis_Obs_X_Turbine_Int_Add_Shrub, type = "state")

# effect OUTSIDE wind farm (reference = 0) 
beta_outside <- coef_est["psi(veg_cover)"]

se_outside <- sqrt(
  vcov_mat["psi(veg_cover)", "psi(veg_cover)"]
)

# effect INSIDE wind farm (reference + interaction)
beta_inside <- coef_est["psi(veg_cover)"] +
  coef_est["psi(as.factor(turbine_interior)1:veg_cover)"]

var_inside <-
  vcov_mat["psi(veg_cover)", "psi(veg_cover)"] +
  vcov_mat["psi(as.factor(turbine_interior)1:veg_cover)",
           "psi(as.factor(turbine_interior)1:veg_cover)"] +
  2 * vcov_mat["psi(veg_cover)",
               "psi(as.factor(turbine_interior)1:veg_cover)"]

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

# 85% CI for interior x bio comm
# Extract coefficients and VCOV matrix
coef_est <- coef(Bio_Com_X_Turbine_Int_Add_Shrub, type = "state")
vcov_mat <- vcov(Bio_Com_X_Turbine_Int_Add_Shrub, type = "state")

# Exact coefficient names from your model
intercept_name   <- "psi(Int)"
turbine_name     <- "psi(as.factor(turbine_interior)1)"
habitat_name     <- "psi(as.factor(biotic_com_2)2)"
interaction_name <- "psi(as.factor(turbine_interior)1:as.factor(biotic_com_2)2)"
shrub_yucca_density_name       <- "psi(shrub_yucca_density)"

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

# AICc Table with all possible turbine interior models 

model_list <- list(
  "Visual obstruction + shrub density" = Habitat,
  "Turbine interior + visual obstruction + shrub density" =
    Habitat_Add_Turbine_Int,
  "Turbine interior x visual obstruction + shrub density" =
    Vis_Obs_X_Turbine_Int_Add_Shrub,
  "Turbine interior x biotic community + shrub density" = 
    Bio_Com_X_Turbine_Int_Add_Shrub,
  "Null" = Null,
  "Turbine interior" = Turbine_Int_Occu,
  "Visual obstruction" = Vis_Obs_Occu,
  "Shrub density" = Shrub_Dense_Occu,
  "Biotic community + shrub density" = Bio_Com_Shrub,
  "Turbine interior + visual obstruction" = Vis_Obs_Add_Turbine_Int,
  "Turbine interior + shrub density" = Shrub_Dense_Add_Turbine_Int
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine interior" = "Null",
  "Turbine interior + visual obstruction + shrub density" = 
    "Visual obstruction + shrub density",
  "Turbine interior + shrub density" = "Shrub density",
  "Turbine interior + visual obstruction" = "Visual obstruction",
  "Turbine interior x visual obstruction + shrub density" = 
    "Visual obstruction + shrub density",
  "Turbine interior x biotic community + shrub density" = 
    "Biotic community + shrub density"
)

# Ensure model_list matches wind_pairs
stopifnot(all(names(wind_pairs) %in% names(model_list)))
stopifnot(all(unlist(wind_pairs) %in% names(model_list)))

#### Turbine Visibility Models ####

Bio_Com_X_Turbine_Vis <- occu( ~ max_trig_dist + NDVI_1_9km 
                               ~ X150cm_turbine_vis * as.factor(biotic_com_2) +
                                 shrub_yucca_density, occu.odhe,
                               starts = c(-0.05, 1, -1,1,1,-3,-0.5,0.5))

Vis_Obs_X_Turbine_Vis <- occu( ~ max_trig_dist + NDVI_1_9km 
                               ~ veg_cover * X150cm_turbine_vis +
                                 shrub_yucca_density, occu.odhe)

Canopy_Cov_X_Turbine_Vis <- occu( ~ max_trig_dist + NDVI_1_9km 
                                  ~ canopy_cov * X150cm_turbine_vis +
                                    shrub_yucca_density, occu.odhe)

Habitat_Add_Turbine_Vis <- occu( ~ max_trig_dist + NDVI_1_9km 
                                 ~ shrub_yucca_density + veg_cover +
                                   X150cm_turbine_vis, occu.odhe)

Turbine_Vis_Occu <- occu( ~ max_trig_dist + NDVI_1_9km 
                          ~ X150cm_turbine_vis, occu.odhe)

Vis_Obs_Add_Turbine_Vis <- occu( ~ max_trig_dist + NDVI_1_9km 
                                 ~ X150cm_turbine_vis + veg_cover,
                                 occu.odhe)

Shrub_Dense_Add_Turbine_Vis <- occu( ~ max_trig_dist + NDVI_1_9km  
                                     ~ X150cm_turbine_vis + shrub_yucca_density, 
                                     occu.odhe)    

# 85% CI for bio community interaction 
# Extract coefficients and VCOV matrix
coef_est <- coef(Bio_Com_X_Turbine_Vis, type = "state")
vcov_mat <- vcov(Bio_Com_X_Turbine_Vis, type = "state")

beta_grassland <- coef_est["psi(X150cm_turbine_vis)"]
se_grassland <- sqrt(vcov_mat["psi(X150cm_turbine_vis)", 
                              "psi(X150cm_turbine_vis)"])

beta_woodland <- coef_est["psi(X150cm_turbine_vis)"] + 
  coef_est["psi(X150cm_turbine_vis:as.factor(biotic_com_2)2)"]

var_woodland <- 
  vcov_mat["psi(X150cm_turbine_vis)", "psi(X150cm_turbine_vis)"] +
  vcov_mat["psi(X150cm_turbine_vis:as.factor(biotic_com_2)2)",
           "psi(X150cm_turbine_vis:as.factor(biotic_com_2)2)"] +
  2 * vcov_mat["psi(X150cm_turbine_vis)", 
               "psi(X150cm_turbine_vis:as.factor(biotic_com_2)2)"]
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
  "Visual obstruction + shrub density" = Habitat,
  "Turbine visibility + visual obstruction + shrub density" = 
    Habitat_Add_Turbine_Vis,
  "Turbine visibility x visual obstruction + shrub density" = 
    Vis_Obs_X_Turbine_Vis,
  "Turbine visibility x biotic community + shrub density" = 
    Bio_Com_X_Turbine_Vis,
  "Null" = Null,
  "Turbine visibility" = Turbine_Vis_Occu,
  "Visual obstruction" = Vis_Obs_Occu,
  "Shrub density" = Shrub_Dense_Occu,
  "Biotic community + shrub density" = Bio_Com_Shrub,
  "Canopy cover + shrub density" = Canopy_Shrub,
  "Turbine visibility + shrub density" = Shrub_Dense_Add_Turbine_Vis,
  "Turbine visibility + visual obstruction" = Vis_Obs_Add_Turbine_Vis,
  "Turbine visibility x canopy cover + shrub density" = Canopy_Cov_X_Turbine_Vis
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine visibility" = "Null",
  "Turbine visibility + visual obstruction + shrub density" = 
    "Visual obstruction + shrub density",
  "Turbine visibility + shrub density" = "Shrub density",
  "Turbine visibility + visual obstruction" = "Visual obstruction",
  "Turbine visibility x visual obstruction + shrub density" = 
    "Visual obstruction + shrub density",
  "Turbine visibility x biotic community + shrub density" = 
    "Biotic community + shrub density",
  "Turbine visibility x canopy cover + shrub density" = 
    "Canopy cover + shrub density"
)

# Ensure model_list matches wind_pairs
stopifnot(all(names(wind_pairs) %in% names(model_list)))
stopifnot(all(unlist(wind_pairs) %in% names(model_list)))

# The beta values (ß), standard errors, and 85% confidence intervals for 
  #parameter estimates within the best-supported model describing the effect of
  #the visibility of wind turbines from 150 cm off the ground on the probability 
  #of habitat selection (ψ) for mule deer are listed in Table S3.22.

#### Turbine Distance Models #####

Bio_Com_X_Turbine_Dist <- occu( ~ max_trig_dist + NDVI_1_9km  
                                ~ turbine_dist * as.factor(biotic_com_2) +
                                  shrub_yucca_density, occu.odhe)

Vis_Obs_X_Turbine_Dist <- occu( ~ max_trig_dist + NDVI_1_9km  
                                ~ veg_cover * turbine_dist +
                                  shrub_yucca_density, occu.odhe)

Habitat_Add_Turbine_Dist <- occu( ~ max_trig_dist + NDVI_1_9km 
                                  ~ shrub_yucca_density + veg_cover +
                                    turbine_dist, occu.odhe)

Turbine_Dist_Occu <- occu( ~ max_trig_dist + NDVI_1_9km 
                           ~ turbine_dist, occu.odhe)

Vis_Obs_Add_Turbine_Dist <- occu( ~ max_trig_dist + NDVI_1_9km 
                                  ~ turbine_dist + veg_cover, 
                                  occu.odhe)

Shrub_Dense_Add_Turbine_Dist <- occu( ~ max_trig_dist + NDVI_1_9km  
                                      ~ turbine_dist + shrub_yucca_density, 
                                      occu.odhe)

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
  "Visual obstruction + shrub density" = Habitat,
  "Turbine distance + visual obstruction + shrub density" = Habitat_Add_Turbine_Dist,
  "Turbine distance x visual obstruction + shrub density" = Vis_Obs_X_Turbine_Dist,
  "Turbine distance x biotic community + shrub density" = Bio_Com_X_Turbine_Dist,
  "Null" = Null,
  "Turbine distance" = Turbine_Dist_Occu,
  "Visual obstruction" = Vis_Obs_Occu,
  "Shrub density" = Shrub_Dense_Occu,
  "Biotic community + shrub density" = Bio_Com_Shrub,
  "Turbine distance + visual obstruction" = Vis_Obs_Add_Turbine_Dist,
  "Turbine distance + shrub density" = Shrub_Dense_Add_Turbine_Dist
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine distance" = "Null",
  "Turbine distance + visual obstruction + shrub density" = 
    "Visual obstruction + shrub density",
  "Turbine distance + shrub density" = "Shrub density",
  "Turbine distance + visual obstruction" = "Visual obstruction",
  "Turbine distance x visual obstruction + shrub density" = 
    "Visual obstruction + shrub density",
  "Turbine distance x biotic community + shrub density" = 
    "Biotic community + shrub density"
)

# Ensure model_list matches wind_pairs
stopifnot(all(names(wind_pairs) %in% names(model_list)))
stopifnot(all(unlist(wind_pairs) %in% names(model_list)))

# The beta values (ß), standard errors, and 85% confidence intervals for 
  #parameter estimates within the best-supported model describing the effect of
  #the distance to nearest wind turbine on the probability of habitat selection 
  #(ψ) for mule deer are listed in Table S3.23.

#### Turbine Density Models ####

Bio_Com_X_Turbine_Dense <- occu( ~ max_trig_dist + NDVI_1_9km 
                                 ~ turbine_density_1_9km * 
                                   as.factor(biotic_com_2) +
                                   shrub_yucca_density, occu.odhe)

Vis_Obs_X_Turbine_Dense <- occu( ~ max_trig_dist + NDVI_1_9km 
                                 ~ veg_cover * turbine_density_1_9km +
                                   shrub_yucca_density, occu.odhe)

Habitat_Add_Turbine_Dense <- occu( ~ max_trig_dist + NDVI_1_9km 
                                   ~ shrub_yucca_density + veg_cover +
                                     turbine_density_1_9km, occu.odhe)

Turbine_Dense_Occu <- occu( ~ max_trig_dist + NDVI_1_9km 
                            ~ turbine_density_1_9km, occu.odhe)

Vis_Obs_Add_Turbine_Dense <- occu( ~ max_trig_dist + NDVI_1_9km 
                                   ~ turbine_density_1_9km + veg_cover,
                                     occu.odhe)

Shrub_Dense_Add_Turbine_Dense <- occu( ~ max_trig_dist + NDVI_1_9km 
                                       ~ turbine_density_1_9km +
                                         shrub_yucca_density, occu.odhe)

# 85% CI for bio community interaction 
# Extract coefficients and VCOV matrix
coef_est <- coef(Bio_Com_X_Turbine_Dense, type = "state")
vcov_mat <- vcov(Bio_Com_X_Turbine_Dense, type = "state")

beta_grassland <- coef_est["psi(turbine_density_1_9km)"]
se_grassland <- sqrt(vcov_mat["psi(turbine_density_1_9km)", 
                              "psi(turbine_density_1_9km)"])

beta_woodland <- coef_est["psi(turbine_density_1_9km)"] + 
  coef_est["psi(turbine_density_1_9km:as.factor(biotic_com_2)2)"]

var_woodland <- 
  vcov_mat["psi(turbine_density_1_9km)", "psi(turbine_density_1_9km)"] +
  vcov_mat["psi(turbine_density_1_9km:as.factor(biotic_com_2)2)",
           "psi(turbine_density_1_9km:as.factor(biotic_com_2)2)"] +
  2 * vcov_mat["psi(turbine_density_1_9km)", 
               "psi(turbine_density_1_9km:as.factor(biotic_com_2)2)"]
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
  "Visual obstruction + shrub density" = Habitat,
  "Turbine density + visual obstruction + shrub density" = Habitat_Add_Turbine_Dense,
  "Turbine density x visual obstruction + shrub density" = Vis_Obs_X_Turbine_Dense,
  "Turbine density x biotic community + shrub density" = Bio_Com_X_Turbine_Dense,
  "Null" = Null,
  "Turbine density" = Turbine_Dense_Occu,
  "Visual obstruction" = Vis_Obs_Occu,
  "Shrub density" = Shrub_Dense_Occu,
  "Biotic community + shrub density" = Bio_Com_Shrub,
  "Turbine density + visual obstruction" = Vis_Obs_Add_Turbine_Dense,
  "Turbine density + shrub density" = Shrub_Dense_Add_Turbine_Dense
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine density" = "Null",
  "Turbine density + visual obstruction + shrub density" = 
    "Visual obstruction + shrub density",
  "Turbine density + shrub density" = "Shrub density",
  "Turbine density + visual obstruction" = "Visual obstruction",
  "Turbine density x visual obstruction + shrub density" = 
    "Visual obstruction + shrub density",
  "Turbine density x biotic community + shrub density" = 
    "Biotic community + shrub density"
)

# Ensure model_list matches wind_pairs
stopifnot(all(names(wind_pairs) %in% names(model_list)))
stopifnot(all(unlist(wind_pairs) %in% names(model_list)))

#### Access Road Distance Models ####

Bio_Com_X_Turbine_Rd_Dist <- occu( ~ max_trig_dist + NDVI_1_9km 
                                   ~ turbine_rd_dist * as.factor(biotic_com_2) +
                                     shrub_yucca_density, occu.odhe)

Vis_Obs_X_Turbine_Rd_Dist <- occu( ~ max_trig_dist + NDVI_1_9km 
                                   ~ veg_cover * turbine_rd_dist +
                                     shrub_yucca_density, occu.odhe)

Habitat_Add_Turbine_Rd_Dist <- occu( ~ max_trig_dist + NDVI_1_9km 
                                     ~ shrub_yucca_density + 
                                       veg_cover + turbine_rd_dist, 
                                       occu.odhe)

Turbine_Rd_Dist_Occu <- occu( ~ max_trig_dist + NDVI_1_9km ~ 
                              ~ turbine_rd_dist, occu.odhe)

Vis_Obs_Add_Turbine_Rd_Dist <- occu( ~ max_trig_dist + NDVI_1_9km 
                                     ~ turbine_rd_dist + veg_cover,
                                       occu.odhe)

Shrub_Dense_Add_Turbine_Rd_Dist <- occu( ~ max_trig_dist + NDVI_1_9km 
                                         ~ turbine_rd_dist +
                                           shrub_yucca_density, occu.odhe)

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

# AICc Table with all possible access rd dist models 

model_list <- list(
  "Visual obstruction + shrub density" = Habitat,
  "Turbine road distance + visual obstruction + shrub density" = Habitat_Add_Turbine_Rd_Dist,
  "Turbine road distance x visual obstruction + shrub density" = Vis_Obs_X_Turbine_Rd_Dist,
  "Turbine road distance x biotic community + shrub density" = Bio_Com_X_Turbine_Rd_Dist,
  "Null" = Null,
  "Turbine road distance" = Turbine_Rd_Dist_Occu,
  "Visual obstruction" = Vis_Obs_Occu,
  "Shrub density" = Shrub_Dense_Occu,
  "Biotic community + shrub density" = Bio_Com_Shrub,
  "Turbine road distance + visual obstruction" = Vis_Obs_Add_Turbine_Rd_Dist,
  "Turbine road distance + shrub density" = Shrub_Dense_Add_Turbine_Rd_Dist
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine road distance" = "Null",
  "Turbine road distance + visual obstruction + shrub density" = 
    "Visual obstruction + shrub density",
  "Turbine road distance + shrub density" = "Shrub density",
  "Turbine road distance + visual obstruction" = "Visual obstruction",
  "Turbine road distance x visual obstruction + shrub density" = 
    "Visual obstruction + shrub density",
  "Turbine road distance x biotic community + shrub density" = 
    "Biotic community + shrub density"
)

# Ensure model_list matches wind_pairs
stopifnot(all(names(wind_pairs) %in% names(model_list)))
stopifnot(all(unlist(wind_pairs) %in% names(model_list)))

#### Access Road Density Models ####

Bio_Com_X_Turbine_Rd_Dense <- occu( ~ max_trig_dist + NDVI_1_9km 
                                    ~ turbine_rd_density_1_9km * 
                                      as.factor(biotic_com_2) +
                                      shrub_yucca_density, occu.odhe)

Vis_Obs_X_Turbine_Rd_Dense <- occu( ~ max_trig_dist + NDVI_1_9km ~ 
                                      ~ veg_cover *
                                      turbine_rd_density_1_9km +
                                      shrub_yucca_density, occu.odhe)

Habitat_Add_Turbine_Rd_Dense <- occu( ~ max_trig_dist + NDVI_1_9km 
                                      ~ shrub_yucca_density + 
                                        veg_cover +
                                        turbine_rd_density_1_9km, occu.odhe)

Turbine_Rd_Dense_Occu <- occu( ~ max_trig_dist + NDVI_1_9km 
                               ~ turbine_rd_density_1_9km, occu.odhe)

Vis_Obs_Add_Turbine_Rd_Dense <- occu( ~ max_trig_dist + NDVI_1_9km 
                                      ~ turbine_rd_density_1_9km + 
                                        veg_cover, occu.odhe)

Shrub_Rd_Dense_Add_Turbine_Rd_Dense <- occu( ~ max_trig_dist + NDVI_1_9km 
                                             ~ turbine_rd_density_1_9km +
                                               shrub_yucca_density, occu.odhe)

# 85% CI for bio community interaction 
# Extract coefficients and VCOV matrix
coef_est <- coef(Bio_Com_X_Turbine_Rd_Dense, type = "state")
vcov_mat <- vcov(Bio_Com_X_Turbine_Rd_Dense, type = "state")

beta_grassland <- coef_est["psi(turbine_rd_density_1_9km)"]
se_grassland <- sqrt(vcov_mat["psi(turbine_rd_density_1_9km)", 
                              "psi(turbine_rd_density_1_9km)"])

beta_woodland <- coef_est["psi(turbine_rd_density_1_9km)"] + 
  coef_est["psi(turbine_rd_density_1_9km:as.factor(biotic_com_2)2)"]

var_woodland <- 
  vcov_mat["psi(turbine_rd_density_1_9km)", "psi(turbine_rd_density_1_9km)"] +
  vcov_mat["psi(turbine_rd_density_1_9km:as.factor(biotic_com_2)2)",
           "psi(turbine_rd_density_1_9km:as.factor(biotic_com_2)2)"] +
  2 * vcov_mat["psi(turbine_rd_density_1_9km)", 
               "psi(turbine_rd_density_1_9km:as.factor(biotic_com_2)2)"]
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
  "Visual obstruction + shrub density" = Habitat,
  "Turbine road density + visual obstruction + shrub density" = Habitat_Add_Turbine_Rd_Dense,
  "Turbine road density x visual obstruction + shrub density" = Vis_Obs_X_Turbine_Rd_Dense,
  "Turbine road density x biotic community + shrub density" = Bio_Com_X_Turbine_Rd_Dense,
  "Null" = Null,
  "Turbine road density" = Turbine_Rd_Dense_Occu,
  "Visual obstruction" = Vis_Obs_Occu,
  "Shrub density" = Shrub_Dense_Occu,
  "Biotic community + shrub density" = Bio_Com_Shrub,
  "Turbine road density + visual obstruction" = Vis_Obs_Add_Turbine_Rd_Dense,
  "Turbine road density + shrub density" = Shrub_Rd_Dense_Add_Turbine_Rd_Dense
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine road density" = "Null",
  "Turbine road density + visual obstruction + shrub density" = 
    "Visual obstruction + shrub density",
  "Turbine road density + shrub density" = "Shrub density",
  "Turbine road density + visual obstruction" = "Visual obstruction",
  "Turbine road density x visual obstruction + shrub density" = 
    "Visual obstruction + shrub density",
  "Turbine road density x biotic community + shrub density" = 
    "Biotic community + shrub density"
)

# Ensure model_list matches wind_pairs
stopifnot(all(names(wind_pairs) %in% names(model_list)))
stopifnot(all(unlist(wind_pairs) %in% names(model_list)))

#### Wind Variable Additive Mods ####

Additive1 <- occu( ~ max_trig_dist + NDVI_1_9km
                   ~ as.factor(biotic_com_2) * X150cm_turbine_vis +
                     shrub_yucca_density + turbine_dist, occu.odhe)  

Additive2 <- occu( ~ max_trig_dist + NDVI_1_9km 
                   ~ X150cm_turbine_vis + shrub_yucca_density + turbine_dist + 
                     veg_cover, occu.odhe)  

Additive3 <- occu( ~ max_trig_dist + NDVI_1_9km 
                   ~ X150cm_turbine_vis + shrub_yucca_density + turbine_dist, 
                     occu.odhe)  

Additive4 <- occu( ~ max_trig_dist + NDVI_1_9km 
                   ~ X150cm_turbine_vis + veg_cover + turbine_dist, 
                     occu.odhe)  

Additive5 <- occu( ~ max_trig_dist + NDVI_1_9km 
                   ~ X150cm_turbine_vis + turbine_dist, occu.odhe)  

# biotic com x turbine vis 85% CI
# Extract coefficients and VCOV matrix
coef_est <- coef(Additive1, type = "state")
vcov_mat <- vcov(Additive1, type = "state")

beta_grassland <- coef_est["psi(X150cm_turbine_vis)"]
se_grassland <- sqrt(vcov_mat["psi(X150cm_turbine_vis)", 
                              "psi(X150cm_turbine_vis)"])

beta_woodland <- coef_est["psi(X150cm_turbine_vis)"] + 
  coef_est["psi(as.factor(biotic_com_2)2:X150cm_turbine_vis)"]

var_woodland <- 
  vcov_mat["psi(X150cm_turbine_vis)", "psi(X150cm_turbine_vis)"] +
  vcov_mat["psi(as.factor(biotic_com_2)2:X150cm_turbine_vis)",
           "psi(as.factor(biotic_com_2)2:X150cm_turbine_vis)"] +
  2 * vcov_mat["psi(X150cm_turbine_vis)", 
               "psi(as.factor(biotic_com_2)2:X150cm_turbine_vis)"]
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

# AICc Table with all possible additive models 

model_list <- list(
  "Turbine visibility x biotic community + shrub density + turbine distance" = Additive1, 
  "Turbine visibility + shrub density + turbine distance + visual obstruction" =
    Additive2, 
  "Turbine visibility + shrub density + turbine distance" = Additive3, 
  "Turbine visibility + visual obstruction + turbine distance" = Additive4,
  "Turbine visibility + turbine distance" = Additive5, 
  "Null Model" = Null,
  "Visual obstruction" = Vis_Obs_Occu,
  "Shrub density" = Shrub_Dense_Occu,
  "Visual obstruction + shrub density" = Habitat,
  "Biotic community + shrub density" = Bio_Com_Shrub,
  "shrub density + visual obstruction" = Habitat 
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine visibility x biotic community + shrub density + turbine distance" =
    "Biotic community + shrub density",
  "Turbine visibility + shrub density + turbine distance + visual obstruction" =
    "Visual obstruction + shrub density", 
  "Turbine visibility + shrub density + turbine distance" = "Shrub density", 
  "Turbine visibility + visual obstruction + turbine distance" = 
    "Visual obstruction",
  "Turbine visibility + turbine distance" = "Null Model"
)

# Ensure model_list matches wind_pairs
stopifnot(all(names(wind_pairs) %in% names(model_list)))
stopifnot(all(unlist(wind_pairs) %in% names(model_list)))

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
  ER <- if(all(c(mod_with_wind, mod_without) 
               %in% names(w))) w[mod_with_wind] / w[mod_without] else NA
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

#write.csv(final_table, "mule_deer_turbine_interior.csv", row.names = FALSE)

#write.csv(final_table, "mule_deer_turbine_vis.csv", row.names = FALSE)

#write.csv(final_table, "mule_deer_turbine_dist.csv", row.names = FALSE)

#write.csv(final_table, "mule_deer_turbine_density.csv", row.names = FALSE)

#write.csv(final_table, "mule_deer_turbine_rd_dist.csv", row.names = FALSE)

#write.csv(final_table, "mule_deer_turbine_rd_density.csv", row.names = FALSE)

#write.csv(final_table, "mule_deer_additive_mods.csv", row.names = FALSE)

########################################################
######################### END ##########################
########################################################
