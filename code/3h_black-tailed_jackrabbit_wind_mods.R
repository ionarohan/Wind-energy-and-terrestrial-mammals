###########################################################
#### Code for: base models for black-tailed jackrabbits ###
###########################################################
################ Script by: Iona Rohan ####################
############ Contact: ionarohan12@gmail.com ###############
###########################################################
########## Date Last Modified: 19-March-2026 ##############
###########################################################

###########################################################
###### NOTE:RUN THIS CODE AFTER "1_pre_model_code.R" ######
###########################################################

#Clear work environment
rm(list=ls())

#Note: If you opened this script through the .Rproj file, the only line you 
  #should need to change for the script to run (assuming packages are installed) 
  #is the homewd directory on line 23.

#Set home working directory
  #e.g. homewd = "C:/Users/ionar/Desktop/R Repository/Wind-energy-and-terrestrial-mammals/"
homewd = "<insert your folder here and end with a forward slash>"

#Set wd to data folder on your local computer 
setwd(paste0(homewd, "data/"))

###########################################################
# SETUP CODE FOR BLACK-TAILED JACKRABBITS OCCUPANCY MODELS #
###########################################################

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

#################################################
# CREATE THE WIND MODELS #
#################################################

#### Models with no change over all turbine variables ####

Habitat <- occu( ~  as.factor(precip.cat) + veg_cover_cam_under_1m
                 ~  woodland_percent_0_9km + slope, occu.leca)

Null <- occu( ~  as.factor(precip.cat) + veg_cover_cam_under_1m
              ~ 1, occu.leca)

Slope_Occu <- occu( ~  as.factor(precip.cat) + veg_cover_cam_under_1m
                    ~ slope, occu.leca)

Wood_Occu <- occu( ~  as.factor(precip.cat) + veg_cover_cam_under_1m
                   ~  woodland_percent_0_9km, occu.leca)

Bio_Com_Occu <- occu( ~  as.factor(precip.cat) + 
                         veg_cover_cam_under_1m
                      ~  as.factor(biotic_com_2) + slope, occu.leca)

Canopy_Occu <- occu( ~  as.factor(precip.cat) +
                        veg_cover_cam_under_1m
                     ~  canopy_cov + slope, occu.leca)

## The table with the comparison of the relative weight of evidence between
  #occupancy models separately ranked for the effect of each wind energy 
  #variable on the probability of habitat selection (ψ) of black-tailed 
  #jackrabbits is found in Table S3.20.

#### Turbine Interior Models ####

Turbine_Int_Occu <- occu( ~  as.factor(precip.cat) +
                             veg_cover_cam_under_1m
                          ~  as.factor(turbine_interior), occu.leca)

Habitat_Add_Turbine_Int <-  occu( ~  as.factor(precip.cat) +
                                     veg_cover_cam_under_1m
                                  ~  woodland_percent_0_9km + slope +
                                     as.factor(turbine_interior), 
                                     occu.leca)

Wood_Add_Turbine_Int <-  occu( ~  as.factor(precip.cat) + 
                                  veg_cover_cam_under_1m
                               ~  woodland_percent_0_9km +
                                  as.factor(turbine_interior), occu.leca)

Slope_Add_Turbine_Int <- occu( ~ as.factor(precip.cat) + 
                                 veg_cover_cam_under_1m
                               ~ slope + as.factor(turbine_interior), 
                                 occu.leca)

Bio_Com_X_Turbine_Int <- occu( ~  as.factor(precip.cat) +
                                  veg_cover_cam_under_1m
                               ~  slope + as.factor(turbine_interior) * 
                                  as.factor(biotic_com_2), occu.leca)

Wood_X_Turbine_Int <- occu( ~  as.factor(precip.cat) +
                               veg_cover_cam_under_1m
                            ~  slope + as.factor(turbine_interior) *
                               woodland_percent_0_9km, occu.leca)

# 85% CI for woodland interaction 
# Extract coefficients and VCOV matrix
coef_est <- coef(Wood_X_Turbine_Int, type = "state")
vcov_mat <- vcov(Wood_X_Turbine_Int, type = "state")

# woodland effect OUTSIDE wind farm (reference = 0) 
beta_outside <- coef_est["psi(woodland_percent_0_9km)"]

se_outside <- sqrt(
  vcov_mat["psi(woodland_percent_0_9km)", "psi(woodland_percent_0_9km)"]
)

# woodland effect INSIDE wind farm (reference + interaction)
beta_inside <- coef_est["psi(woodland_percent_0_9km)"] +
  coef_est["psi(as.factor(turbine_interior)1:woodland_percent_0_9km)"]

var_inside <-
  vcov_mat["psi(woodland_percent_0_9km)", "psi(woodland_percent_0_9km)"] +
  vcov_mat["psi(as.factor(turbine_interior)1:woodland_percent_0_9km)",
           "psi(as.factor(turbine_interior)1:woodland_percent_0_9km)"] +
  2 * vcov_mat["psi(woodland_percent_0_9km)",
               "psi(as.factor(turbine_interior)1:woodland_percent_0_9km)"]

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

# 85% CI for bio community interaction 
# Extract coefficients and VCOV matrix
coef_est <- coef(Bio_Com_X_Turbine_Int, type = "state")
vcov_mat <- vcov(Bio_Com_X_Turbine_Int, type = "state")

# Exact coefficient names from your model
intercept_name   <- "psi(Int)"
turbine_name     <- "psi(as.factor(turbine_interior)1)"
habitat_name     <- "psi(as.factor(biotic_com_2)2)"
interaction_name <- 
  "psi(as.factor(turbine_interior)1:as.factor(biotic_com_2)2)"
slope            <- "psi(slope)"


# Grassland, outside wind farm (reference) 
beta_grassland <- coef_est[intercept_name]
se_grassland <- sqrt(vcov_mat[intercept_name, intercept_name])

# Grassland, inside wind farm
beta_grassland_inside <- coef_est[intercept_name] + coef_est[turbine_name]
var_grassland_inside <- vcov_mat[intercept_name, intercept_name] +
  vcov_mat[turbine_name, turbine_name] +
  2 * vcov_mat[intercept_name, turbine_name]
se_grassland_inside <- sqrt(var_grassland_inside)

# ---- Woodland, outside wind farm 
beta_woodland <- coef_est[intercept_name] + coef_est[habitat_name]
var_woodland <- vcov_mat[intercept_name, intercept_name] +
  vcov_mat[habitat_name, habitat_name] +
  2 * vcov_mat[intercept_name, habitat_name]
se_woodland <- sqrt(var_woodland)

# ---- Woodland, inside wind farm 
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

# ---- z-scores for 85% CI
z <- qnorm(c(0.075, 0.925))

# ---- Confidence intervals 
ci_grassland_out <- beta_grassland + z * se_grassland
ci_grassland_in  <- beta_grassland_inside + z * se_grassland_inside
ci_woodland_out  <- beta_woodland + z * se_woodland
ci_woodland_in   <- beta_woodland_inside + z * se_woodland_inside

# ---- Tidy output 
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
  "Percent wooodland + slope" = Habitat,
  "Null" = Null,
  "Biotic community + slope" = Bio_Com_Occu,
  "Slope" = Slope_Occu,
  "Percent woodland" = Wood_Occu,
  "Turbine interior" = Turbine_Int_Occu,
  "Turbine interior + percent wooodland + slope" = Habitat_Add_Turbine_Int,
  "Turbine interior + percent wooodland" = Wood_Add_Turbine_Int,
  "Turbine interior + slope" = Slope_Add_Turbine_Int,
  "Turbine interior x percent wooodland + slope" = Wood_X_Turbine_Int,
  "Turbine interior x biotic community + slope" = Bio_Com_X_Turbine_Int
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine interior" = "Null",
  "Turbine interior + percent wooodland + slope" = "Percent wooodland + slope",
  "Turbine interior + percent wooodland" = "Percent woodland",
  "Turbine interior + slope" = "Slope",
  "Turbine interior x percent wooodland + slope" = "Percent wooodland + slope",
  "Turbine interior x biotic community + slope" = "Biotic community + slope"
)

# Ensure model_list matches wind_pairs
stopifnot(all(names(wind_pairs) %in% names(model_list)))
stopifnot(all(unlist(wind_pairs) %in% names(model_list)))

# The beta values (ß), standard errors, and 85% confidence intervals for 
  #parameter estimates within the best-supported model describing the effect of
  #the interior of the wind farm on the probability of habitat selection (ψ) for 
  #black-tailed jackrabbits are listed in Table S3.32.

#### Turbine Visibility Models ####

Turbine_Vis_Occu <- occu( ~ as.factor(precip.cat) +
                            veg_cover_cam_under_1m
                          ~ X50cm_turbine_vis, occu.leca)

Habitat_Add_Turbine_Vis <- occu( ~  as.factor(precip.cat) +
                                   veg_cover_cam_under_1m
                                 ~  woodland_percent_0_9km + slope +
                                   X50cm_turbine_vis, occu.leca)

Wood_Add_Turbine_Vis <- occu( ~  as.factor(precip.cat) + 
                                veg_cover_cam_under_1m
                              ~  woodland_percent_0_9km +
                                X50cm_turbine_vis, occu.leca)

Slope_Add_Turbine_Vis <- occu( ~  as.factor(precip.cat) + 
                                 veg_cover_cam_under_1m
                               ~  slope + X50cm_turbine_vis, occu.leca)

Bio_Com_X_Turbine_Vis <- occu( ~  as.factor(precip.cat) +
                                 veg_cover_cam_under_1m
                               ~  slope + X50cm_turbine_vis * 
                                 as.factor(biotic_com_2), occu.leca)

Canopy_X_Turbine_Vis <- occu( ~  as.factor(precip.cat) +
                                veg_cover_cam_under_1m
                              ~  slope + X50cm_turbine_vis * 
                                canopy_cov, occu.leca)

Wood_X_Turbine_Vis <- occu( ~ as.factor(precip.cat) +
                              veg_cover_cam_under_1m
                            ~ slope + X50cm_turbine_vis * 
                              woodland_percent_0_9km, occu.leca)

# 85% CI for bio community interaction 
# Extract coefficients and VCOV matrix
coef_est <- coef(Bio_Com_X_Turbine_Vis, type = "state")
vcov_mat <- vcov(Bio_Com_X_Turbine_Vis, type = "state")

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
  "Percent wooodland + slope" = Habitat,
  "Null" = Null,
  "Biotic community + slope" = Bio_Com_Occu,
  "Canopy cover + slope" = Canopy_Occu,
  "Slope" = Slope_Occu,
  "Percent woodland" = Wood_Occu,
  "Turbine visibility" = Turbine_Vis_Occu,
  "Turbine visibility + percent wooodland + slope" = Habitat_Add_Turbine_Vis,
  "Turbine visibility + percent wooodland" = Wood_Add_Turbine_Vis,
  "Turbine visibility + slope" = Slope_Add_Turbine_Vis,
  "Turbine visibility x percent wooodland + slope" = Wood_X_Turbine_Vis,
  "Turbine visibility x biotic community + slope" = Bio_Com_X_Turbine_Vis,
  "Turbine visibility x canopy cover + slope" = Canopy_X_Turbine_Vis
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine visibility" = "Null",
  "Turbine visibility + percent wooodland + slope" = 
    "Percent wooodland + slope",
  "Turbine visibility + percent wooodland" = "Percent woodland",
  "Turbine visibility + slope" = "Slope",
  "Turbine visibility x percent wooodland + slope" = 
    "Percent wooodland + slope",
  "Turbine visibility x biotic community + slope" = "Biotic community + slope",
  "Turbine visibility x canopy cover + slope" = "Canopy cover + slope" 
)

# Ensure model_list matches wind_pairs
stopifnot(all(names(wind_pairs) %in% names(model_list)))
stopifnot(all(unlist(wind_pairs) %in% names(model_list)))

# The beta values (ß), standard errors, and 85% confidence intervals for 
  #parameter estimates within the best-supported model describing the effect of
  #the visibility of wind turbines from 50 cm off the ground on the probability 
  #of habitat selection (ψ) for black-tailed jackrabbits are listed in 
  #Table S3.33.

#### Turbine Distance Models ####

Turbine_Dist_Occu <- occu( ~  as.factor(precip.cat) +
                              veg_cover_cam_under_1m
                           ~  turbine_dist, occu.leca)

Habitat_Add_Turbine_Dist <- occu( ~  as.factor(precip.cat) + 
                                     veg_cover_cam_under_1m
                                  ~  woodland_percent_0_9km + slope + 
                                     turbine_dist, occu.leca)

Wood_Add_Turbine_Dist <- occu( ~  as.factor(precip.cat) +
                                  veg_cover_cam_under_1m
                               ~  woodland_percent_0_9km + turbine_dist, 
                                  occu.leca)

Slope_Add_Turbine_Dist <- occu( ~  as.factor(precip.cat) + 
                                   veg_cover_cam_under_1m
                                ~  slope + turbine_dist, occu.leca)

Bio_Com_X_Turbine_Dist <- occu( ~  as.factor(precip.cat) + 
                                   veg_cover_cam_under_1m
                                ~  slope + turbine_dist * 
                                   as.factor(biotic_com_2), occu.leca)

Wood_X_Turbine_Dist <- occu( ~  as.factor(precip.cat) + 
                                veg_cover_cam_under_1m
                             ~  slope + woodland_percent_0_9km *
                                turbine_dist, occu.leca)

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
  "Percent wooodland + slope" = Habitat,
  "Null" = Null,
  "Biotic community + slope" = Bio_Com_Occu,
  "Slope" = Slope_Occu,
  "Percent woodland" = Wood_Occu,
  "Turbine distance" = Turbine_Dist_Occu,
  "Turbine distance + percent wooodland + slope" = Habitat_Add_Turbine_Dist,
  "Turbine distance + percent wooodland" = Wood_Add_Turbine_Dist,
  "Turbine distance + slope" = Slope_Add_Turbine_Dist,
  "Turbine distance x percent wooodland + slope" = Wood_X_Turbine_Dist,
  "Turbine distance x biotic community + slope" = Bio_Com_X_Turbine_Dist
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine distance" = "Null",
  "Turbine distance + percent wooodland + slope" = 
    "Percent wooodland + slope",
  "Turbine distance + percent wooodland" = "Percent woodland",
  "Turbine distance + slope" = "Slope",
  "Turbine distance x percent wooodland + slope" = 
    "Percent wooodland + slope",
  "Turbine distance x biotic community + slope" = "Biotic community + slope"
)

# Ensure model_list matches wind_pairs
stopifnot(all(names(wind_pairs) %in% names(model_list)))
stopifnot(all(unlist(wind_pairs) %in% names(model_list)))

#### Turbine Density Models ####

Turbine_Dense_Occu <- occu( ~  as.factor(precip.cat) + 
                               veg_cover_cam_under_1m
                            ~  turbine_density_0_9km, occu.leca)

Habitat_Add_Turbine_Dense <- occu( ~  as.factor(precip.cat) +
                                      veg_cover_cam_under_1m
                                   ~  woodland_percent_0_9km + slope +
                                      turbine_density_0_9km, occu.leca)

Wood_Add_Turbine_Dense <- occu( ~  as.factor(precip.cat) +
                                   veg_cover_cam_under_1m
                                ~  woodland_percent_0_9km + 
                                   turbine_density_0_9km, occu.leca)

Slope_Add_Turbine_Dense <- occu( ~  as.factor(precip.cat) +
                                    veg_cover_cam_under_1m
                                ~   slope + turbine_density_0_9km, occu.leca)

Bio_Com_X_Turbine_Dense <- occu( ~  as.factor(precip.cat) +
                                    veg_cover_cam_under_1m
                                 ~  slope + turbine_density_0_9km * 
                                    as.factor(biotic_com_2), occu.leca)

Wood_X_Turbine_Dense <- occu( ~  as.factor(precip.cat) + 
                                 veg_cover_cam_under_1m
                              ~  slope + woodland_percent_0_9km *
                                 turbine_density_0_9km, occu.leca)

# 85% CI for bio community interaction 
# Extract coefficients and VCOV matrix
coef_est <- coef(Bio_Com_X_Turbine_Dense, type = "state")
vcov_mat <- vcov(Bio_Com_X_Turbine_Dense, type = "state")

beta_grassland <- coef_est["psi(turbine_density_0_9km)"]
se_grassland <- sqrt(vcov_mat["psi(turbine_density_0_9km)", 
                              "psi(turbine_density_0_9km)"])

beta_woodland <- coef_est["psi(turbine_density_0_9km)"] + 
  coef_est["psi(turbine_density_0_9km:as.factor(biotic_com_2)2)"]

var_woodland <- 
  vcov_mat["psi(turbine_density_0_9km)", "psi(turbine_density_0_9km)"] +
  vcov_mat["psi(turbine_density_0_9km:as.factor(biotic_com_2)2)",
           "psi(turbine_density_0_9km:as.factor(biotic_com_2)2)"] +
  2 * vcov_mat["psi(turbine_density_0_9km)", 
               "psi(turbine_density_0_9km:as.factor(biotic_com_2)2)"]
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
  "Percent wooodland + slope" = Habitat,
  "Base Detection Submodel" = Null,
  "Biotic community + slope" = Bio_Com_Occu,
  "Slope" = Slope_Occu,
  "Percent woodland" = Wood_Occu,
  "Turbine density" = Turbine_Dense_Occu,
  "Turbine density + percent wooodland + slope" = Habitat_Add_Turbine_Dense,
  "Turbine density + percent wooodland" = Wood_Add_Turbine_Dense,
  "Turbine density + slope" = Slope_Add_Turbine_Dense,
  "Turbine density x percent wooodland + slope" = Wood_X_Turbine_Dense,
  "Turbine density x biotic community + slope" = Bio_Com_X_Turbine_Dense
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine density" = "Base Detection Submodel",
  "Turbine density + percent wooodland + slope" = 
    "Percent wooodland + slope",
  "Turbine density + percent wooodland" = "Percent woodland",
  "Turbine density + slope" = "Slope",
  "Turbine density x percent wooodland + slope" = 
    "Percent wooodland + slope",
  "Turbine density x biotic community + slope" = "Biotic community + slope"
)

# Ensure model_list matches wind_pairs
stopifnot(all(names(wind_pairs) %in% names(model_list)))
stopifnot(all(unlist(wind_pairs) %in% names(model_list)))

#### Turbine Road Distance Models ####

Turbine_Rd_Dist_Occu <- occu( ~  as.factor(precip.cat) + 
                                veg_cover_cam_under_1m
                              ~  turbine_rd_dist, occu.leca)

Habitat_Add_Turbine_Rd_Dist <- occu( ~  as.factor(precip.cat) +
                                       veg_cover_cam_under_1m
                                     ~  woodland_percent_0_9km + slope +
                                       turbine_rd_dist, occu.leca)

Wood_Add_Turbine_Rd_Dist <- occu( ~  as.factor(precip.cat) +
                                    veg_cover_cam_under_1m
                                  ~  woodland_percent_0_9km +
                                    turbine_rd_dist, occu.leca)

Slope_Add_Turbine_Rd_Dist <- occu( ~  as.factor(precip.cat) +
                                     veg_cover_cam_under_1m
                                   ~  slope + turbine_rd_dist, occu.leca)

Bio_Com_X_Turbine_Rd_Dist <- occu( ~ as.factor(precip.cat) +
                                     veg_cover_cam_under_1m
                                   ~ slope + turbine_rd_dist * 
                                     as.factor(biotic_com_2), occu.leca)

Wood_X_Turbine_Rd_Dist <- occu( ~  as.factor(precip.cat) + 
                                  veg_cover_cam_under_1m
                                ~   slope + woodland_percent_0_9km *
                                  turbine_rd_dist, occu.leca)

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

# AICc Table with all possible turbine rd distance models 

model_list <- list(
  "Percent wooodland + slope" = Habitat,
  "Null" = Null,
  "Biotic community + slope" = Bio_Com_Occu,
  "Slope" = Slope_Occu,
  "Percent woodland" = Wood_Occu,
  "Turbine road distance" = Turbine_Rd_Dist_Occu,
  "Turbine road distance + percent wooodland + slope" = Habitat_Add_Turbine_Rd_Dist,
  "Turbine road distance + percent wooodland" = Wood_Add_Turbine_Rd_Dist,
  "Turbine road distance + slope" = Slope_Add_Turbine_Rd_Dist,
  "Turbine road distance x percent wooodland + slope" = Wood_X_Turbine_Rd_Dist,
  "Turbine road distance x biotic community + slope" = Bio_Com_X_Turbine_Rd_Dist 
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine road distance" = "Null",
  "Turbine road distance + percent wooodland + slope" = 
    "Percent wooodland + slope",
  "Turbine road distance + percent wooodland" = "Percent woodland",
  "Turbine road distance + slope" = "Slope",
  "Turbine road distance x percent wooodland + slope" = 
    "Percent wooodland + slope",
  "Turbine road distance x biotic community + slope" = "Biotic community + slope"
)

# Ensure model_list matches wind_pairs
stopifnot(all(names(wind_pairs) %in% names(model_list)))
stopifnot(all(unlist(wind_pairs) %in% names(model_list)))

# The beta values (ß), standard errors, and 85% confidence intervals for 
  #parameter estimates within the best-supported model describing the effect of
  #the distance to the nearest turbine access road on the probability 
  #of habitat selection (ψ) for black-tailed jackrabbits are listed in 
  #Table S3.33.

#### Turbine Road Density Models ####

Turbine_Rd_Dense_Occu <- occu( ~  as.factor(precip.cat) + 
                                  veg_cover_cam_under_1m
                               ~  turbine_rd_density_0_9km, occu.leca)

Habitat_Add_Turbine_Rd_Dense <- occu( ~  as.factor(precip.cat) +
                                         veg_cover_cam_under_1m
                                      ~  woodland_percent_0_9km + slope +
                                         turbine_rd_density_0_9km, occu.leca)

Wood_Add_Turbine_Rd_Dense <- occu( ~  as.factor(precip.cat) + 
                                      veg_cover_cam_under_1m
                                   ~  woodland_percent_0_9km +
                                      turbine_rd_density_0_9km, occu.leca)

Slope_Add_Turbine_Rd_Dense <- occu( ~  as.factor(precip.cat) + 
                                       veg_cover_cam_under_1m
                                    ~  slope + turbine_rd_density_0_9km,
                                       occu.leca)

Bio_Com_X_Turbine_Rd_Dense <- occu( ~  as.factor(precip.cat) + 
                                       veg_cover_cam_under_1m
                                    ~  slope + turbine_rd_density_0_9km * 
                                       as.factor(biotic_com_2), occu.leca)

Wood_X_Turbine_Rd_Dense <- occu( ~  as.factor(precip.cat) + 
                                    veg_cover_cam_under_1m
                                 ~  slope + woodland_percent_0_9km *
                                    turbine_rd_density_0_9km, occu.leca)

# 85% CI for bio community interaction 
# Extract coefficients and VCOV matrix
coef_est <- coef(Bio_Com_X_Turbine_Rd_Dense, type = "state")
vcov_mat <- vcov(Bio_Com_X_Turbine_Rd_Dense, type = "state")

beta_grassland <- coef_est["psi(turbine_rd_density_0_9km)"]
se_grassland <- sqrt(vcov_mat["psi(turbine_rd_density_0_9km)", 
                              "psi(turbine_rd_density_0_9km)"])

beta_woodland <- coef_est["psi(turbine_rd_density_0_9km)"] + 
  coef_est["psi(turbine_rd_density_0_9km:as.factor(biotic_com_2)2)"]

var_woodland <- 
  vcov_mat["psi(turbine_rd_density_0_9km)", "psi(turbine_rd_density_0_9km)"] +
  vcov_mat["psi(turbine_rd_density_0_9km:as.factor(biotic_com_2)2)",
           "psi(turbine_rd_density_0_9km:as.factor(biotic_com_2)2)"] +
  2 * vcov_mat["psi(turbine_rd_density_0_9km)", 
               "psi(turbine_rd_density_0_9km:as.factor(biotic_com_2)2)"]
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

# AICc Table with all possible turbine rd density models 

model_list <- list(
  "Percent wooodland + slope" = Habitat,
  "Base Detection Submodel" = Null,
  "Biotic community + slope" = Bio_Com_Occu,
  "Slope" = Slope_Occu,
  "Percent woodland" = Wood_Occu,
  "Turbine road density" = Turbine_Rd_Dense_Occu,
  "Turbine road density + percent wooodland + slope" = Habitat_Add_Turbine_Rd_Dense,
  "Turbine road density + percent wooodland" = Wood_Add_Turbine_Rd_Dense,
  "Turbine road density + slope" = Slope_Add_Turbine_Rd_Dense,
  "Turbine road density x percent wooodland + slope" = Wood_X_Turbine_Rd_Dense,
  "Turbine road density x biotic community + slope" = Bio_Com_X_Turbine_Rd_Dense
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine road density" = "Base Detection Submodel",
  "Turbine road density + percent wooodland + slope" = 
    "Percent wooodland + slope",
  "Turbine road density + percent wooodland" = "Percent woodland",
  "Turbine road density + slope" = "Slope",
  "Turbine road density x percent wooodland + slope" = 
    "Percent wooodland + slope",
  "Turbine road density x biotic community + slope" = "Biotic community + slope"
)

# Ensure model_list matches wind_pairs
stopifnot(all(names(wind_pairs) %in% names(model_list)))
stopifnot(all(unlist(wind_pairs) %in% names(model_list)))

#### Additive Wind Variable Models ####

Rd_Dist_Vis_Habitat <- occu( ~ as.factor(precip.cat) + 
                               veg_cover_cam_under_1m
                             ~ slope + woodland_percent_0_9km +
                               turbine_rd_dist + X50cm_turbine_vis, occu.leca)

Rd_Dist_Vis_Wood <- occu( ~  as.factor(precip.cat) +
                             veg_cover_cam_under_1m  
                          ~  woodland_percent_0_9km + turbine_rd_dist + 
                             X50cm_turbine_vis, occu.leca)

Rd_Dist_Vis_Slope <- occu( ~  as.factor(precip.cat) + 
                              veg_cover_cam_under_1m 
                           ~  slope + turbine_rd_dist + X50cm_turbine_vis, 
                              occu.leca)

Rd_Dist_Int_Habitat <- occu( ~ as.factor(precip.cat) + 
                               veg_cover_cam_under_1m
                             ~ slope + woodland_percent_0_9km +
                               turbine_rd_dist + 
                               as.factor(turbine_interior), occu.leca)

Rd_Dist_Int_Wood <- occu( ~  as.factor(precip.cat) + 
                             veg_cover_cam_under_1m
                          ~  woodland_percent_0_9km + turbine_rd_dist + 
                             as.factor(turbine_interior), occu.leca)

Rd_Dist_Int_Slope <- occu( ~  as.factor(precip.cat) + 
                              veg_cover_cam_under_1m 
                           ~  slope + turbine_rd_dist + 
                              as.factor(turbine_interior), occu.leca)

Vis_Int_Habitat <- occu( ~ as.factor(precip.cat) + 
                           veg_cover_cam_under_1m
                         ~ slope + woodland_percent_0_9km + X50cm_turbine_vis +
                           as.factor(turbine_interior), occu.leca)

Vis_Int_Wood <- occu( ~ as.factor(precip.cat) + 
                        veg_cover_cam_under_1m 
                      ~ woodland_percent_0_9km + X50cm_turbine_vis + 
                        as.factor(turbine_interior), occu.leca)

Vis_Int_Slope <- occu( ~ as.factor(precip.cat) + 
                         veg_cover_cam_under_1m 
                       ~ slope + X50cm_turbine_vis + 
                         as.factor(turbine_interior), occu.leca)

Vis_Int_Rd_Dist_Habitat <- occu( ~ as.factor(precip.cat) + 
                                   veg_cover_cam_under_1m
                                 ~ slope + woodland_percent_0_9km + 
                                   turbine_rd_dist + X50cm_turbine_vis + 
                                   as.factor(turbine_interior), occu.leca)

Vis_Int_Rd_Dist_Wood <- occu( ~  as.factor(precip.cat) + 
                                 veg_cover_cam_under_1m 
                              ~  woodland_percent_0_9km + turbine_rd_dist +
                                 X50cm_turbine_vis + 
                                 as.factor(turbine_interior), occu.leca)

Vis_Int_Rd_Dist_Slope <- occu( ~  as.factor(precip.cat) + 
                                  veg_cover_cam_under_1m  
                               ~  slope + turbine_rd_dist + X50cm_turbine_vis + 
                                  as.factor(turbine_interior), occu.leca)

# AICc Table with all possible additive models 

model_list <- list(
  "Null" = Null,
  "Percent wooodland + slope" = Habitat,
  "Slope" = Slope_Occu,
  "Percent woodland" = Wood_Occu,
   "Turbine interior + turbine visibility + percent wooodland + slope" = 
    Vis_Int_Habitat,
  "Turbine interior + access road distance + percent wooodland + slope" = 
    Rd_Dist_Int_Habitat,
  "Turbine interior + turbine visibility + percent wooodland" = Vis_Int_Wood,
  "Turbine interior + turbine visibility + slope" = Vis_Int_Slope,
 "Turbine interior + access road distance + percent wooodland" = 
   Rd_Dist_Int_Wood,
  "Turbine interior + access road distance + slope" = Rd_Dist_Int_Slope,
  "Turbine visibility + access road distance + percent wooodland + slope" =  
   Rd_Dist_Vis_Habitat,
  "Turbine visibility + access road distance + percent wooodland" =  
    Rd_Dist_Vis_Wood,
  "Turbine visibility + access road distance + slope" = Rd_Dist_Vis_Slope,
  "Turbine interior + turbine visibility + access road distance + 
    percent wooodland + slope" = Vis_Int_Rd_Dist_Habitat,
  "Turbine interior + turbine visibility + access road distance + 
   percent wooodland"  = Vis_Int_Rd_Dist_Wood,
  "Turbine interior + turbine visibility + access road distance + slope" =
    Vis_Int_Rd_Dist_Slope
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine interior + turbine visibility + percent wooodland + slope" = 
    "Percent wooodland + slope",
  "Turbine interior + access road distance + percent wooodland + slope" = 
    "Percent wooodland + slope",
  "Turbine interior + turbine visibility + percent wooodland" = 
    "Percent woodland",
  "Turbine interior + turbine visibility + slope" = "Slope",
  "Turbine interior + access road distance + percent wooodland" =
    "Percent woodland",
  "Turbine interior + access road distance + slope" = "Slope",
  "Turbine visibility + access road distance + percent wooodland + slope" =  
    "Percent wooodland + slope",
  "Turbine visibility + access road distance + percent wooodland" =  
    "Percent woodland",
  "Turbine visibility + access road distance + slope" = "Slope",
  "Turbine interior + turbine visibility + access road distance + 
    percent wooodland + slope" = "Percent wooodland + slope",
  "Turbine interior + turbine visibility + access road distance + 
   percent wooodland"  = "Percent woodland",
  "Turbine interior + turbine visibility + access road distance + slope" = 
    "Slope"
)

# Ensure model_list matches wind_pairs
stopifnot(all(names(wind_pairs) %in% names(model_list)))
stopifnot(all(unlist(wind_pairs) %in% names(model_list)))

# The beta values (ß), standard errors, and 85% confidence intervals for 
  #parameter estimates within the best-supported model describing the effect of
  #the visibility of wind turbines from 50 cm off the ground and the wind farm 
  #interior on the probability of habitat selection (ψ) for black-tailed 
  #jackrabbits are listed in Table S3.32 and the visibility of wind turbines
  #from 50 cm off the ground and the distance the the nearest access road on the
  #probability of habitat selection (ψ) for black-tailed jackrabbits are 
  #listed in Table S3.33.

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

#write.csv(final_table, "jack_turbine_interior.csv", row.names = FALSE)

#write.csv(final_table, "jack_turbine_vis.csv", row.names = FALSE)

#write.csv(final_table, "jack_turbine_dist.csv", row.names = FALSE)

#write.csv(final_table, "jack_turbine_density.csv", row.names = FALSE)

#write.csv(final_table, "jack_turbine_rd_dist.csv", row.names = FALSE)

#write.csv(final_table, "jack_turbine_rd_density.csv", row.names = FALSE)

#write.csv(final_table, "jack_turbine_additive_mods.csv", row.names = FALSE)

########################################################
######################### END ##########################
########################################################