###########################################################
########### Code for: wind models for mule deer ###########
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
# e.g. homewd = "C:/Users/ionar/Desktop/R Repository/Wind-energy-and-terrestrial-mammals/"
homewd = "<insert your folder here and end with a forward slash>"

#Set wd to data folder on your local computer 
setwd(paste0(homewd, "data/"))

# Load packages 
library(unmarked)
library(tidyverse)
library(MuMIn)
library(AICcmodavg)

###########################################################
# SETUP CODE FOR MULE DEER OCCUPANCY MODELS #
###########################################################

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

#####################################################################
# CREATE FUNCTION FOR MODEL SELECTION TABLES #
#####################################################################

occu_model_selection <- function(model_list, wind_pairs = NULL, 
                                 file_name = NULL) {
  
  # 1. Create AICc selection table
  
  top_mods <- aictab(
    cand.set = model_list,
    modnames = names(model_list)
  )
  
  top_mods_df <- as.data.frame(top_mods)
  top_mods_df$Model <- as.character(top_mods_df$Modnames)
  
  # Extract weights
  w <- top_mods_df$AICcWt
  names(w) <- top_mods_df$Model
  
  # 2. Create CI function
  
  get_85CI <- function(model, type = "state") {
    est <- try(coef(model, type = type), silent = TRUE)
    if(inherits(est, "try-error") || length(est) == 0) return(NULL)
    
    vc <- try(vcov(model, type = type), silent = TRUE)
    if(inherits(vc, "try-error")) return(NULL)
    
    se <- sqrt(diag(vc))
    z <- qnorm(0.925)
    
    data.frame(
      Parameter = names(est),
      CI_string = sprintf("%.2f [%.2f, %.2f]",
                          est, est - z*se, est + z*se)
    )
  }
  
  # 3. Compute CIs
  ci_psi_list <- lapply(model_list, get_85CI, type = "state")
  ci_det_list <- lapply(model_list, get_85CI, type = "det")
  
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
  
  if(ncol(ci_psi_wide) > 1) {
    names(ci_psi_wide)[-1] <- paste0(names(ci_psi_wide)[-1], "_psi")
  }
  if(ncol(ci_det_wide) > 1) {
    names(ci_det_wide)[-1] <- paste0(names(ci_det_wide)[-1], "_det")
  }
  
  # 4. Calculate LRTs 
  
  llr_df <- NULL
  
  if(!is.null(wind_pairs)) {
    llr_results <- lapply(names(wind_pairs), function(mod_with_wind) {
      mod_without <- wind_pairs[[mod_with_wind]]
      
      if(!(mod_with_wind %in% names(model_list)) ||
         !(mod_without %in% names(model_list))) {
        return(NULL)
      }
      
      m1 <- model_list[[mod_without]]
      m2 <- model_list[[mod_with_wind]]
      
      ll1 <- logLik(m1)
      ll2 <- logLik(m2)
      
      lrt <- 2 * (ll2 - ll1)
      
      k1 <- attr(ll1, "df")
      k2 <- attr(ll2, "df")
      
      if(is.null(k1) | is.null(k2)) {
        k1 <- length(coef(m1))
        k2 <- length(coef(m2))
      }
      
      df_diff <- k2 - k1
      pval <- pchisq(lrt, df = df_diff, lower.tail = FALSE)
      
      data.frame(
        Model = mod_with_wind,
        LRT_stat = as.numeric(lrt),
        LRT_df = df_diff,
        LRT_p = as.numeric(pval)
      )
    })
    
    llr_df <- bind_rows(llr_results)
  }
  
  # 5. Calculate evidence ratios
  
  pairwise_ER_df <- NULL
  
  if(!is.null(wind_pairs)) {
    pairwise_ER <- lapply(names(wind_pairs), function(mod_with_wind) {
      mod_without <- wind_pairs[[mod_with_wind]]
      
      ER <- if(all(c(mod_with_wind, mod_without) %in% names(w))) {
        w[mod_with_wind] / w[mod_without]
      } else NA
      
      data.frame(
        Model = mod_with_wind,
        Evidence_Ratio_vs_NonWind = ER
      )
    })
    
    pairwise_ER_df <- bind_rows(pairwise_ER)
  }
  
  # 6. Merge everything
  
  final_table <- top_mods_df %>%
    left_join(ci_psi_wide, by = "Model") %>%
    left_join(ci_det_wide, by = "Model") %>%
    left_join(llr_df, by = "Model") %>%
    left_join(pairwise_ER_df, by = "Model")
  
  # 7. Write CSV 
  if(!is.null(file_name)) {
    write.csv(final_table, file_name, row.names = FALSE)
  }
  
}

#################################################
# CREATE THE WIND MODELS #
#################################################

#### Models with no change over all turbine variables ####

Habitat <- occu( ~ max_trig_dist + livestock.count 
                 ~ shrub_yucca_density + veg_cover, occu.odhe,
                   starts = c(0, 2, 0, -2, -1, -2))

Veg_Cov_Occu <- occu( ~ max_trig_dist + livestock.count 
                      ~ veg_cover, occu.odhe)

Canopy_Shrub <- occu( ~ max_trig_dist + livestock.count 
                      ~ shrub_yucca_density + canopy_cov, occu.odhe)

Shrub_Dense_Occu <- occu( ~ max_trig_dist + livestock.count 
                          ~ shrub_yucca_density, occu.odhe)

Bio_Com_Occu <- occu( ~ max_trig_dist + livestock.count 
                      ~ as.factor(biotic_com_2), occu.odhe)

Bio_Com_Shrub <- occu( ~ max_trig_dist + livestock.count 
                       ~ as.factor(biotic_com_2) + shrub_yucca_density, 
                       occu.odhe)

Null <- occu( ~ max_trig_dist + livestock.count 
              ~ 1, occu.odhe, starts = c(0, -2, -1, -2))

## The table with the comparison of the relative weight of evidence between
#occupancy models separately ranked for the effect of each wind energy 
#variable on the probability of habitat selection (ψ) of mule deer is found 
#in Table S3.13.

## Check correlations between wind variables and habitat variables 

# Read in site-level covariates
site.covs <- read.csv("site_covs.csv", nrows = 102, header = TRUE)

# site-level variable correlations
wind.hab.cor <- site.covs %>% 
  select( biotic_com_2, shrub_yucca_density, veg_cover, canopy_cov,
          turbine_interior, X150cm_turbine_vis, turbine_dist,
          turbine_density_1_9km, turbine_rd_dist, turbine_rd_density_1_9km)

cors <- cor(wind.hab.cor, method='spearman')  
# none correlated above |0.7|

#### Turbine Interior Models ####

Bio_Com_X_Turbine_Int_Add_Shrub <- occu( ~ max_trig_dist + livestock.count 
                                         ~ as.factor(turbine_interior) *
                                           as.factor(biotic_com_2) +
                                           shrub_yucca_density, occu.odhe)

Veg_Cov_X_Turbine_Int_Add_Shrub <- occu( ~ max_trig_dist + livestock.count 
                                         ~ as.factor(turbine_interior) *
                                           veg_cover +
                                           shrub_yucca_density, occu.odhe)

Habitat_Add_Turbine_Int <- occu( ~ max_trig_dist + livestock.count  
                                 ~ shrub_yucca_density + veg_cover +
                                   as.factor(turbine_interior), occu.odhe,
                                   starts = c(1, 2, 0, -1, -2, -1, -2))

Turbine_Int_Occu <- occu( ~ max_trig_dist + livestock.count 
                          ~ as.factor(turbine_interior), occu.odhe)

Veg_Cov_Add_Turbine_Int <- occu( ~ max_trig_dist + livestock.count 
                                 ~ as.factor(turbine_interior) +
                                   veg_cover, occu.odhe)

Shrub_Dense_Add_Turbine_Int <- occu( ~ max_trig_dist + livestock.count 
                                     ~ as.factor(turbine_interior) +
                                       shrub_yucca_density, occu.odhe)   

# 85% CI for interior x veg cover
# Extract coefficients and VCOV matrix
coef_est <- coef(Veg_Cov_X_Turbine_Int_Add_Shrub, type = "state")
vcov_mat <- vcov(Veg_Cov_X_Turbine_Int_Add_Shrub, type = "state")

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

# 85% CI for interior x bio community
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
  "vegetation cover + shrub density" = Habitat,
  "Turbine interior + vegetation cover + shrub density" =
    Habitat_Add_Turbine_Int,
  "Turbine interior x vegetation cover + shrub density" =
    Veg_Cov_X_Turbine_Int_Add_Shrub,
  "Turbine interior x biotic community + shrub density" = 
    Bio_Com_X_Turbine_Int_Add_Shrub,
  "Null" = Null,
  "Turbine interior" = Turbine_Int_Occu,
  "vegetation cover" = Veg_Cov_Occu,
  "Shrub density" = Shrub_Dense_Occu,
  "Biotic community + shrub density" = Bio_Com_Shrub,
  "Turbine interior + vegetation cover" = Veg_Cov_Add_Turbine_Int,
  "Turbine interior + shrub density" = Shrub_Dense_Add_Turbine_Int
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine interior" = "Null",
  "Turbine interior + vegetation cover + shrub density" = 
    "vegetation cover + shrub density",
  "Turbine interior + shrub density" = "Shrub density",
  "Turbine interior + vegetation cover" = "vegetation cover",
  "Turbine interior x vegetation cover + shrub density" = 
    "vegetation cover + shrub density",
  "Turbine interior x biotic community + shrub density" = 
    "Biotic community + shrub density"
)

# Compile final table and export to .csv
final_table <- occu_model_selection(
  model_list = model_list,
  wind_pairs = wind_pairs,
  file_name = paste0(homewd, "outputs/mule_deer_turbine_interior.csv")
)

#### Turbine Visibility Models ####

Bio_Com_X_Turbine_Vis <- occu( ~ max_trig_dist + livestock.count 
                               ~ X150cm_turbine_vis * as.factor(biotic_com_2) +
                                 shrub_yucca_density, occu.odhe,
                                 starts = c(0,-1,5,2,5,-3,0,-2))
                                 # large SE for interaction term

Veg_Cov_X_Turbine_Vis <- occu( ~ max_trig_dist + livestock.count 
                               ~ veg_cover * X150cm_turbine_vis +
                                 shrub_yucca_density, occu.odhe)

Canopy_Cov_X_Turbine_Vis <- occu( ~ max_trig_dist + livestock.count 
                                  ~ canopy_cov * X150cm_turbine_vis +
                                    shrub_yucca_density, occu.odhe)

Habitat_Add_Turbine_Vis <- occu( ~ max_trig_dist + livestock.count 
                                 ~ shrub_yucca_density + veg_cover +
                                   X150cm_turbine_vis, occu.odhe)

Turbine_Vis_Occu <- occu( ~ max_trig_dist + livestock.count 
                          ~ X150cm_turbine_vis, occu.odhe)

Veg_Cov_Add_Turbine_Vis <- occu( ~ max_trig_dist + livestock.count 
                                 ~ X150cm_turbine_vis + veg_cover,
                                 occu.odhe)

Shrub_Dense_Add_Turbine_Vis <- occu( ~ max_trig_dist + livestock.count  
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
# woodland 85% CI very large

# AICc Table with all possible turbine visibility models 

model_list <- list(
  "vegetation cover + shrub density" = Habitat,
  "Turbine visibility + vegetation cover + shrub density" = 
    Habitat_Add_Turbine_Vis,
  "Turbine visibility x vegetation cover + shrub density" = 
    Veg_Cov_X_Turbine_Vis,
  "Turbine visibility x biotic community + shrub density" = 
    Bio_Com_X_Turbine_Vis,
  "Null" = Null,
  "Turbine visibility" = Turbine_Vis_Occu,
  "vegetation cover" = Veg_Cov_Occu,
  "Shrub density" = Shrub_Dense_Occu,
  "Biotic community + shrub density" = Bio_Com_Shrub,
  "Canopy cover + shrub density" = Canopy_Shrub,
  "Turbine visibility + shrub density" = Shrub_Dense_Add_Turbine_Vis,
  "Turbine visibility + vegetation cover" = Veg_Cov_Add_Turbine_Vis,
  "Turbine visibility x canopy cover + shrub density" = Canopy_Cov_X_Turbine_Vis
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine visibility" = "Null",
  "Turbine visibility + vegetation cover + shrub density" = 
    "vegetation cover + shrub density",
  "Turbine visibility + shrub density" = "Shrub density",
  "Turbine visibility + vegetation cover" = "vegetation cover",
  "Turbine visibility x vegetation cover + shrub density" = 
    "vegetation cover + shrub density",
  "Turbine visibility x biotic community + shrub density" = 
    "Biotic community + shrub density",
  "Turbine visibility x canopy cover + shrub density" = 
    "Canopy cover + shrub density"
)

# Compile final table and export to .csv
final_table <- occu_model_selection(
  model_list = model_list,
  wind_pairs = wind_pairs,
  file_name = paste0(homewd, "outputs/mule_deer_turbine_vis.csv")
)

# The beta values (ß), standard errors, and 85% confidence intervals for 
#parameter estimates within the best-supported model describing the effect of
#the visibility of wind turbines from 150 cm off the ground on the probability 
#of habitat selection (ψ) for mule deer are listed in Table S3.22.

#### Turbine Distance Models #####

Bio_Com_X_Turbine_Dist <- occu( ~ max_trig_dist + livestock.count  
                                ~ turbine_dist * as.factor(biotic_com_2) +
                                  shrub_yucca_density, occu.odhe)

Veg_Cov_X_Turbine_Dist <- occu( ~ max_trig_dist + livestock.count  
                                ~ veg_cover * turbine_dist +
                                  shrub_yucca_density, occu.odhe)

Habitat_Add_Turbine_Dist <- occu( ~ max_trig_dist + livestock.count 
                                  ~ shrub_yucca_density + veg_cover +
                                    turbine_dist, occu.odhe)

Turbine_Dist_Occu <- occu( ~ max_trig_dist + livestock.count 
                           ~ turbine_dist, occu.odhe)

Veg_Cov_Add_Turbine_Dist <- occu( ~ max_trig_dist + livestock.count 
                                  ~ turbine_dist + veg_cover, 
                                  occu.odhe)

Shrub_Dense_Add_Turbine_Dist <- occu( ~ max_trig_dist + livestock.count  
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
  "vegetation cover + shrub density" = Habitat,
  "Turbine distance + vegetation cover + shrub density" = Habitat_Add_Turbine_Dist,
  "Turbine distance x vegetation cover + shrub density" = Veg_Cov_X_Turbine_Dist,
  "Turbine distance x biotic community + shrub density" = Bio_Com_X_Turbine_Dist,
  "Null" = Null,
  "Turbine distance" = Turbine_Dist_Occu,
  "vegetation cover" = Veg_Cov_Occu,
  "Shrub density" = Shrub_Dense_Occu,
  "Biotic community + shrub density" = Bio_Com_Shrub,
  "Turbine distance + vegetation cover" = Veg_Cov_Add_Turbine_Dist,
  "Turbine distance + shrub density" = Shrub_Dense_Add_Turbine_Dist
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine distance" = "Null",
  "Turbine distance + vegetation cover + shrub density" = 
    "vegetation cover + shrub density",
  "Turbine distance + shrub density" = "Shrub density",
  "Turbine distance + vegetation cover" = "vegetation cover",
  "Turbine distance x vegetation cover + shrub density" = 
    "vegetation cover + shrub density",
  "Turbine distance x biotic community + shrub density" = 
    "Biotic community + shrub density"
)

# Compile final table and export to .csv
final_table <- occu_model_selection(
  model_list = model_list,
  wind_pairs = wind_pairs,
  file_name = paste0(homewd, "outputs/mule_deer_turbine_dist.csv")
)

# The beta values (ß), standard errors, and 85% confidence intervals for 
#parameter estimates within the best-supported model describing the effect of
#the distance to nearest wind turbine on the probability of habitat selection 
#(ψ) for mule deer are listed in Table S3.23.

#### Turbine Density Models ####

Bio_Com_X_Turbine_Dense <- occu( ~ max_trig_dist + livestock.count 
                                 ~ turbine_density_1_9km * 
                                   as.factor(biotic_com_2) +
                                   shrub_yucca_density, occu.odhe)

Veg_Cov_X_Turbine_Dense <- occu( ~ max_trig_dist + livestock.count 
                                 ~ veg_cover * turbine_density_1_9km +
                                   shrub_yucca_density, occu.odhe)

Habitat_Add_Turbine_Dense <- occu( ~ max_trig_dist + livestock.count 
                                   ~ shrub_yucca_density + veg_cover +
                                     turbine_density_1_9km, occu.odhe,
                                     starts = c(5,-5,-5,0,-3,0,-2))

Turbine_Dense_Occu <- occu( ~ max_trig_dist + livestock.count 
                            ~ turbine_density_1_9km, occu.odhe)

Veg_Cov_Add_Turbine_Dense <- occu( ~ max_trig_dist + livestock.count 
                                   ~ turbine_density_1_9km + veg_cover,
                                   occu.odhe)

Shrub_Dense_Add_Turbine_Dense <- occu( ~ max_trig_dist + livestock.count 
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
  "vegetation cover + shrub density" = Habitat,
  "Turbine density + vegetation cover + shrub density" = Habitat_Add_Turbine_Dense,
  "Turbine density x vegetation cover + shrub density" = Veg_Cov_X_Turbine_Dense,
  "Turbine density x biotic community + shrub density" = Bio_Com_X_Turbine_Dense,
  "Null" = Null,
  "Turbine density" = Turbine_Dense_Occu,
  "vegetation cover" = Veg_Cov_Occu,
  "Shrub density" = Shrub_Dense_Occu,
  "Biotic community + shrub density" = Bio_Com_Shrub,
  "Turbine density + vegetation cover" = Veg_Cov_Add_Turbine_Dense,
  "Turbine density + shrub density" = Shrub_Dense_Add_Turbine_Dense
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine density" = "Null",
  "Turbine density + vegetation cover + shrub density" = 
    "vegetation cover + shrub density",
  "Turbine density + shrub density" = "Shrub density",
  "Turbine density + vegetation cover" = "vegetation cover",
  "Turbine density x vegetation cover + shrub density" = 
    "vegetation cover + shrub density",
  "Turbine density x biotic community + shrub density" = 
    "Biotic community + shrub density"
)

# Compile final table and export to .csv
final_table <- occu_model_selection(
  model_list = model_list,
  wind_pairs = wind_pairs,
  file_name = paste0(homewd, "outputs/mule_deer_turbine_density.csv")
)

#### Access Road Distance Models ####

Bio_Com_X_Turbine_Rd_Dist <- occu( ~ max_trig_dist + livestock.count 
                                   ~ turbine_rd_dist * as.factor(biotic_com_2) +
                                     shrub_yucca_density, occu.odhe)

Veg_Cov_X_Turbine_Rd_Dist <- occu( ~ max_trig_dist + livestock.count 
                                   ~ veg_cover * turbine_rd_dist +
                                     shrub_yucca_density, occu.odhe)

Habitat_Add_Turbine_Rd_Dist <- occu( ~ max_trig_dist + livestock.count 
                                     ~ shrub_yucca_density + 
                                       veg_cover + turbine_rd_dist, 
                                     occu.odhe)

Turbine_Rd_Dist_Occu <- occu( ~ max_trig_dist + livestock.count ~ 
                              ~ turbine_rd_dist, occu.odhe,
                                starts = c(5,1,-3,0,-2))

Veg_Cov_Add_Turbine_Rd_Dist <- occu( ~ max_trig_dist + livestock.count 
                                     ~ turbine_rd_dist + veg_cover,
                                     occu.odhe)

Shrub_Dense_Add_Turbine_Rd_Dist <- occu( ~ max_trig_dist + livestock.count 
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
  "vegetation cover + shrub density" = Habitat,
  "Turbine road distance + vegetation cover + shrub density" = Habitat_Add_Turbine_Rd_Dist,
  "Turbine road distance x vegetation cover + shrub density" = Veg_Cov_X_Turbine_Rd_Dist,
  "Turbine road distance x biotic community + shrub density" = Bio_Com_X_Turbine_Rd_Dist,
  "Null" = Null,
  "Turbine road distance" = Turbine_Rd_Dist_Occu,
  "vegetation cover" = Veg_Cov_Occu,
  "Shrub density" = Shrub_Dense_Occu,
  "Biotic community + shrub density" = Bio_Com_Shrub,
  "Turbine road distance + vegetation cover" = Veg_Cov_Add_Turbine_Rd_Dist,
  "Turbine road distance + shrub density" = Shrub_Dense_Add_Turbine_Rd_Dist
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine road distance" = "Null",
  "Turbine road distance + vegetation cover + shrub density" = 
    "vegetation cover + shrub density",
  "Turbine road distance + shrub density" = "Shrub density",
  "Turbine road distance + vegetation cover" = "vegetation cover",
  "Turbine road distance x vegetation cover + shrub density" = 
    "vegetation cover + shrub density",
  "Turbine road distance x biotic community + shrub density" = 
    "Biotic community + shrub density"
)

# Compile final table and export to .csv
final_table <- occu_model_selection(
  model_list = model_list,
  wind_pairs = wind_pairs,
  file_name = paste0(homewd, "outputs/mule_deer_turbine_rd_dist.csv")
)

#### Access Road Density Models ####

Bio_Com_X_Turbine_Rd_Dense <- occu( ~ max_trig_dist + livestock.count 
                                    ~ turbine_rd_density_1_9km * 
                                      as.factor(biotic_com_2) +
                                      shrub_yucca_density, occu.odhe)

Veg_Cov_X_Turbine_Rd_Dense <- occu( ~ max_trig_dist + livestock.count ~ 
                                      ~ veg_cover *
                                      turbine_rd_density_1_9km +
                                      shrub_yucca_density, occu.odhe)

Habitat_Add_Turbine_Rd_Dense <- occu( ~ max_trig_dist + livestock.count 
                                      ~ shrub_yucca_density + 
                                        veg_cover +
                                        turbine_rd_density_1_9km, occu.odhe)

Turbine_Rd_Dense_Occu <- occu( ~ max_trig_dist + livestock.count 
                               ~ turbine_rd_density_1_9km, occu.odhe)

Veg_Cov_Add_Turbine_Rd_Dense <- occu( ~ max_trig_dist + livestock.count 
                                      ~ turbine_rd_density_1_9km + 
                                        veg_cover, occu.odhe)

Shrub_Rd_Dense_Add_Turbine_Rd_Dense <- occu( ~ max_trig_dist + livestock.count 
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
  "vegetation cover + shrub density" = Habitat,
  "Turbine road density + vegetation cover + shrub density" = Habitat_Add_Turbine_Rd_Dense,
  "Turbine road density x vegetation cover + shrub density" = Veg_Cov_X_Turbine_Rd_Dense,
  "Turbine road density x biotic community + shrub density" = Bio_Com_X_Turbine_Rd_Dense,
  "Null" = Null,
  "Turbine road density" = Turbine_Rd_Dense_Occu,
  "vegetation cover" = Veg_Cov_Occu,
  "Shrub density" = Shrub_Dense_Occu,
  "Biotic community + shrub density" = Bio_Com_Shrub,
  "Turbine road density + vegetation cover" = Veg_Cov_Add_Turbine_Rd_Dense,
  "Turbine road density + shrub density" = Shrub_Rd_Dense_Add_Turbine_Rd_Dense
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine road density" = "Null",
  "Turbine road density + vegetation cover + shrub density" = 
    "vegetation cover + shrub density",
  "Turbine road density + shrub density" = "Shrub density",
  "Turbine road density + vegetation cover" = "vegetation cover",
  "Turbine road density x vegetation cover + shrub density" = 
    "vegetation cover + shrub density",
  "Turbine road density x biotic community + shrub density" = 
    "Biotic community + shrub density"
)

# Compile final table and export to .csv
final_table <- occu_model_selection(
  model_list = model_list,
  wind_pairs = wind_pairs,
  file_name = paste0(homewd, "outputs/mule_deer_turbine_rd_density.csv")
)

#### Wind Variable Additive Mods ####

Additive1 <- occu( ~ max_trig_dist + livestock.count
                   ~ as.factor(biotic_com_2) * X150cm_turbine_vis +
                     shrub_yucca_density + turbine_dist, occu.odhe,
                     starts = c(0,5,-1,0,0,5,-3,-1,-2))  
               # large SE for biotic community parameter and interaction term

Additive2 <- occu( ~ max_trig_dist + livestock.count 
                   ~ X150cm_turbine_vis + shrub_yucca_density + turbine_dist + 
                     veg_cover, occu.odhe)  

Additive3 <- occu( ~ max_trig_dist + livestock.count 
                   ~ X150cm_turbine_vis + shrub_yucca_density + turbine_dist, 
                   occu.odhe)  

Additive4 <- occu( ~ max_trig_dist + livestock.count 
                   ~ X150cm_turbine_vis + veg_cover + turbine_dist, 
                   occu.odhe)  

Additive5 <- occu( ~ max_trig_dist + livestock.count 
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
  "Turbine visibility + shrub density + turbine distance + vegetation cover" =
    Additive2, 
  "Turbine visibility + shrub density + turbine distance" = Additive3, 
  "Turbine visibility + vegetation cover + turbine distance" = Additive4,
  "Turbine visibility + turbine distance" = Additive5, 
  "Null Model" = Null,
  "vegetation cover" = Veg_Cov_Occu,
  "Shrub density" = Shrub_Dense_Occu,
  "vegetation cover + shrub density" = Habitat,
  "Biotic community + shrub density" = Bio_Com_Shrub
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine visibility x biotic community + shrub density + turbine distance" =
    "Biotic community + shrub density",
  "Turbine visibility + shrub density + turbine distance + vegetation cover" =
    "vegetation cover + shrub density", 
  "Turbine visibility + shrub density + turbine distance" = "Shrub density", 
  "Turbine visibility + vegetation cover + turbine distance" = 
    "vegetation cover",
  "Turbine visibility + turbine distance" = "Null Model"
)

# Compile final table and export to .csv
final_table <- occu_model_selection(
  model_list = model_list,
  wind_pairs = wind_pairs,
  file_name = paste0(homewd, "outputs/mule_deer_additive_mods.csv")
)

########################################################
######################### END ##########################
########################################################