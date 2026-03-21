###########################################################
########### Code for: base models for kit foxes ###########
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
library(AICcmodavg)

###########################################################
# SETUP CODE FOR KIT FOX OCCUPANCY MODELS #
###########################################################

# Read in edited .csv
detHist <- read.csv(file = "kit fox detection hist.csv", row.names = 1)
# Read in site.covs.scaled
site.covs.scaled <- readRDS("site_covs_scaled.RData")
# Read in obsCovs.scaled
obsCovs.scaled <- readRDS("obsCovs_scaled.RData")

# Change from integer to numeric
detHist %>%
  mutate(across(1:51, as.numeric))

# Create detection history 
occu.vuma <- unmarkedFrameOccu(y=detHist, 
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

Habitat <- occu( ~ as.factor(cam_moved) + NDVI_1_9km
                 ~ jackrabbit_count_avg + coy_count_avg, occu.vuma)

Null <- occu( ~ as.factor(cam_moved) + NDVI_1_9km ~ 1, occu.vuma)

Jack_Occu <- occu( ~ as.factor(cam_moved) + NDVI_1_9km  
                   ~ jackrabbit_count_avg, occu.vuma)

Coy_Occu <- occu( ~ as.factor(cam_moved) + NDVI_1_9km  
                  ~ coy_count_avg, occu.vuma)

Bio_Com_Occu <- occu( ~ as.factor(cam_moved) + NDVI_1_9km 
                      ~ as.factor(biotic_com_2), occu.vuma,
                        starts = c(0, -5, -5, -2, -2))
#large SE for biotic comm

Bio_Com_Habitat <- occu( ~ as.factor(cam_moved) + NDVI_1_9km 
                         ~ as.factor(biotic_com_2) + coy_count_avg + 
                           jackrabbit_count_avg, occu.vuma,
                           starts = c(-1, -5, -5, 5, -5, -1, -2))
#large SE for biotic comm

Canopy_Habitat <- occu( ~ as.factor(cam_moved) + NDVI_1_9km 
                        ~ canopy_cov + coy_count_avg + 
                          jackrabbit_count_avg, occu.vuma,
                          starts = c(-5, -5, -5, 5, -5, -2, -2))
#large SE for canopy cov and intercept

## The table with the comparison of the relative weight of evidence between
  #occupancy models separately ranked for the effect of each wind energy 
  #variable on the probability of habitat selection (ψ) of kit foxes is found 
  #in Table S3.17.

## Check correlations between wind variables and habitat variables 

# Read in site-level covariates
site.covs <- read.csv("site_covs.csv", nrows = 102, header = TRUE)

# site-level variable correlations
wind.hab.cor <- site.covs %>% 
  select(jackrabbit_count_avg, biotic_com_2, coy_count_avg, canopy_cov,
         turbine_interior, X50cm_turbine_vis, turbine_dist,
         turbine_density_1_9km, turbine_rd_dist, turbine_rd_density_1_9km)

cors <- cor(wind.hab.cor, method='spearman')  
# none correlated above |0.7|

#### Turbine Interior Models ####

Turbine_Int_Occu <- occu( ~ as.factor(cam_moved) + NDVI_1_9km 
                          ~ as.factor(turbine_interior), occu.vuma)

Habitat_Add_Turbine_Int <- occu( ~ as.factor(cam_moved) + NDVI_1_9km
                                 ~ as.factor(turbine_interior) +
                                   jackrabbit_count_avg + coy_count_avg, 
                                 occu.vuma,
                                 starts = c(0, 0, 0.1, 0.1, 0, 1, 0.1))

Jack_Add_Turbine_Int <- occu( ~ as.factor(cam_moved) + NDVI_1_9km
                              ~ as.factor(turbine_interior) + 
                                jackrabbit_count_avg, occu.vuma)

Coy_Add_Turbine_Int <- occu( ~ as.factor(cam_moved) + NDVI_1_9km
                             ~ as.factor(turbine_interior) + coy_count_avg, 
                             occu.vuma)

Bio_Com_X_Turbine_Int <- occu( ~ as.factor(cam_moved) 
                               ~ as.factor(turbine_interior) *
                                 as.factor(biotic_com_2) + 
                                 jackrabbit_count_avg + coy_count_avg, occu.vuma,
                                 starts = c(-1, 0, -5, 2, -1, 2, -2, -1))
#large SE for biotic comm and interaction

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
slope_name       <- "psi(jackrabbit_count_avg)"
vertical_cover   <- "psi(coy_count_avg)"

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
  "Jackrabbit count + coyote count" = Habitat,
  "Null" = Null,
  "Biotic community" = Bio_Com_Occu,
  "Jackrabbit count" = Jack_Occu,
  "Coyote count" = Coy_Occu,
  "Turbine interior" = Turbine_Int_Occu,
  "Biotic community +jackrabbit count + coyote count" = Bio_Com_Habitat,
  "Turbine interior + jackrabbit count + coyote count" = Habitat_Add_Turbine_Int,
  "Turbine interior + jackrabbit count" = Jack_Add_Turbine_Int,
  "Turbine interior + coyote count" = Coy_Add_Turbine_Int,
  "Turbine interior x biotic community + jackrabbit count + coyote count" =
    Bio_Com_X_Turbine_Int
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine interior" = "Null",
  "Turbine interior + jackrabbit count + coyote count" = 
    "Jackrabbit count + coyote count",
  "Turbine interior + jackrabbit count" = "Jackrabbit count",
  "Turbine interior + coyote count" = "Coyote count",
  "Turbine interior x biotic community + jackrabbit count + coyote count" =
    "Biotic community +jackrabbit count + coyote count" 
)

# Compile final table and export to .csv
final_table <- occu_model_selection(
  model_list = model_list,
  wind_pairs = wind_pairs,
  file_name = paste0(homewd, "outputs/kit_fox_turbine_interior.csv")
)

#### Turbine Visibility Models ####

Turbine_Vis_Occu <- occu( ~ as.factor(cam_moved) + NDVI_1_9km 
                          ~ X50cm_turbine_vis, occu.vuma)

Habitat_Add_Turbine_Vis <- occu( ~ as.factor(cam_moved) + NDVI_1_9km
                                 ~ X50cm_turbine_vis + jackrabbit_count_avg + 
                                   coy_count_avg, occu.vuma,
                                   starts = c(0, 0, 0.1, 0.1, 0, 1, 0.1))

Jack_Add_Turbine_Vis <- occu( ~ as.factor(cam_moved) + NDVI_1_9km
                              ~ X50cm_turbine_vis + jackrabbit_count_avg, 
                                occu.vuma)
 
Coy_Add_Turbine_Vis <- occu( ~ as.factor(cam_moved) + NDVI_1_9km
                             ~ X50cm_turbine_vis + coy_count_avg, occu.vuma) 

Bio_Com_X_Turbine_Vis <- occu( ~ as.factor(cam_moved) 
                               ~ X50cm_turbine_vis * as.factor(biotic_com_2) + 
                                 jackrabbit_count_avg + coy_count_avg, occu.vuma,
                                 starts = c(-2, 0, -5, 2, -1, 5, -2, -1))
#large SE for biotic comm and interaction

Canopy_X_Turbine_Vis <- occu( ~ as.factor(cam_moved) 
                              ~ canopy_cov * X50cm_turbine_vis + 
                                jackrabbit_count_avg + coy_count_avg, occu.vuma,
                                starts = c(-5, -5, 5, 1, -2, 5, -2, -1))
#large SE for occupancy parameters 

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
  "Jackrabbit count + coyote count" = Habitat,
  "Null" = Null,
  "Biotic community + jackrabbit count + coyote count" = Bio_Com_Habitat,
  "Canopy cover + jackrabbit count + coyote count" = Canopy_Habitat,  
  "Jackrabbit count" = Jack_Occu,
  "Coyote count" = Coy_Occu,
  "Turbine visibility" = Turbine_Vis_Occu,
  "Turbine visibility + jackrabbit count + coyote count" = 
    Habitat_Add_Turbine_Vis,
  "Turbine visibility + jackrabbit count" = Jack_Add_Turbine_Vis,
  "Turbine visibility + coyote count" = Coy_Add_Turbine_Vis,
  "Turbine visibility x biotic community + jackrabbit count + coyote count" =
    Bio_Com_X_Turbine_Vis,
  "Turbine visibility x canopy cover + jackrabbit count + coyote count" =
    Canopy_X_Turbine_Vis
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine visibility" = "Null",
  "Turbine visibility + jackrabbit count + coyote count" = 
    "Jackrabbit count + coyote count",
  "Turbine visibility + jackrabbit count" = "Jackrabbit count",
  "Turbine visibility + coyote count" = "Coyote count",
  "Turbine visibility x biotic community + jackrabbit count + coyote count" =
    "Biotic community + jackrabbit count + coyote count",
  "Turbine visibility x canopy cover + jackrabbit count + coyote count" = 
    "Canopy cover + jackrabbit count + coyote count"
)

# Compile final table and export to .csv
final_table <- occu_model_selection(
  model_list = model_list,
  wind_pairs = wind_pairs,
  file_name = paste0(homewd, "outputs/kit_fox_turbine_vis.csv")
)

#### Turbine Distance Models ####

Turbine_Dist_Occu <- occu( ~ as.factor(cam_moved) + NDVI_1_9km 
                           ~ turbine_dist, occu.vuma)

Habitat_Add_Turbine_Dist <- occu( ~ as.factor(cam_moved) + NDVI_1_9km
                                  ~ turbine_dist + jackrabbit_count_avg +
                                    coy_count_avg, occu.vuma)

Jack_Add_Turbine_Dist <- occu( ~ as.factor(cam_moved) + NDVI_1_9km
                               ~ turbine_dist + jackrabbit_count_avg, 
                               occu.vuma)

Coy_Add_Turbine_Dist <- occu( ~ as.factor(cam_moved) + NDVI_1_9km
                              ~ turbine_dist + coy_count_avg, occu.vuma)

Bio_Com_X_Turbine_Dist <- occu( ~ as.factor(cam_moved) 
                                ~ turbine_dist * as.factor(biotic_com_2) +
                                  jackrabbit_count_avg + coy_count_avg, 
                                  occu.vuma,
                                  starts = c(-1, 0, -10, 2, -1, 2, -2, -1))
#large SE for biotic comm and interaction

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
  "Jackrabbit count + coyote count" = Habitat,
  "Null" = Null,
  "Biotic community + jackrabbit count + coyote count" = Bio_Com_Habitat,
  "Jackrabbit count" = Jack_Occu,
  "Coyote count" = Coy_Occu,
  "Turbine distance" = Turbine_Dist_Occu,
  "Turbine distance + jackrabbit count + coyote count" = 
    Habitat_Add_Turbine_Dist,
  "Turbine distance + jackrabbit count" = Jack_Add_Turbine_Dist,
  "Turbine distance + coyote count" = Coy_Add_Turbine_Dist,
  "Turbine distance x biotic community + jackrabbit count + coyote count" =
    Bio_Com_X_Turbine_Dist
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine distance" = "Null",
  "Turbine distance + jackrabbit count + coyote count" = 
    "Jackrabbit count + coyote count",
  "Turbine distance + jackrabbit count" = "Jackrabbit count",
  "Turbine distance + coyote count" = "Coyote count",
  "Turbine distance x biotic community + jackrabbit count + coyote count" =
    "Biotic community + jackrabbit count + coyote count"
)

# Compile final table and export to .csv
final_table <- occu_model_selection(
  model_list = model_list,
  wind_pairs = wind_pairs,
  file_name = paste0(homewd, "outputs/kit_fox_turbine_dist.csv")
)

#### Turbine Density Models ####

Turbine_Dense_Occu <- occu( ~ as.factor(cam_moved) + NDVI_1_9km 
                            ~ turbine_density_1_9km, occu.vuma)

Habitat_Add_Turbine_Dense <- occu( ~ as.factor(cam_moved) + NDVI_1_9km
                                   ~ turbine_density_1_9km +
                                     jackrabbit_count_avg + coy_count_avg, 
                                     occu.vuma)

Jack_Add_Turbine_Dense <- occu( ~ as.factor(cam_moved) + NDVI_1_9km
                                ~ turbine_density_1_9km + jackrabbit_count_avg, 
                                  occu.vuma)

Coy_Add_Turbine_Dense <- occu( ~ as.factor(cam_moved) + NDVI_1_9km
                               ~ turbine_density_1_9km + coy_count_avg, 
                                 occu.vuma)

Bio_Com_X_Turbine_Dense <- occu( ~ as.factor(cam_moved) 
                                 ~ turbine_density_1_9km *
                                   as.factor(biotic_com_2) + 
                                   jackrabbit_count_avg + coy_count_avg, 
                                   occu.vuma,
                                   starts = c(-1, 0, -10, 2, -1, -2, -2, -1))
#large SE for biotic comm and interaction

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
  "Jackrabbit count + coyote count" = Habitat,
  "Null" = Null,
  "Biotic community + jackrabbit count + coyote count" = Bio_Com_Habitat,
  "Jackrabbit count" = Jack_Occu,
  "Coyote count" = Coy_Occu,
  "Turbine density" = Turbine_Dense_Occu,
  "Turbine density + jackrabbit count + coyote count" = Habitat_Add_Turbine_Dense,
  "Turbine density + jackrabbit count" = Jack_Add_Turbine_Dense,
  "Turbine density + coyote count" = Coy_Add_Turbine_Dense,
  "Turbine density x biotic community + jackrabbit count + coyote count" =
    Bio_Com_X_Turbine_Dense
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine density" = "Null",
  "Turbine density + jackrabbit count + coyote count" = 
    "Jackrabbit count + coyote count",
  "Turbine density + jackrabbit count" = "Jackrabbit count",
  "Turbine density + coyote count" = "Coyote count",
  "Turbine density x biotic community + jackrabbit count + coyote count" =
    "Biotic community + jackrabbit count + coyote count"
)


# Compile final table and export to .csv
final_table <- occu_model_selection(
  model_list = model_list,
  wind_pairs = wind_pairs,
  file_name = paste0(homewd, "outputs/kit_fox_turbine_density.csv")
)

#### Access Road Distance Models ####

Turbine_Rd_Dist_Occu <- occu( ~ as.factor(cam_moved) + NDVI_1_9km 
                              ~ turbine_rd_dist, occu.vuma)

Habitat_Add_Turbine_Rd_Dist <- occu( ~ as.factor(cam_moved) + NDVI_1_9km
                                     ~ turbine_rd_dist + jackrabbit_count_avg +
                                       coy_count_avg, occu.vuma,
                                       starts = c(0, -1, 5, -5, -5, -1, -2))

Jack_Add_Turbine_Rd_Dist <- occu( ~ as.factor(cam_moved) + NDVI_1_9km
                                  ~ turbine_rd_dist + jackrabbit_count_avg, 
                                  occu.vuma)

Coy_Add_Turbine_Rd_Dist <- occu( ~ as.factor(cam_moved) + NDVI_1_9km
                                 ~ turbine_rd_dist + coy_count_avg, occu.vuma)

Bio_Com_X_Turbine_Rd_Dist <- occu( ~ as.factor(cam_moved)
                                   ~ turbine_rd_dist * as.factor(biotic_com_2) +
                                     jackrabbit_count_avg + coy_count_avg, 
                                     occu.vuma,
                                     starts = c(-1, 0, -10, 2, -1, 1, -2, -1))
#large SE for biotic comm and interaction

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
  "Jackrabbit count + coyote count" = Habitat,
  "Null" = Null,
  "Biotic community + jackrabbit count + coyote count" = Bio_Com_Habitat,
  "Jackrabbit count" = Jack_Occu,
  "Coyote count" = Coy_Occu,
  "Turbine road distance" = Turbine_Rd_Dist_Occu,
  "Turbine road distance + jackrabbit count + coyote count" = 
    Habitat_Add_Turbine_Rd_Dist,
  "Turbine road distance + jackrabbit count" = Jack_Add_Turbine_Rd_Dist,
  "Turbine road distance + coyote count" = Coy_Add_Turbine_Rd_Dist,
  "Turbine road distance x biotic community + jackrabbit count + coyote count" =
    Bio_Com_X_Turbine_Rd_Dist
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine road distance" = "Null",
  "Turbine road distance + jackrabbit count + coyote count" = 
    "Jackrabbit count + coyote count",
  "Turbine road distance + jackrabbit count" = "Jackrabbit count",
  "Turbine road distance + coyote count" = "Coyote count",
  "Turbine road distance x biotic community + jackrabbit count + coyote count" =
    "Biotic community + jackrabbit count + coyote count"
)


# Compile final table and export to .csv
final_table <- occu_model_selection(
  model_list = model_list,
  wind_pairs = wind_pairs,
  file_name = paste0(homewd, "outputs/kit_fox_turbine_rd_dist.csv")
)

#### Access Road Density Models ####

Turbine_Rd_Dense_Occu <- occu( ~ as.factor(cam_moved) + NDVI_1_9km 
                               ~ turbine_rd_density_1_9km, occu.vuma)

Habitat_Add_Turbine_Rd_Dense <- occu( ~ as.factor(cam_moved) + NDVI_1_9km
                                      ~ turbine_rd_density_1_9km +
                                        jackrabbit_count_avg + coy_count_avg, 
                                        occu.vuma)

Jack_Add_Turbine_Rd_Dense <- occu( ~ as.factor(cam_moved) + NDVI_1_9km
                                   ~ turbine_rd_density_1_9km +
                                     jackrabbit_count_avg, occu.vuma)

Coy_Add_Turbine_Rd_Dense <- occu( ~ as.factor(cam_moved) + NDVI_1_9km
                                  ~ turbine_rd_density_1_9km + coy_count_avg, 
                                    occu.vuma)

Bio_Com_X_Turbine_Rd_Dense <- occu( ~ as.factor(cam_moved) 
                                    ~ turbine_rd_density_1_9km *
                                      as.factor(biotic_com_2) +
                                      jackrabbit_count_avg + coy_count_avg, 
                                      occu.vuma,
                                      starts = c(-1, 0, -10, 2, -2, -5, -2, -1))
#large SE for biotic comm and interaction

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

# AICc Table with all possible turbine rd density models 

model_list <- list(
  "Jackrabbit count + coyote count" = Habitat,
  "Null" = Null,
  "Biotic community + jackrabbit count + coyote count" = Bio_Com_Habitat,
  "Jackrabbit count" = Jack_Occu,
  "Coyote count" = Coy_Occu,
  "Turbine road density" = Turbine_Rd_Dense_Occu,
  "Turbine road density + jackrabbit count + coyote count" = 
    Habitat_Add_Turbine_Rd_Dense,
  "Turbine road density + jackrabbit count" = Jack_Add_Turbine_Rd_Dense,
  "Turbine road density + coyote count" = Coy_Add_Turbine_Rd_Dense,
  "Turbine road density x biotic community + jackrabbit count + coyote count" =
    Bio_Com_X_Turbine_Rd_Dense
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine road density" = "Null",
  "Turbine road density + jackrabbit count + coyote count" = 
    "Jackrabbit count + coyote count",
  "Turbine road density + jackrabbit count" = "Jackrabbit count",
  "Turbine road density + coyote count" = "Coyote count",
  "Turbine road density x biotic community + jackrabbit count + coyote count" =
    "Biotic community + jackrabbit count + coyote count"
)


# Compile final table and export to .csv
final_table <- occu_model_selection(
  model_list = model_list,
  wind_pairs = wind_pairs,
  file_name = paste0(homewd, "outputs/kit_fox_turbine_rd_density.csv")
)

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
  ER <- if(all(c(mod_with_wind, mod_without) %in% 
               names(w))) w[mod_with_wind] / w[mod_without] else NA
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

#write.csv(final_table, "kit_fox_turbine_interior.csv", row.names = FALSE)

#write.csv(final_table, "kit_fox_turbine_vis.csv", row.names = FALSE)

#write.csv(final_table, "kit_fox_turbine_dist.csv", row.names = FALSE)

#write.csv(final_table, "kit_fox_turbine_density.csv", row.names = FALSE)

#write.csv(final_table, "kit_fox_turbine_rd_dist.csv", row.names = FALSE)

#write.csv(final_table, "kit_fox_turbine_rd_density.csv", row.names = FALSE)


########################################################
######################### END ##########################
########################################################