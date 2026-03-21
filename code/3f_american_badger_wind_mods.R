###########################################################
######## Code for: base models for American badgers #######
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
# SETUP CODE FOR AMERICAN BADGER OCCUPANCY MODELS #
###########################################################

# Read in edited .csv
detHist <- read.csv(file = "badger detection hist.csv", row.names = 1)
# Read in site.covs.scaled
site.covs.scaled <- readRDS("site_covs_scaled.RData")
# Read in obsCovs.scaled
obsCovs.scaled <- readRDS("obsCovs_scaled.RData")

# Change from integer to numeric

detHist <- detHist %>%
  mutate(across(1:51, as.numeric))

# Create detection history 
occu.tata <- unmarkedFrameOccu(y=detHist, 
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

Habitat <- occu( ~ detection_angle + as.factor(precip.cat)
                 ~ vertical_cover + slope, occu.tata)

Null <- occu( ~ detection_angle + as.factor(precip.cat) ~ 1, occu.tata)

Slope_Occu <- occu( ~ detection_angle + as.factor(precip.cat)
                    ~ slope, occu.tata)

Canopy_Habitat <- occu( ~ detection_angle + as.factor(precip.cat)
                        ~ slope + vertical_cover + canopy_cov, occu.tata)

Cov_Occu <- occu( ~ detection_angle + as.factor(precip.cat)
                  ~ vertical_cover, occu.tata)

Bio_Com_Cov_Slope <- occu( ~ detection_angle + as.factor(precip.cat)
                           ~  vertical_cover + as.factor(biotic_com_2) + slope,
                              occu.tata)

## The table with the comparison of the relative weight of evidence between
  #occupancy models separately ranked for the effect of each wind energy 
  #variable on the probability of habitat selection (ψ) of American badgers 
  #is found in Table S3.18.

## Check correlations between wind variables and habitat variables 

# Read in site-level covariates
site.covs <- read.csv("site_covs.csv", nrows = 102, header = TRUE)

# site-level variable correlations
wind.hab.cor <- site.covs %>% 
  select( biotic_com_2, vertical_cover, canopy_cov, slope,
          turbine_interior, X50cm_turbine_vis, turbine_dist,
          turbine_density_1_6km, turbine_rd_dist, turbine_rd_density_1_6km)

cors <- cor(wind.hab.cor, method='spearman')  
# none correlated above |0.7|

#### Turbine Interior Models ####

Turbine_Int_Occu <- occu( ~ detection_angle + as.factor(precip.cat) 
                          ~ as.factor(turbine_interior), occu.tata)

Habitat_Add_Turbine_Int <- occu( ~ detection_angle + as.factor(precip.cat)
                                 ~  as.factor(turbine_interior) +
                                   vertical_cover + slope, occu.tata)

Cov_Add_Turbine_Int <- occu( ~ detection_angle + as.factor(precip.cat)
                             ~  as.factor(turbine_interior) +
                               vertical_cover, occu.tata)

Slope_Add_Turbine_Int <- occu( ~ detection_angle + as.factor(precip.cat)
                               ~  as.factor(turbine_interior) + slope, 
                               occu.tata)

Bio_Com_X_Turbine_Int_Slope <- occu( ~ detection_angle + as.factor(precip.cat)
                                     ~ as.factor(turbine_interior) *
                                       as.factor(biotic_com_2) + slope, occu.tata)

Bio_Com_X_Turbine_Int_Cov <- occu( ~ detection_angle + as.factor(precip.cat)
                                   ~ as.factor(turbine_interior) *
                                     as.factor(biotic_com_2) + vertical_cover, 
                                   occu.tata)

Bio_Com_X_Turbine_Int_Habitat <- occu( ~ detection_angle + as.factor(precip.cat)
                                       ~ as.factor(turbine_interior) *
                                         as.factor(biotic_com_2) + 
                                         vertical_cover + slope, occu.tata)

# 85% CI for bio community interaction 
# Extract coefficients and VCOV matrix
coef_est <- coef(Bio_Com_X_Turbine_Int_Habitat, type = "state")
vcov_mat <- vcov(Bio_Com_X_Turbine_Int_Habitat, type = "state")

# Exact coefficient names from your model
intercept_name   <- "psi(Int)"
turbine_name     <- "psi(as.factor(turbine_interior)1)"
habitat_name     <- "psi(as.factor(biotic_com_2)2)"
interaction_name <- 
  "psi(as.factor(turbine_interior)1:as.factor(biotic_com_2)2)"
slope_name       <- "psi(slope)"
vertical_cover   <- "psi(vertical_cover)"

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
  "Vertical cover + slope" = Habitat,
  "Null" = Null,
  "Slope" = Slope_Occu,
  "Vertical cover" = Cov_Occu,
  "Biotic community + slope + vertical cover" = Bio_Com_Cov_Slope,
  "Turbine interior" = Turbine_Int_Occu,
  "Turbine interior + vertical cover + slope" = Habitat_Add_Turbine_Int,
  "Turbine interior + vertical cover" = Cov_Add_Turbine_Int,
  "Turbine interior + slope" = Slope_Add_Turbine_Int,
  "Turbine interior x biotic community + slope + vertical cover" =
    Bio_Com_X_Turbine_Int_Habitat
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine interior" = "Null",
  "Turbine interior + vertical cover + slope" = "Vertical cover + slope",
  "Turbine interior + vertical cover" = "Vertical cover",
  "Turbine interior + slope" = "Slope",
  "Turbine interior x biotic community + slope + vertical cover" =
    "Biotic community + slope + vertical cover"
)

# Compile final table and export to .csv
final_table <- occu_model_selection(
  model_list = model_list,
  wind_pairs = wind_pairs,
  file_name = paste0(homewd, "outputs/badger_turbine_interior.csv")
)

#### Turbine Visibility Models ####

Turbine_Vis_Occu <- occu( ~ detection_angle + as.factor(precip.cat) 
                          ~ X50cm_turbine_vis, occu.tata)

Habitat_Add_Turbine_Vis <- occu( ~ detection_angle + as.factor(precip.cat)
                                 ~ X50cm_turbine_vis + vertical_cover + slope,
                                 occu.tata)

Cov_Add_Turbine_Vis <- occu( ~ detection_angle + as.factor(precip.cat)
                             ~ X50cm_turbine_vis + vertical_cover, occu.tata)

Slope_Add_Turbine_Vis <- occu( ~ detection_angle + as.factor(precip.cat)
                               ~ X50cm_turbine_vis + slope, occu.tata)

Bio_Com_X_Turbine_Vis <- occu( ~ detection_angle + as.factor(precip.cat)
                               ~ X50cm_turbine_vis *
                                 as.factor(biotic_com_2) + vertical_cover +
                                 slope, occu.tata, 
                                 starts = c(0,0,-1,0,-1,-3,-3,0,0))
#large SE for biotic comm and interaction 

Canopy_Cov_X_Turbine_Vis <- occu( ~ detection_angle + as.factor(precip.cat)
                                  ~ X50cm_turbine_vis * canopy_cov + 
                                    vertical_cover + slope, occu.tata)

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
  "Vertical cover + slope" = Habitat,
  "Null" = Null,
  "Biotic community + slope + vertical cover" = Bio_Com_Cov_Slope,
  "Slope" = Slope_Occu,
  "Vertical cover" = Cov_Occu,
  "Turbine visibility" = Turbine_Vis_Occu,
  "Canopy cover + slope + vertical cover" = Canopy_Habitat,
  "Turbine visibility + vertical cover + slope" = Habitat_Add_Turbine_Vis,
  "Turbine visibility + vertical cover" = Cov_Add_Turbine_Vis,
  "Turbine visibility + slope" = Slope_Add_Turbine_Vis,
  "Turbine visibility x biotic community + slope + vertical cover" = 
    Bio_Com_X_Turbine_Vis,
  "Turbine visibility x canopy cover + slope + vertical cover" = 
    Canopy_Cov_X_Turbine_Vis
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine visibility" = "Null",
  "Turbine visibility + vertical cover + slope" = "Vertical cover + slope",
  "Turbine visibility + vertical cover" = "Vertical cover",
  "Turbine visibility + slope" = "Slope",
  "Turbine visibility x canopy cover + slope + vertical cover" = 
    "Canopy cover + slope + vertical cover",
  "Turbine visibility x biotic community + slope + vertical cover" =
    "Biotic community + slope + vertical cover"
)

# Compile final table and export to .csv
final_table <- occu_model_selection(
  model_list = model_list,
  wind_pairs = wind_pairs,
  file_name = paste0(homewd, "outputs/badger_turbine_vis.csv")
)

#### Turbine Distance Models ####

Turbine_Dist_Occu <- occu( ~ detection_angle + as.factor(precip.cat) 
                           ~ turbine_dist, occu.tata)

Habitat_Add_Turbine_Dist <- occu( ~ detection_angle + as.factor(precip.cat)
                                  ~ turbine_dist + vertical_cover + slope,
                                  occu.tata)

Cov_Add_Turbine_Dist <- occu( ~ detection_angle + as.factor(precip.cat)
                              ~ turbine_dist + vertical_cover, occu.tata)

Slope_Add_Turbine_Dist <- occu( ~ detection_angle + as.factor(precip.cat)
                                ~ turbine_dist + slope, occu.tata)

Bio_Com_X_Turbine_Dist <- occu( ~ detection_angle + as.factor(precip.cat)
                                ~ turbine_dist * as.factor(biotic_com_2) +
                                  vertical_cover + slope, occu.tata)

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
  "Vertical cover + slope" = Habitat,
  "Null" = Null,
  "Biotic community + slope + vertical cover" = Bio_Com_Cov_Slope,
  "Slope" = Slope_Occu,
  "Vertical cover" = Cov_Occu,
  "Turbine distance" = Turbine_Dist_Occu,
  "Turbine distance + vertical cover + slope" = Habitat_Add_Turbine_Dist,
  "Turbine distance + vertical cover" = Cov_Add_Turbine_Dist,
  "Turbine distance + slope" = Slope_Add_Turbine_Dist,
  "Turbine distance x biotic community + slope + vertical cover" =
    Bio_Com_X_Turbine_Dist
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine distance" = "Null",
  "Turbine distance + vertical cover + slope" = "Vertical cover + slope",
  "Turbine distance + vertical cover" = "Vertical cover",
  "Turbine distance + slope" = "Slope",
  "Turbine distance x biotic community + slope + vertical cover" =
    "Biotic community + slope + vertical cover"
)

# Compile final table and export to .csv
final_table <- occu_model_selection(
  model_list = model_list,
  wind_pairs = wind_pairs,
  file_name = paste0(homewd, "outputs/badger_turbine_dist.csv")
)

#### Turbine Density Models ####

Turbine_Dense_Occu <- occu( ~ detection_angle + as.factor(precip.cat) 
                            ~ turbine_density_1_6km, occu.tata)

Habitat_Add_Turbine_Dense <- occu( ~ detection_angle + as.factor(precip.cat)
                                   ~  turbine_density_1_6km +
                                      vertical_cover + slope, occu.tata)

Cov_Add_Turbine_Dense <- occu( ~ detection_angle + as.factor(precip.cat)
                               ~  turbine_density_1_6km +
                                  vertical_cover, occu.tata)

Slope_Add_Turbine_Dense <- occu( ~ detection_angle + as.factor(precip.cat)
                                 ~  turbine_density_1_6km + slope, occu.tata)

Bio_Com_X_Turbine_Dense <- occu( ~ detection_angle + as.factor(precip.cat)
                                 ~ turbine_density_1_6km *
                                   as.factor(biotic_com_2) + vertical_cover +
                                   slope, occu.tata)

# 85% CI for bio community interaction 
# Extract coefficients and VCOV matrix
coef_est <- coef(Bio_Com_X_Turbine_Dense, type = "state")
vcov_mat <- vcov(Bio_Com_X_Turbine_Dense, type = "state")

beta_grassland <- coef_est["psi(turbine_density_1_6km)"]
se_grassland <- sqrt(vcov_mat["psi(turbine_density_1_6km)", 
                              "psi(turbine_density_1_6km)"])

beta_woodland <- coef_est["psi(turbine_density_1_6km)"] + 
  coef_est["psi(turbine_density_1_6km:as.factor(biotic_com_2)2)"]

var_woodland <- 
  vcov_mat["psi(turbine_density_1_6km)", "psi(turbine_density_1_6km)"] +
  vcov_mat["psi(turbine_density_1_6km:as.factor(biotic_com_2)2)",
           "psi(turbine_density_1_6km:as.factor(biotic_com_2)2)"] +
  2 * vcov_mat["psi(turbine_density_1_6km)", 
               "psi(turbine_density_1_6km:as.factor(biotic_com_2)2)"]
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
  "Vertical cover + slope" = Habitat,
  "Null" = Null,
  "Biotic community + slope + vertical cover" = Bio_Com_Cov_Slope,
  "Slope" = Slope_Occu,
  "Vertical cover" = Cov_Occu,
  "Turbine density" = Turbine_Dense_Occu,
  "Turbine density + vertical cover + slope" = Habitat_Add_Turbine_Dense,
  "Turbine density + vertical cover" = Cov_Add_Turbine_Dense,
  "Turbine density + slope" = Slope_Add_Turbine_Dense,
  "Turbine density x biotic community + slope + vertical cover" =
    Bio_Com_X_Turbine_Dense
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine density" = "Null",
  "Turbine density + vertical cover + slope" = "Vertical cover + slope",
  "Turbine density + vertical cover" = "Vertical cover",
  "Turbine density + slope" = "Slope",
  "Turbine density x biotic community + slope + vertical cover" =
    "Biotic community + slope + vertical cover"
)

# Compile final table and export to .csv
final_table <- occu_model_selection(
  model_list = model_list,
  wind_pairs = wind_pairs,
  file_name = paste0(homewd, "outputs/badger_turbine_density.csv")
)

#### Access Road Distance Models ####

Turbine_Rd_Dist_Occu <- occu( ~ detection_angle + as.factor(precip.cat) 
                              ~ turbine_rd_dist, occu.tata)

Habitat_Add_Turbine_Rd_Dist <- occu( ~ detection_angle + as.factor(precip.cat)
                                     ~ turbine_rd_dist + vertical_cover + slope,
                                     occu.tata)

Cov_Add_Turbine_Rd_Dist <- occu( ~ detection_angle + as.factor(precip.cat)
                                 ~  turbine_rd_dist + vertical_cover, occu.tata)

Slope_Add_Turbine_Rd_Dist <- occu( ~ detection_angle + as.factor(precip.cat)
                                   ~  turbine_rd_dist + slope, occu.tata)

Bio_Com_X_Turbine_Rd_Dist <- occu( ~ detection_angle + as.factor(precip.cat)
                                   ~ turbine_rd_dist * as.factor(biotic_com_2) +
                                     vertical_cover + slope, occu.tata)

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
  "Vertical cover + slope" = Habitat,
  "Null" = Null,
  "Biotic community + slope + vertical cover" = Bio_Com_Cov_Slope,
  "Slope" = Slope_Occu,
  "Vertical cover" = Cov_Occu,
  "Turbine road distance" = Turbine_Rd_Dist_Occu,
  "Turbine road distance + vertical cover + slope" = Habitat_Add_Turbine_Rd_Dist,
  "Turbine road distance + vertical cover" = Cov_Add_Turbine_Rd_Dist,
  "Turbine road distance + slope" = Slope_Add_Turbine_Rd_Dist,
  "Turbine road distance x biotic community + slope + vertical cover" =
    Bio_Com_X_Turbine_Rd_Dist
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine road distance" = "Null",
  "Turbine road distance + vertical cover + slope" = "Vertical cover + slope",
  "Turbine road distance + vertical cover" = "Vertical cover",
  "Turbine road distance + slope" = "Slope",
  "Turbine road distance x biotic community + slope + vertical cover" =
    "Biotic community + slope + vertical cover"
)

# Compile final table and export to .csv
final_table <- occu_model_selection(
  model_list = model_list,
  wind_pairs = wind_pairs,
  file_name = paste0(homewd, "outputs/badger_turbine_rd_dist.csv")
)

#### Access Road Density Models ####

Turbine_Rd_Dense_Occu <- occu( ~ detection_angle + as.factor(precip.cat) 
                               ~ turbine_rd_density_1_6km, occu.tata)

Habitat_Add_Turbine_Rd_Dense <- occu( ~ detection_angle + as.factor(precip.cat)
                                      ~  turbine_rd_density_1_6km +
                                         vertical_cover + slope, occu.tata)

Cov_Add_Turbine_Rd_Dense <- occu( ~ detection_angle + as.factor(precip.cat)
                                  ~ turbine_rd_density_1_6km + vertical_cover,
                                    occu.tata)

Slope_Add_Turbine_Rd_Dense <- occu( ~ detection_angle + as.factor(precip.cat)
                                    ~ turbine_rd_density_1_6km + slope, 
                                      occu.tata)

Bio_Com_X_Turbine_Rd_Dense <- occu( ~ detection_angle + as.factor(precip.cat)
                                    ~ turbine_rd_density_1_6km *
                                      as.factor(biotic_com_2) + vertical_cover +
                                      slope, occu.tata)

# 85% CI for bio community interaction 
# Extract coefficients and VCOV matrix
coef_est <- coef(Bio_Com_X_Turbine_Rd_Dense, type = "state")
vcov_mat <- vcov(Bio_Com_X_Turbine_Rd_Dense, type = "state")

beta_grassland <- coef_est["psi(turbine_rd_density_1_6km)"]
se_grassland <- sqrt(vcov_mat["psi(turbine_rd_density_1_6km)", 
                              "psi(turbine_rd_density_1_6km)"])

beta_woodland <- coef_est["psi(turbine_rd_density_1_6km)"] + 
  coef_est["psi(turbine_rd_density_1_6km:as.factor(biotic_com_2)2)"]

var_woodland <- 
  vcov_mat["psi(turbine_rd_density_1_6km)", "psi(turbine_rd_density_1_6km)"] +
  vcov_mat["psi(turbine_rd_density_1_6km:as.factor(biotic_com_2)2)",
           "psi(turbine_rd_density_1_6km:as.factor(biotic_com_2)2)"] +
  2 * vcov_mat["psi(turbine_rd_density_1_6km)", 
               "psi(turbine_rd_density_1_6km:as.factor(biotic_com_2)2)"]
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

# calculate SE for woodland 
mod <- Bio_Com_X_Turbine_Rd_Dense

b <- coef(mod)
V <- vcov(mod)

main  <- "psi(turbine_rd_density_1_6km)"
int   <- "psi(turbine_rd_density_1_6km:as.factor(biotic_com_2)2)"

# Grassland (reference)
est_g <- b[main]
se_g  <- sqrt(V[main, main])

# Woodland (main + interaction)
est_w <- b[main] + b[int]
se_w  <- sqrt(V[main, main] + V[int, int] + 2*V[main, int])

data.frame(
  Biotic_Community = c("Grassland","Woodland"),
  Estimate = c(est_g, est_w),
  SE = c(se_g, se_w)
)

# AICc Table with all possible turbine rd density models 

model_list <- list(
  "Vertical cover + slope" = Habitat,
  "Null" = Null,
  "Biotic community + slope + vertical cover" = Bio_Com_Cov_Slope,
  "Slope" = Slope_Occu,
  "Vertical cover" = Cov_Occu,
  "Turbine road density" = Turbine_Rd_Dense_Occu,
  "Turbine road density + vertical cover + slope" = Habitat_Add_Turbine_Rd_Dense,
  "Turbine road density + vertical cover" = Cov_Add_Turbine_Rd_Dense,
  "Turbine road density + slope" = Slope_Add_Turbine_Rd_Dense,
  "Turbine road density x biotic community + slope + vertical cover" =
    Bio_Com_X_Turbine_Rd_Dense
)

# Pair wind models with their matching non wind model

wind_pairs <- list(
  "Turbine road density" = "Null",
  "Turbine road density + vertical cover + slope" = "Vertical cover + slope",
  "Turbine road density + vertical cover" = "Vertical cover",
  "Turbine road density + slope" = "Slope",
  "Turbine road density x biotic community + slope + vertical cover" =
    "Biotic community + slope + vertical cover"
)

# Compile final table and export to .csv
final_table <- occu_model_selection(
  model_list = model_list,
  wind_pairs = wind_pairs,
  file_name = paste0(homewd, "outputs/badger_turbine_rd_density.csv")
)

# The beta values (ß), standard errors, and 85% confidence intervals for 
  #parameter estimates within the best-supported model describing the effect of
  #the density of access roads within a 1.6 km squared radius surrounding the 
  #site on the probability of habitat selection (ψ) for American badgers
  #are listed in Table S3.30.

########################################################
######################### END ##########################
########################################################