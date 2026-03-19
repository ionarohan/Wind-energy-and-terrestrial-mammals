###########################################################
##################### Pre-model code ######################
###########################################################
################ Script by: Iona Rohan ####################
############ Contact: ionarohan12@gmail.com ###############
###########################################################
########## Date Last Modified: 19-March-2026 ##############
###########################################################

###########################################################
##### NOTE: RUN THIS CODE BEFORE RUNNING ALL SPECIES' #####
########### BASE MODELS, WIND MODELS, OR PLOTS ############
###########################################################

#Clear work environment
rm(list=ls())

#Note: If you opened this script through the .Rproj file, the only line you 
  #should need to change for the script to run (assuming packages are installed) 
  #is the homewd directory on line 24

#Set home working directory
# e.g. homewd = "C:/Users/ionar/Desktop/R Repository/Wind-energy-and-terrestrial-mammals/"
homewd = "<insert your folder here and end with a forward slash>"

#Set wd to data folder on your local computer 
setwd(paste0(homewd, "data/"))

# Load packages 
library(dplyr)
library(lubridate)
library(AICcmodavg) 
library(MuMIn)
library(ggpubr)
library(unmarked)
library(lmtest)
library(xlsx)
library(ggplot2)
library(tidyverse)
library(scales)
library(patchwork)
library(cowplot)

# Input site-level covariates
site.covs <- read.csv("site_covs.csv", nrows = 102, header=TRUE)

# Function to scale site-level covariates 
scale_site_covs <- function(site_covs,
                            skip_covs = character(),
                            center = TRUE,
                            scale  = TRUE) {
  if (!is.data.frame(site_covs))
    stop("`site_covs` must be a data‑frame.")
  
  # function to decides whether to scale a column
  should_scale <- function(x, nm) {
    is.numeric(x) && !(nm %in% skip_covs)
  }
  
  # allocate vectors for means/SDs
  means <- numeric(length(site_covs))
  sds   <- numeric(length(site_covs))
  
  # scale in place
  site_covs_scaled <- site_covs
  for (i in seq_along(site_covs)) {
    col <- site_covs[[i]]
    nm  <- names(site_covs)[i]
    
    if (should_scale(col, nm)) {
      m <- if (center) mean(col, na.rm = TRUE) else 0
      s <- if (scale)  sd(col,   na.rm = TRUE) else 1
      site_covs_scaled[[i]] <- (col - m) / s
      means[i] <- m
      sds[i]   <- s
    } else {         
      means[i] <- NA
      sds[i]   <- NA
    }
  }
  
  stats <- data.frame(variable = names(site_covs),
                      mean     = means,
                      sd       = sds,
                      scaled   = !is.na(means),
                      row.names = NULL)
  
  list(site_covs_scaled = site_covs_scaled,
       stats            = stats)
}

# Skip categorical variables in scaling 
skip_vars <- c("site", "cam_number", "biotic_com_2", "cam_moved", "turbine_interior")

scaled_out <- scale_site_covs(site.covs, skip_covs = skip_vars)

site.covs.scaled <- scaled_out$site_covs_scaled   
scale.stats      <- scaled_out$stats             

# Create quadratic wind terms

site.covs.scaled$woodland_percent_2_4km2 <-
  site.covs.scaled$woodland_percent_2_4km^2

site.covs.scaled$woodland_percent_1_3km2 <- 
  site.covs.scaled$woodland_percent_1_3km^2

site.covs.scaled$shrub_yucca_density2 <- site.covs.scaled$shrub_yucca_density^2

site.covs.scaled$NDVI_2_4km2 <- site.covs.scaled$NDVI_2_4km^2

site.covs.scaled$NDVI_1_3km2 <- site.covs.scaled$NDVI_1_3km^2

site.covs.scaled$obs_dist2 <- site.covs.scaled$obs_dist^2

site.covs.scaled$veg_cover_cam_2 <- site.covs.scaled$veg_cover_cam^2

site.covs.scaled$veg_cover2 <- site.covs.scaled$veg_cover^2

site.covs.scaled$vertical_cover2 <- site.covs.scaled$vertical_cover^2

site.covs.scaled$detection_angle2 <- site.covs.scaled$detection_angle^2

site.covs.scaled$max_trig_dist2 <- site.covs.scaled$max_trig_dist^2

site.covs.scaled$tree_density_5_70_2 <- site.covs.scaled$tree_density_5_70^2

site.covs.scaled$herbaceous_cov2 <- site.covs.scaled$herbaceous_cov^2

#write site.covs.scaled
saveRDS(site.covs.scaled, "site_covs_scaled.RData")


# Insert observation-level covariates 

# Climate variables 
max.temp <- read.csv("max_temp.csv")
precip.cat <- read.csv("precip_cat.csv")  # Categorical
precip.cm <- read.csv("precip_cm.csv")
humid <- read.csv("humidity.csv")
wind <- read.csv("wind_speed.csv")
rain.month <- read.csv("rain.month.csv")
rain.week <- read.csv("rain.week.csv")
days.since.rain <- read.csv("days.since.rain.csv")

# Human variables 
vehicle.count <- read.csv("vehicle_count.csv")
vehicle.active <- read.csv("vehicle_hours.csv")
people.active <- read.csv("people_hours.csv")  # Only walkers (no vehicles)
human.active <- read.csv("human_hours.csv")   # Combo of vehicles and people

# Predator variables 
bobcat.active <- read.csv("bobcat_hours.csv")
bobcat.count <- read.csv("bobcat_count.csv")
coyote.count <- read.csv("coyote_count.csv")
coyote.active <- read.csv("coyote_hours.csv")
coy.bob.active <- read.csv("predator_hours.csv")
coy.bob.count <- read.csv("predator_count.csv")
all.pred.active <- read.csv("all_pred_hours.csv")
all.pred.count <- read.csv("all_pred_count.csv")
meso.active <- read.csv("meso_hours.csv") #foxes, badger, skunk
meso.count <- read.csv("meso_count.csv")  #foxes, badger, skunk

# Livestock variables 
cow.count <- read.csv("cow_count.csv")
cow.active <- read.csv("cow_hours.csv")       
sheep.count <- read.csv("sheep_count.csv")
sheep.active <- read.csv("sheep_hours.csv")
livestock.count <- read.csv("livestock_count.csv")
livestock.active <- read.csv("livestock_hours.csv")

# Prey Variables
lagomorph.active <- read.csv("lagomorph_hours.csv") 
lagomorph.count <- read.csv("lagomorph_count.csv")  
ungulate.active <- read.csv("ungulate_hours.csv")
ungulate.count <- read.csv("ungulate_count.csv")
jackrabbit.count <- read.csv("jackrabbit_count.csv")
jackrabbit.active <- read.csv("jackrabbit_hours.csv")
cottontail.count <- read.csv("cottontail_count.csv")
cottontail.active <- read.csv("cottontail_hours.csv")
rodent.count <- read.csv("rodent_count.csv")
rodent.active <- read.csv("rodent_hours.csv")
lago.rodent.active <- read.csv("lago_rodent_hours.csv") 
lago.rodent.count <- read.csv("lago_rodent_count.csv")  

# Camera Type
jackrabbit.cam <- read.csv("cam.type.jackrabbit.csv")
badger.cam <- read.csv("cam.type.null.csv")

# Create observation matrix  
obsCovs = list(max.temp = max.temp[,3:53], 
               precip.cm = precip.cm[,3:53],
               precip.cat = precip.cat[,3:53], 
               humid = humid[,3:53],
               wind = wind[,3:53], 
               rain.month = rain.month[,2:52], 
               rain.week = rain.week[,2:52],
               days.since.rain = days.since.rain[,2:52], 
               cow.count = cow.count[,3:53], 
               all.pred.count = all.pred.count[,3:53], 
               cow.active = cow.active[,3:53],
               all.pred.active = all.pred.active[,3:53], 
               vehicle.count = vehicle.count[,3:53],  
               vehicle.active = vehicle.active[,3:53], 
               people.active = people.active[,3:53],
               human.active = human.active[,3:53], 
               bobcat.active = bobcat.active[,3:53],
               bobcat.count = bobcat.count[,3:53], 
               coyote.count = coyote.count[,3:53],
               coyote.active = coyote.active[,3:53], 
               sheep.count = sheep.count[,3:53],
               sheep.active = sheep.active[,3:53], 
               livestock.count = livestock.count[,3:53],
               livestock.active = livestock.active[,3:53], 
               coy.bob.active = coy.bob.active[,3:53],
               coy.bob.count = coy.bob.count[,3:53], 
               meso.active = meso.active[,3:53],
               meso.count = meso.count[,3:53], 
               lagomorph.active = lagomorph.active[,3:53],
               lagomorph.count = lagomorph.count[,3:53], 
               lago.rodent.active = lago.rodent.active[,3:53],
               lago.rodent.count = lago.rodent.count[,3:53],
               ungulate.active = ungulate.active[,3:53], 
               ungulate.count = ungulate.count[,3:53], 
               jackrabbit.count = jackrabbit.count[,3:53], 
               jackrabbit.active = jackrabbit.active[,3:53], 
               cottontail.count = cottontail.count[,3:53], 
               cottontail.active = cottontail.active[,3:53],
               rodent.count = rodent.count[,3:53], 
               rodent.active = rodent.active[,3:53],
               jackrabbit.cam = jackrabbit.cam[,2:52], 
               badger.cam = badger.cam[,2:52]
)

#write obsCovs
saveRDS(obsCovs, "obsCovs.RData")

# Skip categorical observation-level covariates in scaling
categorical_covs <- c(
  "precip.cat", "jackrabbit.cam", "badger.cam"
)

# Function to scale a numeric site x occasion matrix
scale_matrix_manual_skip <- function(name, mat) {
  if (name %in% categorical_covs) return(mat)  # skip categorical variables
  m <- mean(unlist(mat), na.rm = TRUE)
  s <- sd(unlist(mat), na.rm = TRUE)
  return((mat - m) / s)
}

# Apply to all obsCovs using names
obsCovs.scaled <- mapply(scale_matrix_manual_skip, names(obsCovs), obsCovs, 
                         SIMPLIFY = FALSE)

#write  obsCovs.scaled
saveRDS(obsCovs.scaled, "obsCovs_scaled.RData")

# Store means and SDs for scaled covariates (exclude categorical ones)
obsCovs.stats <- lapply(obsCovs[setdiff(names(obsCovs), categorical_covs)], 
                        function(mat) {
                          list(mean = mean(unlist(mat), na.rm = TRUE),
                               sd   = sd(unlist(mat), na.rm = TRUE))
                        })

########################################################
######################### END ##########################
########################################################

