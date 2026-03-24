###########################################################
## Code for: Turbine visibility is a strong predictor of ##
## altered habitat selection by terrestrial mammals at a ##
######## wind energy facility in central New Mexico #######
###########################################################
##### This script is for the creation of the gray fox #####
#################### base models ##########################
###########################################################
################ Script by: Iona Rohan ####################
############ Contact: ionarohan12@gmail.com ###############
###########################################################
########## Date Last Modified: 23-March-2026 ##############
###########################################################

#Clear work environment
rm(list=ls())

#Note: If you opened this script through the .Rproj file, the only line you 
  #should need to change for the script to run (assuming packages are installed) 
  #is the homewd directory on line 24.

#Set home working directory
  #e.g. homewd = "C:/Users/ionar/Desktop/R Repository/Wind-energy-and-terrestrial-mammals/"
homewd = "<insert your folder here and end with a forward slash>"

#Set wd to data folder on your local computer 
setwd(paste0(homewd, "data/"))

# Load packages 
library(tidyverse)
library(unmarked)
library(MuMIn)
library(xlsx)
library(AICcmodavg)

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
detHist %>%
  mutate(across(1:51, as.numeric))

# Create detection history 
occu.urci <- unmarkedFrameOccu(y=detHist, 
                               siteCovs = site.covs.scaled, 
                               obsCovs = obsCovs.scaled)

##### Null model ####
urci.null <- occu( ~ 1 ~ 1, occu.urci, linkPsi="logit", se=TRUE,
                   control = list(maxit = 10000))
null.aicc <- AICc(urci.null)

# Null detection probability (Table 2)
backTransform(urci.null['det'])

### Calculating overall detection probability if K = 38
1-(1-0.104)^38

# 98% probability of detecting a gray fox at least once given the average
#sampling period length at a site and the null detection model (see Table 2)

# Null occupancy probability 
backTransform(urci.null['state'])

#################################################
# CREATE THE DETECTION MODEL #
#################################################

## Gray fox detection hypotheses are listed in Table S2.1.

## Table S2.3 lists which detection hypothesis groups proceeded (i.e. performed 
  # better than the null detection model p(.), did not have uninformative 
  # parameters, and supported our original hypothesis).

## Candidate detection models are found in Table S2.5.
  # Models were limited to 4 parameters total and were tested with no variables 
  # on the probability of habitat selection (ψ(.)).

##### Water Hypothesis Group ####

water.tank.dense <- occu( ~ water_tank_density_1_5_km 
                          ~ 1, occu.urci, starts = c(-1, 0, 0))    
null.aicc - AICc(water.tank.dense) #worse than null

water.tank.dist <- occu( ~ dist_water_tank
                         ~ 1, occu.urci, starts = c(-1, 0, 0))    
null.aicc - AICc(water.tank.dist) #better than null
confint(water.tank.dist, level=0.85, type="det") 
# 85% CI does not overlap zero but is positive (no support for hypothesis)

# none proceed

##### Precipitation Hypothesis Group ####

# Daily precipitation 

precip.cat <- occu( ~ as.factor(precip.cat) 
                    ~ 1, occu.urci, starts = c(-1, 0, 0))
null.aicc - AICc(precip.cat) #worse than null

precip.cm <- occu( ~ precip.cm ~ 1, occu.urci, starts = c(-1, 0, 0))
null.aicc - AICc(precip.cm) #worse than null

# Long-term precipitation

day.since.rain <- occu( ~ days.since.rain ~ 1, occu.urci, starts = c(-1, 0, 0)) 
null.aicc - AICc(day.since.rain, k=2) #worse than null

rain.month <- occu( ~ rain.month ~ 1, occu.urci, starts = c(-1, 0, 0))   
null.aicc - AICc(rain.month, k=2) #better than null
confint(rain.month, level=0.85, type="det") 
# 85% CI does not overlaps zero but is negative (no support for hypothesis)

rain.week <- occu( ~ rain.week ~ 1, occu.urci, starts = c(-1, 0, 0))    
null.aicc - AICc(rain.week, k=2) #worse than null

# Water source and long-term precipitation interaction

precip.water.dist <- occu( ~ water_tank_density_1_5_km * rain.month 
                           ~ 1, occu.urci)    
null.aicc - AICc(precip.water.dist, k=2) #worse than null 

# none proceed

#### Humidity Hypothesis Group ####

humidity <- occu( ~ humid ~ 1, occu.urci, starts = c(-1, 0, 0))    
null.aicc - AICc(humidity, k=2)  #worse than null

# none proceed

##### Temperature Hypothesis Group ####

temp.max <- occu( ~ max.temp ~ 1, occu.urci, starts = c(-1, 0, 0))       
null.aicc - AICc(temp.max, k=2) #worse than null

# Temperature interactions

temp.water.dense <- occu( ~ max.temp * water_tank_density_1_5_km ~ 1, occu.urci) 
null.aicc - AICc(temp.water.dense, k=2)  #worse than null

temp.canopy <- occu( ~ max.temp * canopy_cov ~ 1, occu.urci)      
null.aicc - AICc(temp.canopy, k=2)  #worse than null

temp.humid <- occu( ~ max.temp * humid ~ 1, occu.urci)      
null.aicc - AICc(temp.humid, k=2)  #worse than null

# none proceed

#### Wind Hypothesis Group ####

wind <- occu( ~ wind ~ 1, occu.urci, starts = c(-1, 0, 0))       
null.aicc - AICc(wind, k=2)  #worse than null

# none proceed

##### Predator Activity Hypothesis Group ####

# Related Univariates - choose best
coy.hours <- occu( ~ coyote.active ~ 1, occu.urci, starts = c(-1, 0, 0)) 
null.aicc - AICc(coy.hours,k=2)  #worse than null

coy.total <- occu( ~ coyote.count ~ 1, occu.urci, starts = c(-1, -1, 0)) 
null.aicc - AICc(coy.total,k=2)  #worse than null

coy.bob.hours <- occu( ~ coy.bob.active ~ 1, occu.urci, starts = c(-1, 0, 0)) 
null.aicc - AICc(coy.bob.hours,k=2)  #worse than null

coy.bob.total <- occu( ~ coy.bob.count ~ 1, occu.urci, starts = c(-1, -1, 0)) 
null.aicc - AICc(coy.bob.total,k=2) #better than null
confint(coy.bob.total, level=0.85, type="det") #85% CI overlaps zero

bob.hours <- occu( ~ bobcat.active ~ 1, occu.urci, starts = c(-1, -2, 0))     
#large SE for bobcat active
null.aicc - AICc(bob.hours,k=2)  #worse than null

bob.total <- occu( ~ bobcat.count ~ 1, occu.urci, starts = c(-1, -2, -2))  
#large SE for bobcat count
null.aicc - AICc(bob.total,k=2)  #worse than null

# Water sources and predator presence interactions (use best univariates)
coy.bob.total.water <- occu( ~ coy.bob.count * dist_water_tank ~ 1, occu.urci)
null.aicc - AICc(coy.bob.total.water) #better than null
confint(coy.bob.total.water, level=0.85, type="det") #85% CI overlaps zero

# none proceed

##### Livestock Activity Hypothesis Group ####

sheep.hours <- occu( ~ sheep.active ~ 1, occu.urci, starts = c(-1, 0, 0))  
null.aicc - AICc(sheep.hours, k=2) #worse than null

sheep.total <- occu( ~ sheep.count ~ 1, occu.urci, starts = c(-1, -2, -1))   
#large SE for sheep count
null.aicc - AICc(sheep.total,k=2) #worse than null

# none proceed

#### Prey Activity Hypothesis Group ####

lagomorph.total <- occu( ~ lagomorph.count ~ 1, occu.urci, starts = c(-1, 0, 0))
null.aicc - AICc(lagomorph.total,k=2) #worse than null

lagomorph.hours <- occu( ~ lagomorph.active ~ 1, occu.urci, 
                           starts = c(-1, -1, -1))     
null.aicc - AICc(lagomorph.hours,k=2) #worse than null

rodent.total <- occu( ~ rodent.count ~ 1, occu.urci, starts = c(-1, -2, -5))
#large SE for rodent count
null.aicc - AICc(rodent.total,k=2) #worse than null

rodent.hours <- occu( ~ rodent.active ~ 1, occu.urci, starts = c(-1, -2, -1))   
#large SE for rodent hours
null.aicc - AICc(rodent.hours,k=2) #worse than null

lago.rodent.total <- occu( ~ lago.rodent.count ~ 1, occu.urci, 
                             starts = c(-1, 0, 0))
null.aicc - AICc(lago.rodent.total,k=2) #worse than null

lago.rodent.hours <- occu( ~ lago.rodent.active ~ 1, occu.urci,  
                             starts = c(-1, -1, -1))     
null.aicc - AICc(lago.rodent.hours,k=2) #worse than null

cottontail.total <- occu( ~ cottontail.count ~ 1, occu.urci, 
                            starts = c(-1, 0, 0))
null.aicc - AICc(cottontail.total,k=2) #worse than null

cottontail.hours <- occu( ~ cottontail.active ~ 1, occu.urci,  
                            starts = c(-1, -1, -1))     
null.aicc - AICc(cottontail.hours,k=2) #worse than null

jackrabbit.total <- occu( ~ jackrabbit.count ~ 1, occu.urci, 
                            starts = c(-1, 0, 0))
null.aicc - AICc(jackrabbit.total,k=2) #worse than null

jackrabbit.hours <- occu( ~ jackrabbit.active ~ 1, occu.urci,  
                            starts = c(-1, -1, -1))     
null.aicc - AICc(jackrabbit.hours,k=2) #worse than null

#none proceed

##### Human Activity Hypothesis Group #####

vehicle.total <- occu( ~ vehicle.count ~ 1, occu.urci, starts = c(-1, 0, 0))
null.aicc - AICc(vehicle.total,k=2) #worse than null

vehicle.hours <- occu( ~ vehicle.active ~ 1, occu.urci, starts = c(-1, -1, -1))
null.aicc - AICc(vehicle.hours,k=2) #worse than null

human.hours <- occu( ~ human.active ~ 1, occu.urci, starts = c(-1, 0, 0)) 
null.aicc - AICc(human.hours,k=2) #worse than null

people.hours <- occu( ~ people.active ~ 1, occu.urci, starts = c(-1, 0, 0)) 
null.aicc - AICc(people.hours,k=2) #worse than null

##### Biotic Community Type Hypothesis Group ####

ndvi <- occu( ~ NDVI_1_5km ~ 1, occu.urci)
null.aicc - AICc(ndvi)  #worse than null

wood <- occu( ~ woodland_percent_1_5km ~ 1, occu.urci)          
null.aicc - AICc(wood, k=2) #worse than null

bio.com <- occu( ~ as.factor(biotic_com_2) ~ 1, occu.urci)        
null.aicc - AICc(bio.com, k=2)  #worse than null

tree <- occu( ~  tree_density_5_70 ~ 1, occu.urci)
null.aicc - AICc(tree)  #better than null
confint(tree, level=0.85, type="det") 
# 85% CI does not overlap zero but is negative (no support for hypothesis)

#none proceed

##### Vegetation Concealment Cover Hypothesis Group ####

obs.distance <- occu( ~ obs_dist ~ 1, occu.urci)       
null.aicc - AICc(obs.distance, k=2) #worse than null

veg.cov.cam <- occu( ~ veg_cover_cam ~ 1, occu.urci)       
null.aicc - AICc(veg.cov.cam, k=2) # worse than null

detect.angle <- occu( ~ detection_angle ~ 1, occu.urci, starts = c(-1, 0, 0))  
null.aicc - AICc(detect.angle, k=2) # worse than null

max.trig.dist <- occu( ~ max_trig_dist ~ 1, occu.urci, starts = c(-1, 0, 0))   
null.aicc - AICc(max.trig.dist, k=2) # worse than null

# none proceed

#### Sampling Period Hypothesis Group ####

### Sampling Period Length ###
samp.period <- occu( ~ sampling_period_length ~ 1, occu.urci)          
null.aicc - AICc(samp.period,k=2) # worse than null

### Camera Moved ###
cam.moved <- occu( ~ as.factor(cam_moved) ~ 1, occu.urci,starts = c(-1, 0, 0))
null.aicc - AICc(cam.moved, k=2) # better than null
confint(cam.moved, level=0.85, type="det") # 85% CI does not include zero

# cam moved proceeds

#################################################
# DETECTION MODEL SELECTION #
#################################################

# Variables that proceed (1)
  # cam_moved

# Null model
urci.null <- occu( ~ 1 ~ 1, occu.urci, linkPsi="logit", starts = c(-1, -1))

# Single variable model
mod1 <- occu( ~ as.factor(cam_moved) ~ 1, occu.urci, starts = c(-1, 0, 0))   

# Model selection
cand.models <- list(mod1, urci.null)

modnames <- c("mod1", "urci.null")

aicc_table <- aictab(cand.set = cand.models, modnames = modnames, sort = TRUE)
print(aicc_table)

# Candidate detection models are found in Table S2.5.

# Top model diagnostics

summary(mod1)
AICc(urci.null) - AICc(mod1)

# Calculate the 85% confidence intervals for variables 
confint(mod1, type = "det", level = 0.85)    
#CIs for detection variables - none overlap zero

confint(mod1, type = "state", level = 0.85)  
#CIs for occupancy intercept - does not overlap zero

#################################################
# CREATE THE HABITAT SELECTION (OCCUPANCY) MODEL#
#################################################

## Gray fox habitat selection hypotheses are listed in Table S2.2.

## Table S2.4 lists which occupancy hypothesis groups proceeded (i.e. performed 
# better than the null occupancy model ψ(.), did not have uninformative 
# parameters, and supported our original hypothesis).

## Candidate habitat selection (occupancy) models are found in Table S2.6.
  # Models were limited to 6 parameters and were derived from the best-supported
  # detection model (mod1 above).

occu.null.aicc <- AICc(mod1)

#### Water Hypothesis Group ####

water.tank.dense <- occu( ~ as.factor(cam_moved)
                          ~ water_tank_density_1_5_km, occu.urci)       
occu.null.aicc - AICc(water.tank.dense) #worse than null

water.dist.tank <- occu( ~ as.factor(cam_moved)
                         ~ dist_water_tank, occu.urci)      
occu.null.aicc - AICc(water.dist.tank, k=2) #worse than null 

#none proceed

#### Predator Activity Hypothesis Group ####

coy.count.avg <- occu( ~ as.factor(cam_moved)
                       ~ coy_count_avg, occu.urci)
occu.null.aicc - AICc(coy.count.avg, k=2) #better than null

bobcat.count.avg <- occu( ~ as.factor(cam_moved)
                          ~ bobcat_count_avg, occu.urci)
occu.null.aicc - AICc(bobcat.count.avg, k=2) #worse than null

# none proceed

#### Prey Activity Hypothesis Group ####

prey.count.avg <- occu( ~  as.factor(cam_moved)
                        ~  prey_count_avg, occu.urci)
occu.null.aicc - AICc(prey.count.avg, k=2) #worse than null

lagomorph.count.avg <- occu( ~  as.factor(cam_moved)
                             ~  lagomorph_count_avg, occu.urci)
occu.null.aicc - AICc(lagomorph.count.avg, k=2)  #worse than null

rodent.count.avg <- occu( ~  as.factor(cam_moved)
                          ~  rodent_count_avg, occu.urci)
occu.null.aicc - AICc(rodent.count.avg, k=2) #worse than null

jackrabbit.count.avg <- occu( ~ as.factor(cam_moved)
                              ~ jackrabbit_count_avg, occu.urci)
occu.null.aicc - AICc(jackrabbit.count.avg, k=2)  #worse than null

cottontail.count.avg <- occu( ~  as.factor(cam_moved)
                              ~  cottontail_count_avg, occu.urci)
occu.null.aicc - AICc(cottontail.count.avg, k=2)  #worse than null

# none proceed

#### Livestock Activity Hypothesis Group ####

sheep.count.avg <- occu( ~  as.factor(cam_moved)
                         ~  sheep_count_avg, occu.urci, 
                            starts = c(-2, -5, -2, -2)) 
                            # large SE for sheep count
occu.null.aicc - AICc(sheep.count.avg, k=2) #worse than null

# none proceed 

#### Biotic Community Type Hypothesis Group ####

tree <- occu( ~  as.factor(cam_moved) ~  tree_density_5_70, occu.urci)
occu.null.aicc - AICc(tree) #better than null
confint(tree, level=0.85, type="state") # 85% CI does not overlap zero

ndvi <- occu( ~ as.factor(cam_moved) ~ NDVI_1_5km, occu.urci)
occu.null.aicc - AICc(ndvi ) #better than null
confint(ndvi , level=0.85, type="state") # 85% CI does not overlap zero

woodland3 <- occu( ~ as.factor(cam_moved)
                   ~ woodland_percent_1_5km, occu.urci)          
occu.null.aicc - AICc(woodland3, k=2) #better than null
confint(woodland3, level=0.85, type="state") # 85% CI does not overlap zero

shrubland <- occu( ~  as.factor(cam_moved)
                   ~  shrub_yucca_density, occu.urci)
occu.null.aicc - AICc(shrubland) #better than null
confint(shrubland, level=0.85, type="state") # 85% CI does not overlap zero

grassland <- occu( ~ as.factor(cam_moved)
                   ~ herbaceous_cov, occu.urci)
occu.null.aicc - AICc(grassland) #better than null
confint(grassland, level=0.85, type="state")  # 85% CI does not overlap zero

bio.com <- occu(  ~ as.factor(cam_moved) 
                  ~ as.factor(biotic_com_2), occu.urci)        
occu.null.aicc - AICc(bio.com, k=2)  #better than null
confint(bio.com, level=0.85, type="state") # 85% CI does not overlap zero

# NDVI proceeds

#### Vegetation Concealment Cover Hypothesis Group  ####

vert.cover <- occu( ~  as.factor(cam_moved)
                    ~  vertical_cover, occu.urci)
occu.null.aicc - AICc(vert.cover)  #worse than null 

veg.cov <- occu( ~  as.factor(cam_moved)
                 ~  veg_cover, occu.urci)
occu.null.aicc - AICc(veg.cov)  #better than null
confint(veg.cov, level=0.85, type="state") # 85% CI does not overlap zero

# veg cover proceeds 

##### Topography Hypothesis Group ####

topo.pos <- occu( ~  as.factor(cam_moved)
                  ~  topo_pos, occu.urci)
occu.null.aicc - AICc(topo.pos) #better than null
confint(topo.pos, level=0.85, type="state") 
# 85% CI does not overlap zero but is negative (no support for hypothesis)

slope <- occu( ~ as.factor(cam_moved) ~ slope, occu.urci)
occu.null.aicc - AICc(slope) #better than null
confint(slope, level=0.85, type="state") # 85% CI does not overlap zero

elev <- occu( ~ as.factor(cam_moved) ~ elevation, occu.urci)
occu.null.aicc - AICc(elev) #better than null
confint(elev, level=0.85, type="state") # 85% CI does not overlap zero

# slope proceeds 

##### Coolness of Site Hypothesis Group ####

aspect <- occu( ~ as.factor(cam_moved) ~ aspect, occu.urci)
occu.null.aicc - AICc(aspect) #worse than null

heat.load <- occu( ~ as.factor(cam_moved) ~ heat_load, occu.urci)
occu.null.aicc - AICc(heat.load) #worse than null

#none proceed

#################################################
# OCCUPANCY MODEL SELECTION #
#################################################

# Variables that proceed (3)
  # NDVI 1.5 KM
  # slope
  # veg cover

# Check correlations between occupancy variables 

# Read in site-level covariates
site.covs <- read.csv("site_covs.csv", nrows = 102, header = TRUE)

site.covs.cor <- site.covs %>% 
  select(veg_cover, slope, NDVI_1_5km)

cor_site_covs <- cor(site.covs.cor, method='spearman') 
# none correlated above |0.7|

# Detection only model 
det.mod <- occu( ~ as.factor(cam_moved) ~ 1, occu.urci) 

# 1-variable occupancy models
mod1 <- occu( ~ as.factor(cam_moved) ~ slope, occu.urci)
mod2 <- occu( ~ as.factor(cam_moved) ~ NDVI_1_5km, occu.urci)
mod3 <- occu( ~ as.factor(cam_moved) ~ veg_cover, occu.urci)

# 2-variable occupancy models
mod4 <- occu( ~ as.factor(cam_moved) ~ slope + NDVI_1_5km, occu.urci)
mod5 <- occu( ~ as.factor(cam_moved) ~ slope + veg_cover, occu.urci)
mod6 <- occu( ~ as.factor(cam_moved) ~ NDVI_1_5km + veg_cover, occu.urci)

# Model selection
cand.models <- list(mod1, mod2, mod3, mod4, mod5, mod6, det.mod)

modnames <- c("mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "det.mod")

aicc_table <- aictab(cand.set = cand.models, modnames = modnames, sort = TRUE)
print(aicc_table)

# Candidate occupancy models are found in Table S2.6.

# Top model diagnostics

AICc(det.mod) - AICc(mod4) 
summary(mod4)

# The parameter estimates, standard errors, and 85% confidence intervals for
  #the best-supported model describing the probability of detection (p) and 
  #habitat selection (ψ) of gray foxes are listed in Table S2.10.

# Calculate the 85% confidence intervals for variables 
confint(mod4, type = "det", level = 0.85)    
#CIs for detection variables do not overlap zero

confint(mod4, type = "state", level = 0.85)  
#CIs for occupancy variables do not overlap zero

# Calculate multicollinearity
unmarked::vif(mod4, type = "state") # all below 2

# Calculate mean daily detection probability with the top model

backTransform(linearComb(mod4,                           
                         coefficients= c(1,0),  
                         type = 'det'))     
# 0.11 OR 11%%


# Calculate the occupancy probability

backTransform(linearComb(mod4,                
                         coefficients=c(1,0,0),   
                         type = 'state')) 
#  0.124 OR 12.4%

# Calculate the overall detection probability if K = 29
1-(1-0.11)^29

# 97% probability of detecting a gray fox at least once given the minimum
  #sampling period length at a site

# Calculate the overall detection probability if K = 38 
1-(1-0.11)^38

# 99% probability of detecting a gray fox at least once given the average
  #sampling period length at a site

# Predict detection probability for survey sites 

pred_det <- predict(mod4,          
                    type = "det",                 
                    newdata = occu.urci@siteCovs)[c("Predicted",
                                                    "SE",
                                                    "lower",
                                                    "upper")]

pred_det_df <- data.frame(Predicted = pred_det$Predicted,
                          StandardError = pred_det$SE,
                          lower = pred_det$lower,
                          upper = pred_det$upper,
                          site.covs)

# Predict occupancy probability for survey sites 

pred_occu <- predict(mod4,          
                     type = "state",                 
                     newdata = occu.urci@siteCovs)[c("Predicted",
                                                     "SE",
                                                     "lower",
                                                     "upper")]

pred_occu_df <- data.frame(Predicted = pred_occu$Predicted,
                           StandardError = pred_occu$SE,
                           lower = pred_occu$lower,
                           upper = pred_occu$upper,
                           site.covs)

########################################################
######################### END ##########################
########################################################