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
library(tidyverse)
library(unmarked)
library(MuMIn)
library(xlsx)

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

#### Null model ####
vuma.null <- occu( ~ 1 ~ 1, occu.vuma, linkPsi="logit", se=TRUE,
                  control = list(maxit = 10000))
null.aicc <- AICc(vuma.null)

# Null detection probability 
backTransform(vuma.null['det'])

# Null occupancy probability 
backTransform(vuma.null['state'])


#################################################
# CREATE THE DETECTION MODEL #
#################################################

## Kit fox detection hypotheses are listed in Table S2.1.

## Table S2.3 lists which detection hypothesis groups proceeded (i.e. performed 
  # better than the null detection model p(.), did not have uninformative 
  # parameters, and supported our original hypothesis).

## Candidate detection models are found in Table S2.5.
  # Models were limited to 4 parameters total and were tested with no variables 
  # on the probability of habitat selection (ψ(.)).

##### Water Hypothesis Group ####

water.tank.dense <- occu( ~ water_tank_density_1_9_km 
                          ~ 1, occu.vuma, starts = c(-1, 0, 0))       
null.aicc - AICc(water.tank.dense) #worse than null

water.dist.tank <- occu( ~ dist_water_tank ~ 1, occu.vuma, starts = c(-1, 0, 0))  
null.aicc - AICc(water.dist.tank, k=2) #worse than null

# none proceeds 

##### Precipitation Hypothesis Group####

# Daily precipitation 

precip.cat <- occu( ~ as.factor(precip.cat) ~ 1, occu.vuma, starts = c(-1, 0, 0))
null.aicc - AICc(precip.cat) #better than null
confint(precip.cat, level=0.85, type="det") #85% CI does not overlap zero

precip.cm <- occu( ~ precip.cm ~ 1, occu.vuma, starts = c(-1, 0, 0))
null.aicc - AICc(precip.cm) #worse than null

# Long-term precipitation

days.since.rain <- occu( ~ days.since.rain ~ 1, occu.vuma, starts = c(-1, 0, 0))    
null.aicc - AICc(days.since.rain, k=2) #better than null
confint(days.since.rain, level=0.85, type="det") #85% CI does not overlap zero

rain.month <- occu( ~ rain.month ~ 1, occu.vuma, starts = c(-1, 0, 0))    
null.aicc - AICc(rain.month, k=2) #better than null
confint(rain.month, level=0.85, type="det") #85% CI does not overlap zero

rain.week <- occu( ~ rain.week ~ 1, occu.vuma, starts = c(-1, 0, 0))    
null.aicc - AICc(rain.week, k=2) #worse than null

# Water source and long-term precipitation interaction

precip.water.dense <- occu( ~ water_tank_density_1_9_km * rain.month 
                            ~ 1, occu.vuma)    
null.aicc - AICc(precip.water.dense, k=2) #worse than null 

# precipitation categorical proceeds 

#### Humidity Hypothesis Group ####
humidity <- occu( ~ humid ~ 1, occu.vuma, starts = c(-1, 0, 0))       
null.aicc - AICc(humidity, k=2) #worse than null

# none proceed

##### Temperature Hypothesis Group ####

temp.max <- occu( ~ max.temp ~ 1, occu.vuma, starts = c(-1, 0, 0))       
null.aicc - AICc(temp.max, k=2) #better than null
confint(temp.max, level=0.85, type="det")  #85% CI does not overlap zero

# Temperature interactions

water.dense <- occu( ~ max.temp * water_tank_density_1_9_km ~ 1, occu.vuma)       
null.aicc - AICc(water.dense, k=2) 
confint(water.dense, level=0.85, type="det")  #85% CI overlaps zero

canopy <- occu( ~ max.temp * canopy_cov ~ 1, occu.vuma)      
null.aicc - AICc(canopy, k=2) #better than null
confint(canopy, level=0.85, type="det") #85% CI overlaps zero

temp.humid <- occu( ~ max.temp * humid ~ 1, occu.vuma)      
null.aicc - AICc(temp.humid, k=2) #worse than null

# max temp proceeds  

#### Wind Hypothesis Group ####

wind <- occu( ~ wind ~ 1, occu.vuma, starts = c(-1, 0, 0))       
null.aicc - AICc(wind, k=2) #better than null
confint(wind, level=0.85, type="det") #85% CI does not overlap zero

# wind proceeds

##### Predator Activity Hypothesis Group ####

coy.hours <- occu( ~ coyote.active ~ 1, occu.vuma, starts = c(0,0,0)) 
null.aicc - AICc(coy.hours,k=2) #worse than null

coy.total <- occu( ~ coyote.count ~ 1, occu.vuma, starts = c(-1, -1, 0)) 
null.aicc - AICc(coy.total,k=2) #worse than null

coy.bob.hours <- occu( ~ coy.bob.active ~ 1, occu.vuma, starts = c(-1, 0, 0)) 
null.aicc - AICc(coy.bob.hours,k=2) #worse than null

coy.bob.total <- occu( ~ coy.bob.count ~ 1, occu.vuma) 
null.aicc - AICc(coy.bob.total,k=2) #worse than null

bob.hours <- occu( ~ bobcat.active ~ 1, occu.vuma)        
null.aicc - AICc(bob.hours,k=2)  #worse than null 

bob.total <- occu( ~ bobcat.count ~ 1, occu.vuma)        
null.aicc - AICc(bob.total,k=2) #worse than null

# Predator interactions

bob.hours.water.dense <- occu( ~ coy.bob.count * water_tank_density_1_9_km
                               ~ 1, occu.vuma, starts = c(-1, 0, 0, 0, 0)) 
null.aicc - AICc(bob.hours.water.dense,k=2) #worse than null

# none proceed

##### Livestock Activity Hypothesis Group####

sheep.hours <- occu( ~ sheep.active ~ 1, occu.vuma, starts = c(0, 0, 0))       
null.aicc - AICc(sheep.hours, k=2) #better than null
confint(sheep.hours, level=0.85, type="det") 
# 85% CI does not overlaps zero but is positive (no support for hypothesis)

sheep.total <- occu( ~ sheep.count ~ 1, occu.vuma,  starts = c(0,0,0))       
null.aicc - AICc(sheep.total,k=2) #worse than null

# none proceeds

#### Prey Activity Hypothesis Group ####

lagomorph.total <- occu( ~ lagomorph.count ~ 1, occu.vuma, starts = c(-1, 0, 0))
null.aicc - AICc(lagomorph.total,k=2) #worse than null

lagomorph.hours <- occu( ~ lagomorph.active 
                         ~ 1, occu.vuma, starts = c(-1, -1, -1))     
null.aicc - AICc(lagomorph.hours,k=2) #worse than null

rodent.total <- occu( ~ rodent.count ~ 1, occu.vuma,  starts = c(-1, -1, 0))
null.aicc - AICc(rodent.total,k=2) #worse than null 

rodent.hours <- occu( ~ rodent.active ~ 1, occu.vuma, starts = c(-1, -1, -1))     
null.aicc - AICc(rodent.hours,k=2) #worse than null

lago.rodent.total <- occu( ~ lago.rodent.count 
                           ~ 1, occu.vuma, starts = c(-1, 0, 0))
null.aicc - AICc(lago.rodent.total,k=2) #worse than null

lago.rodent.hours <- occu( ~ lago.rodent.active 
                           ~ 1, occu.vuma, starts = c(-1, -1, -1))     
null.aicc - AICc(lago.rodent.hours,k=2) #worse than null

cottontail.total <- occu( ~ cottontail.count
                          ~ 1, occu.vuma, starts = c(-1, 0, 0))
null.aicc - AICc(cottontail.total,k=2) #worse than null

cottontail.hours <- occu( ~ cottontail.active 
                          ~ 1, occu.vuma, starts = c(-1, -1, -1))     
null.aicc - AICc(cottontail.hours,k=2) #worse than null

jackrabbit.total <- occu( ~ jackrabbit.count 
                          ~ 1, occu.vuma,  starts = c(-1, 0, 0))
null.aicc - AICc(jackrabbit.total,k=2) #worse than null

jackrabbit.hours <- occu( ~ jackrabbit.active 
                          ~ 1, occu.vuma, starts = c(-1, -1, -1))     
null.aicc - AICc(jackrabbit.hours,k=2) #worse than null

# none proceed

##### Human Activity Hypothesis Group#####

vehicle.total <- occu( ~ vehicle.count ~ 1, occu.vuma, starts = c(-1, 0, 0))
null.aicc - AICc(vehicle.total,k=2) # worse than null

vehicle.hours <- occu( ~ vehicle.active ~ 1, occu.vuma, starts = c(-1, -1, -1))
null.aicc - AICc(vehicle.hours,k=2) # worse than null

human.hours <- occu( ~ human.active ~ 1, occu.vuma, starts = c(-1, 0, 0)) 
null.aicc - AICc(human.hours,k=2) #better than null
confint(human.hours, level=0.85, type="det") # 85% CI overlaps zero

people.hours <- occu( ~ people.active ~ 1, occu.vuma, starts = c(-1, 0, 0)) 
null.aicc - AICc(people.hours,k=2) #better than null
confint(people.hours, level=0.85, type="det") # 85% CI overlaps zero

# none proceeds

#### Biotic Community Type Hypothesis Group ####

ndvi.site <- occu( ~ NDVI_1_9km ~ 1, occu.vuma, starts = c(-1, 0, 0))           
null.aicc - AICc(ndvi.site,k=2) #better than null
confint(ndvi.site, level=0.85, type="det") # 85% CI does not overlap zero 

wood <- occu( ~ woodland_percent_1_9km ~ 1, occu.vuma, starts = c(-1, 0, 0))
null.aicc - AICc(wood, k=2) #better than null
confint(wood, level=0.85, type="det") # 85% CI does not overlap zero 

bio.com <- occu( ~ as.factor(biotic_com_2) ~ 1, occu.vuma, 
                   starts = c(0, -4, -10))        
null.aicc - AICc(bio.com, k=2)  #better than null
confint(bio.com, level=0.85, type="det") # 85% CI overlaps zero 

tree <- occu( ~  tree_density_5_70 ~ 1, occu.vuma, starts = c(0, -8, -7))
null.aicc - AICc(tree)  #better than null
confint(tree, level=0.85, type="det") # 85% CI does not overlap zero 

## NDVI proceeds

##### Vegetation Concealment Cover Hypothesis Group ####

obs.distance <- occu( ~ obs_dist ~ 1, occu.vuma)       
null.aicc - AICc(obs.distance, k=2) # worse than null

veg.cov.cam <- occu( ~ veg_cover_cam ~ 1, occu.vuma)       
null.aicc - AICc(veg.cov.cam, k=2) #better than null
confint(veg.cov.cam, level=0.85, type="det") 
# 85% CI does not overlaps zero but is positive (no support for hypothesis)

detect.angle <- occu( ~ detection_angle ~ 1, occu.vuma, starts = c(-1, 0, 0))  
null.aicc - AICc(detect.angle, k=2) #better than null
confint(detect.angle, level=0.85, type="det") # 85% CI does not overlap zero 

max.trig.dist <- occu( ~ max_trig_dist ~ 1, occu.vuma, starts = c(-1, 0, 0))   
null.aicc - AICc(max.trig.dist, k=2) # worse than null

# detection angle proceeds

#### Sampling Period Hypothesis Group ####

### Sampling Period Length ###
samp.period <- occu( ~ sampling_period_length ~ 1, occu.vuma)          
null.aicc - AICc(samp.period,k=2) #better than null
confint(samp.period, level=0.85, type="det") 
# 85% CI does not overlaps zero but is negative (no support for hypothesis)

### Camera Moved ###
cam.moved <- occu( ~ as.factor(cam_moved) ~ 1, occu.vuma,starts = c(-1, 0, 0))
null.aicc - AICc(cam.moved, k=2) #better than null
confint(cam.moved, level=0.85, type="det") # 85% CI does not overlap zero 

# cam moved proceeds

#################################################
# DETECTION MODEL SELECTION #
#################################################

# Variables that proceed (5)
  #NDVI
  #detection angle
  #camera moved
  #precipitation
  #max temp

# Check correlations between detection variables 

# Read in site-level covariates
site.covs <- read.csv("site_covs.csv", nrows = 102, header = TRUE)

# site-level variables 
site.covs.cor <- site.covs %>% 
  select(NDVI_1_9km, detection_angle)

cor_site_covs <- cor(site.covs.cor, method='spearman') 
# none correlated above |0.7|

# Read in observation-level covariates
obs.covs <- readRDS("obsCovs.RData")

#  observation-level variables
cor_obs <- cor(obs.covs $precip.cat, obs.covs$max.temp)

# Run this code to get correlation matrix in Excel 
  #write.xlsx(cor_obs, file="Correlations kit fox.xlsx")

# Null model
vuma.null <- occu( ~ 1 ~ 1, occu.vuma)

# single-variable models
mod1 <- occu( ~ detection_angle ~ 1, occu.vuma)
mod2 <- occu( ~ as.factor(precip.cat) ~ 1, occu.vuma)
mod3 <- occu( ~ max.temp ~ 1, occu.vuma)
mod4 <- occu( ~ as.factor(cam_moved) ~ 1, occu.vuma)
mod5 <- occu( ~ NDVI_1_9km ~ 1, occu.vuma)
mod6 <- occu( ~ wind ~ 1, occu.vuma)

#2-variable models
mod7 <- occu( ~ detection_angle + max.temp ~ 1, occu.vuma)
mod8 <- occu( ~ detection_angle + as.factor(cam_moved) ~ 1, occu.vuma)
mod9 <- occu( ~ detection_angle + NDVI_1_9km ~ 1, occu.vuma)
mod10 <- occu( ~ as.factor(precip.cat) + max.temp ~ 1, occu.vuma)
mod11 <- occu( ~ as.factor(precip.cat) + as.factor(cam_moved) ~ 1, occu.vuma)
mod12 <- occu( ~ as.factor(precip.cat) + NDVI_1_9km ~ 1, occu.vuma)
mod13 <- occu( ~ max.temp + as.factor(cam_moved) ~ 1, occu.vuma)
mod14 <- occu( ~ max.temp + NDVI_1_9km ~ 1, occu.vuma)
mod15 <- occu( ~ as.factor(cam_moved) + NDVI_1_9km ~ 1, occu.vuma)
mod16 <- occu( ~ detection_angle + as.factor(precip.cat) ~ 1, occu.vuma)
mod17 <- occu( ~ wind + as.factor(precip.cat) ~ 1, occu.vuma)
mod18 <- occu( ~ wind + max.temp ~ 1, occu.vuma)
mod19 <- occu( ~ wind + as.factor(cam_moved) ~ 1, occu.vuma)
mod20 <- occu( ~ wind + NDVI_1_9km ~ 1, occu.vuma)
mod21 <- occu( ~ detection_angle + wind ~ 1, occu.vuma)

# Model selection
top_mods <- model.sel(
  mod1, mod2, mod3, mod4, mod5, mod6, mod7,
  mod8, mod9, mod10, mod11, mod12, mod13, mod14, mod15, 
  mod16, mod17, mod18, mod19, mod20, mod21, vuma.null
)

# Candidate detection models are found in Table S2.5.

# Run this code to see candidate detection models 
#write.xlsx(top_mods, file="Kit Fox Base Models.xlsx", 
       #  sheetName="Detection Models", append=T)  

# Top model diagnostics

summary(mod15)
AICc(vuma.null) - AICc(mod15)

# Calculate the 85% confidence intervals for variables 
confint(mod15, type = "det", level = 0.85)    
#CIs for detection variables - none overlap zero

confint(mod15, type = "state", level = 0.85)  
#CIs for occupancy intercept - overlaps zero

# Calculate multicollinearity
unmarked::vif(mod15, type = "det") # no colliniarity - all below 2

#################################################
# CREATE THE HABITAT SELECTION (OCCUPANCY) MODEL#
#################################################

## Kit fox habitat selection hypotheses are listed in Table S2.2.

## Table S2.4 lists which occupancy hypothesis groups proceeded (i.e. performed 
  # better than the null occupancy model ψ(.), did not have uninformative 
  # parameters, and supported our original hypothesis).

## Candidate habitat selection (occupancy) models are found in Table S2.6.
  # Models were limited to 6 parameters and were derived from the best-supported
  # detection model (mod1 above).

occu.null.aicc <- AICc(mod15)

####Water Hypothesis Group####

water.tank.dense <- occu( ~  as.factor(cam_moved) + NDVI_1_9km 
                          ~ water_tank_density_1_9_km, occu.vuma)       
occu.null.aicc - AICc(water.tank.dense) #worse than null

water.dist.tank <- occu( ~ as.factor(cam_moved) + NDVI_1_9km 
                         ~ dist_water_tank, occu.vuma)      
occu.null.aicc - AICc(water.dist.tank, k=2) #worse than null

# none proceed 

#### Predator Activity Hypothesis Group####

coyote <- occu( ~  as.factor(cam_moved) + NDVI_1_9km 
                ~  coy_count_avg, occu.vuma)
occu.null.aicc - AICc(coyote, k=2) #better than null
confint(coyote, level=0.85, type="state") # 85% CI does not overlap zero

#  coy count proceeds 

#### Prey Activity Hypothesis Group####

prey <- occu( ~ as.factor(cam_moved) + NDVI_1_9km 
              ~ prey_count_avg, occu.vuma)
occu.null.aicc - AICc(prey, k=2) #better than null
confint(prey, level=0.85, type="state")  # 85% CI does not overlap zero

lagomorph <- occu( ~ as.factor(cam_moved) + NDVI_1_9km 
                   ~ lagomorph_count_avg, occu.vuma)
occu.null.aicc - AICc(lagomorph, k=2) #better than null
confint(lagomorph, level=0.85, type="state")  # 85% CI does not overlap zero

rodent <- occu( ~ as.factor(cam_moved) + NDVI_1_9km 
                ~ rodent_count_avg, occu.vuma)
occu.null.aicc - AICc(rodent, k=2) #worse than null

jackrabbit <- occu( ~ as.factor(cam_moved) + NDVI_1_9km 
                    ~ jackrabbit_count_avg, occu.vuma)
occu.null.aicc - AICc(jackrabbit, k=2) #better than null
confint(jackrabbit, level=0.85, type="state")  # 85% CI does not overlap zero

cottontail <- occu( ~ as.factor(cam_moved) + NDVI_1_9km 
                    ~ cottontail_count_avg, occu.vuma) 
occu.null.aicc - AICc(cottontail, k=2) #worse than null

# jackrabbit count proceeds

#### Livestock Activity Hypothesis Group####

sheep <- occu( ~ as.factor(cam_moved) + NDVI_1_9km 
               ~ sheep_count_avg, occu.vuma, starts = c(-1, -1, -1, -1, -1))
occu.null.aicc - AICc(sheep, k=2) #better than null
confint(sheep, level=0.85, type="state") # 85% CI overlaps zero

# none proceed 

#### Biotic Community Type Hypothesis Group ####

tree <- occu( ~ as.factor(cam_moved) + NDVI_1_9km 
              ~ tree_density_5_70, occu.vuma)
occu.null.aicc - AICc(tree) #better than null
confint(tree, level=0.85, type="state") # 85% CI overlaps zero

ndvi <- occu( ~  as.factor(cam_moved) + NDVI_1_9km 
              ~  NDVI_1_9km, occu.vuma)
occu.null.aicc - AICc(ndvi) #worse than null

shrubland <- occu( ~ as.factor(cam_moved) + NDVI_1_9km 
                   ~ shrub_yucca_density, occu.vuma)
occu.null.aicc - AICc(shrubland) #worse than null

grassland <- occu( ~  as.factor(cam_moved) + NDVI_1_9km 
                   ~  herbaceous_cov, occu.vuma)
occu.null.aicc - AICc(grassland) #worse than null

wood <- occu( ~ as.factor(cam_moved) + NDVI_1_9km 
              ~ woodland_percent_1_9km, occu.vuma)          
occu.null.aicc - AICc(wood, k=2)  #worse than null

bio.com <- occu( ~  as.factor(cam_moved) + NDVI_1_9km 
                 ~  as.factor(biotic_com_2), occu.vuma, starts = c(0,0,0,1,0.5))        
occu.null.aicc - AICc(bio.com)  #better than null
confint(bio.com, level=0.85, type="state") # 85% CI overlaps zero

# none proceed 

#### Vegetation Concealment Cover Hypothesis Group  ####

vert.cover <- occu( ~ as.factor(cam_moved) + NDVI_1_9km 
                    ~ vertical_cover, occu.vuma)
occu.null.aicc - AICc(vert.cover) #worse than null

veg.cov <- occu( ~  as.factor(cam_moved) + NDVI_1_9km 
                 ~  veg_cover, occu.vuma)
occu.null.aicc - AICc(veg.cov) #better than null
confint(veg.cov, level=0.85, type="state")  # 85% CI does not overlap zero

# veg cover proceeds

##### Topography Hypothesis Group ####

topo.pos <- occu( ~  as.factor(cam_moved) + NDVI_1_9km 
                  ~  topo_pos, occu.vuma)
occu.null.aicc - AICc(topo.pos) #worse than null

slope <- occu( ~  as.factor(cam_moved) + NDVI_1_9km 
               ~  slope, occu.vuma)
occu.null.aicc - AICc(slope) #worse than null

elev <- occu( ~  as.factor(cam_moved) + NDVI_1_9km 
              ~  elevation, occu.vuma)
occu.null.aicc - AICc(elev) #better than null
confint(elev, level=0.85, type="state")  
# 85% CI does not overlap zero but is positive (no support for hypothesis)

# none proceeds

##### Coolness of Site Hypothesis Group ####

aspect <- occu( ~ as.factor(cam_moved) + NDVI_1_9km ~ aspect, occu.vuma)
occu.null.aicc - AICc(aspect) #better than null
confint(aspect, level=0.85, type="state") 
# 85% CI does not overlap zero but is positive (no support for hypothesis)

heat.load <- occu( ~  as.factor(cam_moved) + NDVI_1_9km ~  heat_load, occu.vuma)
occu.null.aicc - AICc(heat.load) #worse than null

#none proceeds

#################################################
# OCCUPANCY MODEL SELECTION #
#################################################

# Variables that proceed (3)
  # veg cover
  # jackrabbit count
  # coyote count

# Check correlations between occupancy variables 
site.covs.cor <- site.covs %>% 
  select(coy_count_avg, jackrabbit_count_avg, veg_cover)

cor_site_covs <- cor(site.covs.cor, method='spearman')
# none correlated above |0.7|

# Detection only model 
det.mod <- occu( ~ as.factor(cam_moved) + NDVI_1_9km ~ 1, occu.vuma) 

#single variable mods
mod1 <- occu( ~ as.factor(cam_moved) + NDVI_1_9km 
              ~ coy_count_avg, occu.vuma)
mod2 <- occu( ~ as.factor(cam_moved) + NDVI_1_9km 
              ~ jackrabbit_count_avg, occu.vuma)
mod3 <- occu( ~ as.factor(cam_moved) + NDVI_1_9km 
              ~ veg_cover, occu.vuma)

#2-variable mods
mod4 <- occu( ~  as.factor(cam_moved) + NDVI_1_9km 
              ~ coy_count_avg + jackrabbit_count_avg, occu.vuma, 
                starts = c(-1, -1, -1, -1, -1, -1))
mod5 <- occu( ~  as.factor(cam_moved) + NDVI_1_9km 
              ~ coy_count_avg + veg_cover, occu.vuma)
mod6 <- occu( ~  as.factor(cam_moved) + NDVI_1_9km 
              ~ jackrabbit_count_avg + veg_cover, occu.vuma)

# Combine models in model selection table
top_mods <- model.sel(
  mod1, mod2, mod3, mod4, mod5, mod6, det.mod)

# Candidate occupancy models are found in Table S2.6.

# Run this code to see candidate occupancy models 
  #write.xlsx(top_mods, file="Kit Fox Base Models.xlsx", 
        #  sheetName="Occupancy Models", append=T)  

# Top model diagnostics

summary(mod4)
AICc(det.mod) - AICc(mod4)

# The parameter estimates, standard errors, and 85% confidence intervals for
  #the best-supported model describing the probability of detection (p) and 
  #habitat selection (ψ) of kit foxes are listed in Table S2.11.

### Calculate the 85% confidence intervals for variables 
confint(mod4, type = "det", level = 0.85)    
#CIs for detection variables do not overlap zero

confint(mod4, type = "state", level = 0.85) 
#CIs for occupancy variables do not overlap zero

# Calculate mean daily detection probability with the top detection model

backTransform(linearComb(mod4,                           
                         coefficients= c(1,0,0),  
                         type = 'det'))     
# 0.00533 0.53%

# Calculate the overall detection probability if K = 29
1-(1-0.00533)^29

# 14% probability of detecting a kit fox at least once given the minimum
  #sampling period length at a site

# Calculate the overall detection probability if K = 38 
1-(1-0.00533)^38

# 18% probability of detecting a kit fox at least once given the average
  #sampling period length at a site

### Occupancy probability

backTransform(linearComb(mod4,                
                         coefficients=c(1,0,0),   
                         type = 'state')) 
# 0.181   OR 18.1%

# Calculate multicollinearity
unmarked::vif(mod4, type = "det") # no colliniarity - all below 2
unmarked::vif(mod4, type = "state") # no colliniarity - all below 4

#Predict detection probability for survey sites 

pred_det <- predict(mod4,          
                    type = "det",                 
                    newdata = occu.vuma@siteCovs)[c("Predicted",
                                                    "SE",
                                                    "lower",
                                                    "upper")]

pred_det_df <- data.frame(Predicted = pred_det$Predicted,
                          StandardError = pred_det$SE,
                          lower = pred_det$lower,
                          upper = pred_det$upper,
                          site.covs)

#Predict occupancy probability for survey sites 

pred_occu <- predict(mod4,          
                     type = "state",                 
                     newdata = occu.vuma@siteCovs)[c("Predicted",
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