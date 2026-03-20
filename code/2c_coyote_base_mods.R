###########################################################
############ Code for: base models for coyotes ############
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
# SETUP CODE FOR COYOTE OCCUPANCY MODELS #
###########################################################

# Read in edited .csv
detHist <- read.csv(file = "coyote detection hist.csv", row.names = 1)
# Read in site.covs.scaled
site.covs.scaled <- readRDS("site_covs_scaled.RData")
# Read in obsCovs.scaled
obsCovs.scaled <- readRDS("obsCovs_scaled.RData")

# Change from integer to numeric
detHist %>%
  mutate(across(1:51, as.numeric))

# Create detection history 
occu.cala <- unmarkedFrameOccu(y=detHist, 
                               siteCovs = site.covs.scaled, 
                               obsCovs = obsCovs.scaled)

#### Null model ####
cala.null <- occu( ~ 1 ~ 1, occu.cala, linkPsi="logit", starts = c(2, -3),
                   se = TRUE, control = list(maxit = 10000))

null.aicc <- AICc(cala.null)

# Null detection probability 
backTransform(cala.null['det'])

# Null occupancy probability 
backTransform(cala.null['state'])

#################################################
# CREATE THE DETECTION MODEL #
#################################################

## Coyote detection hypotheses are listed in Table S2.1.

## Table S2.3 lists which detection hypothesis groups proceeded (i.e. performed 
  # better than the null detection model p(.), did not have uninformative 
  # parameters, and supported our original hypothesis).

## Candidate detection models are found in Table S2.5.
  # Models were limited to 4 parameters total and were tested with no variables 
  # on the probability of habitat selection (ψ(.)).

##### Water Hypothesis Group ####

water.tank.dense <- occu( ~ water_tank_density_2_4_km 
                          ~ 1, occu.cala, starts = c(-1, 0, 0))       
null.aicc - AICc(water.tank.dense) #worse than null

water.dist.tank <- occu( ~ dist_water_tank ~ 1, occu.cala, 
                         starts = c(-1, -1, 0))      
null.aicc - AICc(water.dist.tank, k=2) #worse than null

# none proceed

##### Precipitation Hypothesis Group####

# Daily precipitation 

precip.cat <- occu( ~ as.factor(precip.cat) 
                    ~ 1, occu.cala, starts = c(-1, -1, 0))
null.aicc - AICc(precip.cat) #better than null
confint(precip.cat, level=0.85, type="det") 
# 85% CI does not overlaps zero but is positive (no support for hypothesis)

precip.cm <- occu( ~ precip.cm ~ 1, occu.cala, starts = c(-1, -1, -1))
null.aicc - AICc(precip.cm) #worse than null

# Long-term precipitation

days.since.rain <- occu( ~ days.since.rain ~ 1, occu.cala, starts = c(-1, 0, 0))    
null.aicc - AICc(days.since.rain, k=2) #worse than null

rain.month <- occu( ~ rain.month ~ 1, occu.cala, starts = c(2, 0, -3))    
null.aicc - AICc(rain.month, k=2) #worse than null

rain.week <- occu( ~ rain.week ~ 1, occu.cala, starts = c(-1, -1, -1))    
null.aicc - AICc(rain.week, k=2) #worse than null

# Water source and long-term precipitation interaction

precip.water.dense <- occu( ~ water_tank_density_2_4_km * days.since.rain 
                            ~ 1, occu.cala)       
null.aicc - AICc(precip.water.dense, k=2) #better than null
confint(precip.water.dense, level=0.85, type="det") #85% CI overlaps zero

# none proceeds

#### Humidity Hypothesis Group ####

humidity <- occu( ~ humid ~ 1, occu.cala, starts = c(2, 0, -3))       
null.aicc - AICc(humidity, k=2) #worse than null

# none proceed

##### Temperature Hypothesis Group ####

temp.max <- occu( ~ max.temp ~ 1, occu.cala, starts = c(-1, 0, 0))       
null.aicc - AICc(temp.max, k=2) #worse than null

# Temperature interactions

temp.water <- occu( ~ max.temp * dist_water_tank ~ 1, occu.cala)       
null.aicc - AICc(temp.water, k=2) #worse than null

temp.humid <- occu( ~ max.temp * humid ~ 1, occu.cala)      
null.aicc - AICc(temp.humid, k=2) #worse than null

temp.canopy <- occu( ~ max.temp * canopy_cov ~ 1, occu.cala)      
null.aicc - AICc(temp.canopy, k=2)  #worse than null

#none proceed

#### Wind Hypothesis Group ####

wind <- occu( ~ wind ~ 1, occu.cala, starts = c(-1, 0, 0))       
null.aicc - AICc(wind, k=2) #worse than null

#none proceed

##### Livestock Activity Hypothesis Group####

cow.hours <- occu( ~ cow.active ~ 1, occu.cala, starts = c(-1, 0, 0))       
null.aicc - AICc(cow.hours) #better than null
confint(cow.hours, level=0.85, type="det") # 85% CI does not overlap zero

cow.total <- occu( ~ cow.count ~ 1, occu.cala, starts = c(-1, 0, 0))       
null.aicc - AICc(cow.total,k=2) #worse than null

sheep.hours <- occu( ~ sheep.active ~ 1, occu.cala, starts = c(2, 0, -3))       
null.aicc - AICc(sheep.hours, k=2) #worse than null

sheep.total <- occu( ~ sheep.count ~ 1, occu.cala, starts = c(0, -1, 0))       
null.aicc - AICc(sheep.total,k=2) #worse than null

stock.total <- occu( ~ livestock.count ~ 1, occu.cala, starts = c(-1, 0, 0))
null.aicc - AICc(stock.total,k=2)  #worse than null

stock.hours <- occu( ~ livestock.active ~ 1, occu.cala, starts = c(-1, 0, 0))     
null.aicc - AICc(stock.hours,k=2) #better than null
confint(stock.hours, level=0.85, type="det")  # 85% CI does not overlap zero

# Livestock activity x water tank interaction
stock.water.dense <- occu( ~ cow.active * water_tank_density_2_4_km 
                           ~ 1, occu.cala)          
null.aicc - AICc(stock.water.dense,k=2)  #worse than null
confint(stock.water.dense, level=0.85, type="det") # 85% CI overlaps zero

# cow active proceeds 

#### Prey Activity Hypothesis Group ####

lagomorph.total <- occu( ~ lagomorph.count 
                         ~ 1, occu.cala, starts = c(-1, 0, 0))
null.aicc - AICc(lagomorph.total,k=2) #worse than null

lagomorph.hours <- occu( ~ lagomorph.active 
                         ~ 1, occu.cala, starts = c(-1, -1, -1))     
null.aicc - AICc(lagomorph.hours,k=2) #better than null
confint(lagomorph.hours, level=0.85, type="det") # 85% CI does not overlap zero

cottontail.total <- occu( ~ cottontail.count 
                          ~ 1, occu.cala, starts = c(-1, 0, 0))
null.aicc - AICc(cottontail.total,k=2) #worse than null

cottontail.hours <- occu( ~ cottontail.active 
                          ~ 1, occu.cala, starts = c(-1, -1, -1))     
null.aicc - AICc(cottontail.hours,k=2) #worse than null

jackrabbit.total <- occu( ~ jackrabbit.count
                          ~ 1, occu.cala, starts = c(-1, 0, 0))
null.aicc - AICc(jackrabbit.total,k=2) #worse than null

jackrabbit.hours <- occu( ~ jackrabbit.active 
                          ~ 1, occu.cala, starts = c(-1, -1, -1))     
null.aicc - AICc(jackrabbit.hours,k=2) #worse than null

ungulate.total <- occu( ~ ungulate.count 
                        ~ 1, occu.cala, starts = c(-1, -1, -1))
null.aicc - AICc(ungulate.total,k=2) #worse than null

ungulate.hours <- occu( ~ ungulate.active
                        ~ 1, occu.cala, starts = c(-1, -1, -1))     
null.aicc - AICc(ungulate.hours,k=2) #better than null
confint(ungulate.hours, level=0.85, type="det") # 85% CI does not overlap zero

# lagomorph hours proceeds

##### Human Activity Hypothesis Group#####

vehicle.total <- occu( ~ vehicle.count ~ 1, occu.cala, starts = c(0, -1, -1))
null.aicc - AICc(vehicle.total,k=2) #worse than null

vehicle.hours <- occu( ~ vehicle.active ~ 1, occu.cala, starts = c(0, -1, -1))     
null.aicc - AICc(vehicle.hours,k=2) #worse than null

people.hours <- occu( ~ people.active ~ 1, occu.cala, starts = c(0, -1, -1))
null.aicc - AICc(people.hours,k=2) #worse than null

human.hours <- occu( ~ human.active ~ 1, occu.cala, starts = c(0, -1, -1)) 
null.aicc - AICc(human.hours,k=2) #better than null
confint(human.hours, level=0.85, type="det") # 85% CI does not overlap zero

# human hours proceeds

#### Road Hypothesis Group #####

small.rd <- occu( ~ small_rd_dist ~ 1, occu.cala, starts = c(-1, -1, -1))     
null.aicc - AICc(small.rd,k=2) #worse than null

##### Biotic Community Type Hypothesis Group ####

wood <- occu( ~ woodland_percent_2_4km ~ 1, occu.cala, starts = c(-1, -1, -1))          
null.aicc - AICc(wood) #worse than null

wood2 <- occu( ~ woodland_percent_2_4km + woodland_percent_2_4km2 
               ~ 1, occu.cala)          
null.aicc - AICc(wood2 ,k=2) #better than null
confint(wood2, level=0.85, type="det") 
# 85% CI does not overlaps zero but is positive (no support for hypothesis)

ndvi <- occu( ~ NDVI_2_4km ~ 1, occu.cala, starts = c(-1, -1, -1))       
null.aicc - AICc(ndvi, k=2) #worse than null

ndvi2 <- occu( ~ NDVI_2_4km + NDVI_2_4km2 
               ~ 1, occu.cala, starts = c(-1, -1, -1, -1))       
null.aicc - AICc(ndvi2, k=2) #better than null
confint(ndvi2, level=0.85, type="det") 
# 85% CI does not overlaps zero but is positive (no support for hypothesis)

bio.com2 <- occu( ~ as.factor(biotic_com_2) ~ 1, occu.cala, starts = c(1, 0, 0))
null.aicc - AICc(bio.com2, k=2) #worse than null

tree <- occu( ~ tree_density_5_70 ~ 1, occu.cala)
null.aicc - AICc(tree)  #better than null
confint(tree, level=0.85, type="det")
# 85% CI does not overlaps zero but is negative (no support for hypothesis)

tree2 <- occu( ~ tree_density_5_70 + tree_density_5_70_2 ~ 1, occu.cala)
null.aicc - AICc(tree2)  #worse than null

# none proceed

##### Vegetation Concealment Cover Hypothesis Group ####

obs.distance <- occu( ~ obs_dist ~ 1, occu.cala, starts = c(2, 1, 1))       
null.aicc - AICc(obs.distance, k=2) #worse than null

obs.distance2 <- occu( ~ obs_dist + obs_dist2 
                       ~ 1, occu.cala, starts = c(2, 1, 1, -3))       
null.aicc - AICc(obs.distance2, k=2) #worse than null

veg.cov.cam <- occu( ~ veg_cover_cam ~ 1, occu.cala, starts = c(2, 1, 1))       
null.aicc - AICc(veg.cov.cam, k=2) #worse than null

veg.cov.cam2 <- occu( ~ veg_cover_cam + veg_cover_cam_2 
                      ~ 1, occu.cala, starts = c(2, 1, 1, -3))       
null.aicc - AICc(veg.cov.cam2, k=2) #worse than null

detect.angle <- occu( ~ detection_angle ~ 1, occu.cala, starts = c(2, 1, 1))  
null.aicc - AICc(detect.angle, k=2) #worse than null

detect.angle2 <- occu( ~ detection_angle + detection_angle2 ~ 1, occu.cala, 
                         starts = c(2, 1, 1, -3))  
null.aicc - AICc(detect.angle2, k=2) #worse than null

max.trig.dist <- occu( ~ max_trig_dist ~ 1, occu.cala, starts = c(2, 1, 1))   
null.aicc - AICc(max.trig.dist, k=2) #worse than null

max.trig.dist2 <- occu( ~ max_trig_dist + max_trig_dist2 
                        ~ 1, occu.cala, starts = c(2, 1, 1, -3))   
null.aicc - AICc(max.trig.dist2, k=2) #worse than null

# none proceed

#### Sampling Period Hypothesis Group ####

### Sampling Period Length ###
samp.period <- occu( ~ sampling_period_length 
                     ~ 1, occu.cala, starts = c(2, 1, -3))          
null.aicc - AICc(samp.period,k=2) #worse than null

### Camera Moved ###
cam.moved <- occu( ~ as.factor(cam_moved) 
                   ~ 1, occu.cala, starts = c(-1, -1, -1))   
null.aicc - AICc(cam.moved, k=2) #worse than null

#none proceed

#################################################
# DETECTION MODEL SELECTION #
#################################################

# Variables that proceed (3)
  # cow active
  # human active
  # lagomorph active 

# Check correlations between detection variables 

# Read in observation-level covariates
obs.covs <- readRDS("obsCovs.RData")

obs_cor <- cor(data.frame(
  cow = as.vector(obs.covs$cow.active),
  human = as.vector(obs.covs$human.active),
  lagomorph = as.vector(obs.covs$lagomorph.active)
))
# none correlated above |0.7|

# Run this code to get correlation matrix in Excel 
  #write.xlsx(obs_cor, file="Correlations coyote.xlsx")

# Null model
mod_null <- occu( ~ 1 ~ 1, occu.cala, starts = c(2, -3))

# single variable models 
mod1 <- occu( ~ human.active ~ 1, occu.cala, starts = c(2, 0, -3))
mod2 <- occu( ~ lagomorph.active ~ 1, occu.cala, starts = c(2, 0, -3))
mod3 <- occu( ~ cow.active ~ 1, occu.cala, starts = c(2, 0, -3))

#2-variable models
mod4 <- occu( ~ human.active + lagomorph.active 
              ~ 1, occu.cala, starts = c(2, 0, -3, 0))
mod5 <- occu( ~ human.active + cow.active 
              ~ 1, occu.cala, starts = c(2, 0, -3, 0))
mod6 <- occu( ~ lagomorph.active + cow.active 
              ~ 1, occu.cala, starts = c(2, 0, -3, 0))

# Model selection
top_mods <- model.sel(mod1, mod2, mod3, mod4, mod5, mod6, mod_null)

# Candidate detection models are found in Table S2.5.

# Run this code to see candidate detection models 
  #write.xlsx(top_mods, file="Coyote Base Models.xlsx", 
          # sheetName="Detection Models", append=T)  

# Top model diagnostics

summary(mod5)
AICc(mod_null) - AICc(mod5) 

### Calculate the 85% confidence intervals for variables 
confint(mod5, type = "det", level = 0.85)    
#CIs for detection variables - none overlap zero

confint(mod5, type = "state", level = 0.85) 
#CIs for occupancy intercept - does not overlap zero

# Calculate multicollinearity
unmarked::vif(mod5, type = "det") # all below 2

# Calculate mean daily detection probability with the top detection model

backTransform(linearComb(mod5,                           
                         coefficients= c(1,0,0),  
                         type = 'det'))     
# 0.0543  or 5.4%

# Occupancy probability

backTransform(linearComb(mod5,                
                         coefficients=c(1),   
                         type = 'state')) 
#  0.763 OR 76.3%

### Calculating overall detection probability if K = 29
1-(1-0.0543)^29

# 80% probability of detecting a mule deer at least once given the minimum
  #sampling period length at a site

### Calculating overall detection probability if K = 38
1-(1-0.0543)^38

# 88% probability of detecting a mule deer at least once given the average
  #sampling period length at a site

#################################################
# CREATE THE HABITAT SELECTION (OCCUPANCY) MODEL#
#################################################

## Mule deer habitat selection hypotheses are listed in Table S2.2.

## Table S2.4 lists which occupancy hypothesis groups proceeded (i.e. performed 
  # better than the null occupancy model ψ(.), did not have uninformative 
  # parameters, and supported our original hypothesis).

## Candidate habitat selection (occupancy) models are found in Table S2.6.
  # Models were limited to 6 parameters and were derived from the best-supported
  # detection model (mod1 above).

occu.null.aicc <- AICc(mod5)

#### Water Hypothesis Group #####

water.tank.dense <- occu(  ~ human.active + cow.active
                           ~ water_tank_density_2_4_km, occu.cala)       
occu.null.aicc - AICc(water.tank.dense)  #worse than null 

water.dist.tank <- occu( ~ human.active + cow.active
                         ~ dist_water_tank, occu.cala)      
occu.null.aicc - AICc(water.dist.tank, k=2) #worse than null 

#none proceed

#### Prey Activity Hypothesis Group####

prey <- occu( ~ human.active + cow.active
              ~ prey_count_avg, occu.cala)
occu.null.aicc - AICc(prey, k=2) #worse than null 

lagomorph <- occu( ~ human.active + cow.active
                   ~ lagomorph_count_avg, occu.cala)
occu.null.aicc - AICc(lagomorph, k=2) #worse than null 

rodent <- occu( ~ human.active + cow.active
                ~ rodent_count_avg, occu.cala)
occu.null.aicc - AICc(rodent, k=2) #worse than null 

jackrabbit <- occu( ~ human.active + cow.active
                    ~ jackrabbit_count_avg, occu.cala)
occu.null.aicc - AICc(jackrabbit, k=2) #worse than null 

cottontail <- occu( ~ human.active + cow.active
                    ~ cottontail_count_avg, occu.cala)
occu.null.aicc - AICc(cottontail, k=2) #worse than null 

# none proceed 

#### Livestock Activity Hypothesis Group####

stock <- occu( ~ human.active + cow.active
               ~ livestock_count_avg, occu.cala)
occu.null.aicc - AICc(stock, k=2) #worse than null 

cow <- occu( ~ human.active + cow.active
             ~ cow_count_avg, occu.cala)
occu.null.aicc - AICc(cow, k=2) #worse than null 

# none proceed 

#### Biotic Community Type Hypothesis Group ####

tree <- occu( ~ human.active + cow.active
              ~ tree_density_5_70, occu.cala)
occu.null.aicc - AICc(tree) #worse than null 

tree2 <- occu( ~ human.active + cow.active
               ~ tree_density_5_70 + tree_density_5_70_2, occu.cala)
occu.null.aicc - AICc(tree2) #worse than null 

shrubland <- occu( ~ human.active + cow.active
                   ~ shrub_yucca_density, occu.cala)
occu.null.aicc - AICc(shrubland) #worse than null

grassland <- occu( ~ human.active + cow.active
                   ~ herbaceous_cov, occu.cala)
occu.null.aicc - AICc(grassland) #worse than nu

grassland2 <- occu( ~ human.active + cow.active
                    ~ herbaceous_cov + herbaceous_cov2, occu.cala)
occu.null.aicc - AICc(grassland2) #worse than null 

wood2 <- occu(~ human.active + cow.active
              ~ woodland_percent_2_4km + woodland_percent_2_4km2 , occu.cala)          
occu.null.aicc - AICc(wood2,k=2)  #worse than null 

wood <- occu(~ human.active + cow.active
             ~ woodland_percent_2_4km, occu.cala)          
occu.null.aicc - AICc(wood,k=2) #worse than null 

ndvi <- occu( ~ human.active + cow.active ~
                NDVI_2_4km, occu.cala, starts = c(-1, -1, -1, -1, -1))          
occu.null.aicc - AICc(ndvi, k=2)  #worse than null 

ndvi2 <- occu( ~ human.active + cow.active ~
                 NDVI_2_4km + NDVI_2_4km2,
                 occu.cala, starts = c(-1, -1, -1, -1, -1, -1))          
occu.null.aicc - AICc(ndvi2, k=2)  #worse than null 

bio.com2 <- occu( ~ human.active + cow.active 
                  ~ as.factor(biotic_com_2), 
                    occu.cala, starts = c(0, 0.5, 1, 0.1, 0.1))        
occu.null.aicc - AICc(bio.com2, k=2) #worse than null 

# none proceeds

#### Vegetation Concealment Cover Hypothesis Group  ####

vert.cover <- occu( ~ human.active + cow.active
                    ~ vertical_cover, occu.cala)
occu.null.aicc - AICc(vert.cover) #worse than null 

vert.cover2 <- occu( ~ human.active + cow.active
                     ~ vertical_cover + vertical_cover2, occu.cala)
occu.null.aicc - AICc(vert.cover2) #worse than null 

veg.cov <- occu( ~ human.active + cow.active
                 ~ veg_cover, occu.cala)
occu.null.aicc - AICc(veg.cov) #worse than null 

veg.cov2 <- occu( ~ human.active + cow.active
                  ~ veg_cover + veg_cover2, occu.cala)
occu.null.aicc - AICc(veg.cov2) #worse than null 

# none proceed 

##### Topography Hypothesis Group ####

topo.pos <- occu( ~ human.active + cow.active ~ topo_pos, occu.cala)
occu.null.aicc - AICc(topo.pos) #worse than null

slope <- occu( ~ human.active + cow.active ~ slope, occu.cala)
occu.null.aicc - AICc(slope)  #worse than null 

elev <- occu( ~ human.active + cow.active ~ elevation, occu.cala)
occu.null.aicc - AICc(elev)  #worse than null 

# none proceeds

##### Coolness of Site Hypothesis Group ####

aspect <- occu( ~ human.active + cow.active ~ aspect, occu.cala)
occu.null.aicc - AICc(aspect) #worse than null 

heat.load <- occu( ~  human.active + cow.active ~ heat_load, occu.cala)
occu.null.aicc - AICc(heat.load) #worse than null 

#none proceeds

#################################################
# OCCUPANCY MODEL SELECTION #
#################################################

# Variables that proceed (0)

# Detection only model 
hab.mod <- occu( ~ human.active + cow.active ~ 1, occu.cala,
                   starts = c(2, 0, -3, 0))

# The parameter estimates, standard errors, and 85% confidence intervals for
  #the best-supported model describing the probability of detection (p) and 
  #habitat selection (ψ) of coyotes are listed in Table S2.9.

# Predict detection probability for survey sites 

det_p_full <- predict(hab.mod, type = "det") 

# Create site index
n_sites <- numSites(occu.cala)
n_occasions <- ncol(getY(occu.cala))
site_index <- rep(1:n_sites, each = n_occasions)

# Add site ID to prediction output
det_p_full$site <- site_index

#  Aggregate by site
det_p_by_site <- det_p_full %>%
  group_by(site) %>%
  summarise(
    mean.p     = mean(Predicted, na.rm = TRUE),
    mean.SE    = sqrt(sum(SE^2, na.rm = TRUE)) / n(),  
    lower85    = mean.p - qnorm(0.925) * mean.SE,
    upper85    = mean.p + qnorm(0.925) * mean.SE
  )

# Read in site-level covariates
site.covs <- read.csv("site_covs.csv", nrows = 102, header = TRUE)

#Predict occupancy probability for survey sites 

pred_occu <- predict(hab.mod ,          
                     type = "state",                 
                     newdata = occu.cala@siteCovs)[c("Predicted",
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