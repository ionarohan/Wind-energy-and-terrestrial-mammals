###########################################################
## Code for: Turbine visibility is a strong predictor of ##
## altered habitat selection by terrestrial mammals at a ##
######## wind energy facility in central New Mexico #######
###########################################################
## This script is for the creation of the striped skunk ###
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
# SETUP CODE FOR STRIPED SKUNK OCCUPANCY MODELS #
###########################################################

# Read in edited .csv
detHist <- read.csv(file = "skunk detection hist.csv", row.names = 1)
# Read in site.covs.scaled
site.covs.scaled <- readRDS("site_covs_scaled.RData")
# Read in obsCovs.scaled
obsCovs.scaled <- readRDS("obsCovs_scaled.RData")

# Change from integer to numeric
detHist %>%
  mutate(across(1:51, as.numeric))

# Create detection history 
occu.meme <- unmarkedFrameOccu(y=detHist, 
                               siteCovs = site.covs.scaled, 
                               obsCovs = obsCovs.scaled)

#### Null model #####
meme.null <- occu( ~ 1 ~ 1, occu.meme, linkPsi="logit", se=TRUE,
                  control = list(maxit = 10000))
null.aicc <- AICc(meme.null)

# Null detection probability (Table 2)
backTransform(meme.null['det'])

### Calculating overall detection probability if K = 38
1-(1-0.0801)^38

# 96% probability of detecting a striped skunk at least once given the average
#sampling period length at a site and the null detection model (see Table 2)

# Null occupancy probability 
backTransform(meme.null['state'])

#################################################
# CREATE THE DETECTION MODEL #
#################################################

## Striped skunk detection hypotheses are listed in Table S2.1.

## Table S2.3 lists which detection hypothesis groups proceeded (i.e. performed 
  # better than the null detection model p(.), did not have uninformative 
  # parameters, and supported our original hypothesis).

## Candidate detection models are found in Table S2.5.
  # Models were limited to 4 parameters total and were tested with no variables 
  # on the probability of habitat selection (ψ(.)).

##### Water Hypothesis Group ####

water.tank.dense <- occu( ~ water_tank_density_1_3km 
                    ~ 1, occu.meme, starts = c(-1, 0, 0))       
null.aicc - AICc(water.tank.dense) #better than null
confint(water.tank.dense, level=0.85, type="det")  #85% CI does not overlap zero

water.dist.tank <- occu( ~ dist_water_tank 
                         ~ 1, occu.meme, starts = c(-1, 0, 0))  
null.aicc - AICc(water.dist.tank, k=2)  # worse than null 

# water tank density proceeds

##### Precipitation Hypothesis Group####

# Daily precipitation 

precip.cat <- occu( ~ as.factor(precip.cat) 
                    ~ 1, occu.meme, starts = c(-1, 0, 0))
null.aicc - AICc(precip.cat) # worse than null 

precip.cm <- occu( ~ precip.cm ~ 1, occu.meme, starts = c(-1, 0, 0))
null.aicc - AICc(precip.cm) # worse than null 

# Long-term precipitation

days.since.rain <- occu( ~ days.since.rain ~ 1 , occu.meme)    
null.aicc - AICc(days.since.rain, k=2) # worse than null 

rain.month <- occu( ~ rain.month ~ 1, occu.meme)    
null.aicc - AICc(rain.month, k=2) # worse than null 

rain.week <- occu( ~ rain.week ~ 1, occu.meme)    
null.aicc - AICc(rain.week, k=2) # worse than null 

# Water source and long-term precipitation interaction

precip.water.dense <-  occu( ~ water_tank_density_1_3km * 
                               rain.month ~ 1, occu.meme)   
null.aicc - AICc(precip.water.dense, k=2) #better than null
confint(precip.water.dense, level=0.85, type="det") 
#85% CI does not overlap zero

# Interaction lower AICc than water tank density, does not proceed

#### Humidity Hypothesis Group ####

humidity <- occu( ~ humid ~ 1, occu.meme)       
null.aicc - AICc(humidity, k=2) #better than null
confint(humidity, level=0.85, type="det")
# 85% CI does not overlaps zero but is negative (no support for hypothesis)

# none proceed

##### Temperature Hypothesis Group ####

temp.max <- occu( ~ max.temp ~ 1, occu.meme, starts = c(-1, 0, 0))       
null.aicc - AICc(temp.max, k=2) #better than null
confint(temp.max, level=0.85, type="det")
# 85% CI does not overlaps zero but is positive (no support for hypothesis)

# Temperature interactions

temp.water.dense <- occu( ~ max.temp * water_tank_density_1_3km
                          ~ 1, occu.meme)       
null.aicc - AICc(temp.water.dense, k=2) #better than null
confint(temp.water.dense, level=0.85, type="det") #85% CI overlaps zero

temp.canopy <- occu( ~ max.temp * canopy_cov ~ 1, occu.meme)      
null.aicc - AICc(temp.canopy, k=2) #better than null
confint(temp.canopy, level=0.85, type="det") 
# 85% CI does not overlaps zero but is negative (no support for hypothesis)

temp.humid <- occu( ~ max.temp * humid ~ 1, occu.meme)      
null.aicc - AICc(temp.humid, k=2) #worse than null

#none proceed

#### Wind Hypothesis Group ####

wind <- occu( ~ wind ~ 1, occu.meme, starts = c(-1, 0, 0))       
null.aicc - AICc(wind, k=2) #worse than null

# none proceed

##### Human Activity Hypothesis Group#####

vehicle.total <- occu( ~ vehicle.count ~ 1, occu.meme, starts = c(-2, -2, -5))
#large SE for vehicle count
null.aicc - AICc(vehicle.total,k=2) #worse than null

vehicle.hours <- occu( ~ vehicle.active ~ 1, occu.meme, starts = c(-2, -2, -1))
#large SE for vehicle hours
null.aicc - AICc(vehicle.hours,k=2) #worse than null

people.hours <- occu( ~ people.active ~ 1, occu.meme, starts = c(-1, 0, 0))
null.aicc - AICc(people.hours,k=2) #worse than null

human.hours <- occu( ~ human.active ~ 1, occu.meme, starts = c(-1, 0, 0)) 
null.aicc - AICc(human.hours,k=2) #worse than null

# none proceeds

##### Biotic Community Type Hypothesis Group ####

ndvi <- occu( ~ NDVI_1_3km ~ 1, occu.meme, starts = c(-1, 0, 0))          
null.aicc - AICc(ndvi, k=2) #worse than null 

wood <- occu( ~ woodland_percent_1_3km ~ 1, occu.meme)          
null.aicc - AICc(wood, k=2) #better than null
confint(wood, level=0.85, type="det") #85% CI does not overlap zero

bio.com <- occu( ~ as.factor(biotic_com_2) ~ 1, occu.meme, 
                 starts = c(-2, -2, -1))        
null.aicc - AICc(bio.com, k=2)  #better than null
confint(bio.com, level=0.85, type="det") #85% CI does not overlap zero

tree <- occu( ~  tree_density_5_70 ~ 1, occu.meme)
null.aicc - AICc(tree)  #better than null
confint(tree, level=0.85, type="det") #85% CI overlaps zero

# woodland percent proceeds

##### Vegetation Concealment Cover Hypothesis Group ####

obs.distance <- occu( ~ obs_dist ~ 1, occu.meme)       
null.aicc - AICc(obs.distance, k=2) #better than null
confint(obs.distance, level=0.85, type="det") # 85% CI includes zero

veg.cov.cam <- occu( ~ veg_cover_cam ~ 1, occu.meme)       
null.aicc - AICc(veg.cov.cam, k=2) #better than null
confint(veg.cov.cam, level=0.85, type="det") # 85% CI does not include zero

detect.angle <- occu( ~ detection_angle ~ 1, occu.meme, starts = c(-1, 0, 0))  
null.aicc - AICc(detect.angle, k=2) #worse than null

max.trig.dist <- occu( ~ max_trig_dist ~ 1, occu.meme, starts = c(-1, 0, 0))   
null.aicc - AICc(max.trig.dist, k=2) #better than null
confint(max.trig.dist, level=0.85, type="det") # 85% CI does not include zero

# veg cover cam proceeds

#### Sampling Period Hypothesis Group ####

### Sampling Period Length ###
samp.period <- occu( ~ sampling_period_length ~ 1, occu.meme)          
null.aicc - AICc(samp.period,k=2) #better than null
confint(samp.period, level=0.85, type="det") # 85% CI does not include zero

### Camera Moved ###
cam.moved <- occu( ~ as.factor(cam_moved) ~ 1, occu.meme,starts = c(-1, 0, 0))  
null.aicc - AICc(cam.moved, k=2) #worse than null

# sampling period proceeds

#################################################
# DETECTION MODEL SELECTION #
#################################################

# Variables that proceed (4)
  #sampling period
  #water tank
  #veg cover cam
  #woodland percent 

##### check correlations between top variables ####

# Read in site-level covariates
site.covs <- read.csv("site_covs.csv", nrows = 102, header = TRUE)

site.covs.cor <- site.covs %>% 
  select(sampling_period_length, 
         water_tank_density_1_3km, veg_cover_cam,
         woodland_percent_1_3km) 

cor_site_covs <- cor(site.covs.cor, method='spearman') 
# none correlated above |0.7|

# Null model
meme.null <- occu( ~ 1 ~ 1, occu.meme)

# 1-variable models
mod1 <- occu( ~ sampling_period_length ~ 1, occu.meme)
mod2 <- occu( ~ veg_cover_cam ~ 1, occu.meme)
mod3 <- occu( ~ water_tank_density_1_3km ~ 1, occu.meme)
mod4 <- occu( ~ woodland_percent_1_3km ~ 1, occu.meme)

# 2-variable models
mod5 <- occu( ~ sampling_period_length + veg_cover_cam 
              ~ 1, occu.meme)
mod6 <- occu( ~ sampling_period_length + water_tank_density_1_3km 
              ~ 1, occu.meme)
mod7 <- occu( ~ sampling_period_length + woodland_percent_1_3km
              ~ 1, occu.meme)
mod8 <- occu( ~ veg_cover_cam + water_tank_density_1_3km 
              ~ 1, occu.meme)
mod9 <- occu( ~ veg_cover_cam + woodland_percent_1_3km
              ~ 1, occu.meme)
mod10 <- occu( ~ water_tank_density_1_3km + woodland_percent_1_3km
              ~ 1, occu.meme)

# Model selection
cand.models <- list(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9, 
                    mod10, meme.null)

modnames <- c("mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7", "mod8", 
              "mod9", "mod10", "meme.null")

aicc_table <- aictab(cand.set = cand.models, modnames = modnames, sort = TRUE)
print(aicc_table)
# Candidate detection models are found in Table S2.5.

# Top model diagnostics

summary(mod9)
AICc(meme.null) - AICc(mod9)

# Calculate the 85% confidence intervals for variables 
confint(mod9, type = "det", level = 0.85)    
#CIs for detection variables - none overlap zero

confint(mod9, type = "state", level = 0.85)  
#CIs for occupancy intercept - does not overlap zero

# Calculate multicollinearity
unmarked::vif(mod9, type = "det") # all below 2

#################################################
# CREATE THE HABITAT SELECTION (OCCUPANCY) MODEL#
#################################################

## Striped skunk habitat selection hypotheses are listed in Table S2.2.

## Table S2.4 lists which occupancy hypothesis groups proceeded (i.e. performed 
  # better than the null occupancy model ψ(.), did not have uninformative 
  # parameters, and supported our original hypothesis).

## Candidate habitat selection (occupancy) models are found in Table S2.6.
  # Models were limited to 6 parameters and were derived from the best-supported
  # detection model (mod1 above).

occu.null.aicc <- AICc(mod9)

####Water Hypothesis Group#####

water.tank.dense <- occu( ~ veg_cover_cam + woodland_percent_1_3km
                          ~ water_tank_density_1_3km, occu.meme)       
occu.null.aicc - AICc(water.tank.dense) #worse than null

water.dist.tank <- occu( ~ veg_cover_cam + woodland_percent_1_3km
                         ~ dist_water_tank, occu.meme)      
occu.null.aicc - AICc(water.dist.tank, k=2) #worse than null

#none proceed 

#### Biotic Community Type Hypothesis Group ####

tree <- occu( ~  veg_cover_cam + woodland_percent_1_3km
              ~  tree_density_5_70, occu.meme)
occu.null.aicc - AICc(tree) #worse than null 

tree2 <- occu( ~  veg_cover_cam + woodland_percent_1_3km
               ~  tree_density_5_70 + tree_density_5_70_2, occu.meme)
occu.null.aicc - AICc(tree2) #worse than null 

ndvi <- occu( ~ veg_cover_cam + woodland_percent_1_3km
              ~ NDVI_1_3km, occu.meme)          
occu.null.aicc - AICc(ndvi, k=2) #worse than null 

ndvi2 <- occu( ~ veg_cover_cam + woodland_percent_1_3km
               ~ NDVI_1_3km + NDVI_1_3km2, occu.meme)          
occu.null.aicc - AICc(ndvi2) #better than null 
confint(ndvi2, level=0.85, type="state") 
# 85% CI does not overlaps zero but is positive (no support for hypothesis)

shrubland <- occu( ~ veg_cover_cam + woodland_percent_1_3km
                   ~ shrub_yucca_density, occu.meme)
occu.null.aicc - AICc(shrubland) #worse than null 

shrubland2 <- occu( ~ veg_cover_cam + woodland_percent_1_3km
                    ~ shrub_yucca_density2, occu.meme)
occu.null.aicc - AICc(shrubland2) #worse than null 

grassland <- occu( ~ veg_cover_cam + woodland_percent_1_3km
                   ~ herbaceous_cov, occu.meme)
occu.null.aicc - AICc(grassland) #better than null
confint(grassland, level=0.85, type="state") #85% CI does not overlap zero

grassland2 <- occu( ~ veg_cover_cam + woodland_percent_1_3km
                   ~ herbaceous_cov2, occu.meme)
occu.null.aicc - AICc(grassland2) #worse than null

wood <- occu(  ~ veg_cover_cam + woodland_percent_1_3km
               ~ woodland_percent_1_3km, occu.meme)          
occu.null.aicc - AICc(wood, k=2) #worse than null

wood2 <- occu(  ~ veg_cover_cam + woodland_percent_1_3km
                ~ woodland_percent_1_3km + woodland_percent_1_3km, occu.meme)          
occu.null.aicc - AICc(wood2, k=2) #worse than null

bio.com <- occu( ~ veg_cover_cam + woodland_percent_1_3km
                 ~ as.factor(biotic_com_2), occu.meme, 
                  starts = c(-2, -2, -1, 0, 1))        
occu.null.aicc - AICc(bio.com, k=2) #worse than null

# herbaceous cov proceeds

#### Vegetation Concealment Cover Hypothesis Group  ####

vert.cover <- occu( ~  veg_cover_cam + woodland_percent_1_3km
                    ~  vertical_cover, occu.meme)
occu.null.aicc - AICc(vert.cover) #worse than null

veg.cov <- occu( ~  veg_cover_cam + woodland_percent_1_3km
                 ~  veg_cover, occu.meme)
occu.null.aicc - AICc(veg.cov ) #worse than null

veg.cov2 <- occu( ~  veg_cover_cam + woodland_percent_1_3km
                  ~  veg_cover + veg_cover2, occu.meme)
occu.null.aicc - AICc(veg.cov2) #worse than null

# none proceed 

##### Topography Hypothesis Group ####

topo.pos <- occu( ~  veg_cover_cam + woodland_percent_1_3km
                  ~  topo_pos, occu.meme)
occu.null.aicc - AICc(topo.pos) #better than null
confint(topo.pos, level=0.85, type="state") #85% CI overlaps zero

slope <- occu( ~ veg_cover_cam + woodland_percent_1_3km
               ~  slope, occu.meme)
occu.null.aicc - AICc(slope) #worse than null

elev <- occu( ~ veg_cover_cam + woodland_percent_1_3km
              ~  elevation, occu.meme)
occu.null.aicc - AICc(elev) #worse than null

#none proceed

##### Coolness of Site Hypothesis Group ####

aspect <- occu( ~  veg_cover_cam + woodland_percent_1_3km
                ~  aspect, occu.meme)
occu.null.aicc - AICc(aspect) #better than null
confint(aspect, level=0.85, type="state") 
# 85% CI does not overlaps zero but is positive (no support for hypothesis)

heat.load <- occu( ~ veg_cover_cam + woodland_percent_1_3km
                   ~  heat_load, occu.meme)
occu.null.aicc - AICc(heat.load) #worse than null

#none proceed

#################################################
# OCCUPANCY MODEL SELECTION #
#################################################

# Variables that proceed (1)
  # herbaceous cover

# Detection only model 
null <-  occu( ~ veg_cover_cam + woodland_percent_1_3km 
               ~ 1, occu.meme)

# 1-variable occupancy models
mod1 <-  occu( ~ veg_cover_cam + woodland_percent_1_3km 
               ~ herbaceous_cov, occu.meme)

# Model selection
cand.models <- list(mod1, null)

modnames <- c("mod1", "null")

aicc_table <- aictab(cand.set = cand.models, modnames = modnames, sort = TRUE)
print(aicc_table)

# Candidate occupancy models are found in Table S2.6.

# Top model diagnostics

summary(mod1)
AICc(null)-AICc(mod1)

# The parameter estimates, standard errors, and 85% confidence intervals for
  #the best-supported model describing the probability of detection (p) and 
  #habitat selection (ψ) of striped skunks are listed in Table S2.13.

# Calculate 85% confidence intervals for variables 
confint(mod1, type = "det", level = 0.85)    
#CIs for detection variables do not overlap zero

confint(mod1, type = "state", level = 0.85)  
#CIs for occupancy variables do not overlap zero

# Calculate mean daily detection probability with the top detection model

backTransform(linearComb(mod1,                           
                         coefficients= c(1,0,0),  
                         type = 'det'))     
#  0.0147  or 1.47%

# Occupancy probability

backTransform(linearComb(mod1,                
                         coefficients=c(1,0),   
                         type = 'state')) 
# 0.2  OR 20%

# Calculate the overall detection probability if K = 29
1-(1-0.0147)^29

# 35% probability of detecting a striped skunk at least once given the minimum
  #sampling period length at a site

# Calculate the overall detection probability if K = 38 
1-(1-0.0147)^38

# 43% probability of detecting a striped skunk at least once given the average
  #sampling period length at a site

#Predict detection probability for survey sites 

pred_det <- predict(mod1,          
                    type = "det",                 
                    newdata = occu.meme@siteCovs)[c("Predicted",
                                                    "SE",
                                                    "lower",
                                                    "upper")]

pred_det_df <- data.frame(Predicted = pred_det$Predicted,
                          StandardError = pred_det$SE,
                          lower = pred_det$lower,
                          upper = pred_det$upper,
                          site.covs)

#Predict occupancy probability for survey sites 

pred_occu <- predict(mod1,          
                     type = "state",                 
                     newdata = occu.meme@siteCovs)[c("Predicted",
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