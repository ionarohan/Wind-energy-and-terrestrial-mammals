###########################################################
########### Code for: base models for pronghorns ##########
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
# SETUP CODE FOR PRONGHORN OCCUPANCY MODELS #
###########################################################

# Read in edited .csv
detHist <- read.csv(file = "pronghorn detection hist.csv", row.names = 1)

# Change from integer to numeric
detHist %>%
  mutate(across(1:51, as.numeric))

# Create detection history 
occu.anam <- unmarkedFrameOccu(y=detHist, 
                               siteCovs = site.covs.scaled, 
                               obsCovs = obsCovs.scaled)

#### Null model ####
anam.null <- occu( ~ 1 ~ 1, occu.anam, linkPsi="logit", se=TRUE,
                   control = list(maxit = 10000))
null.aicc <- AICc(anam.null)

# Null detection probability 
backTransform(anam.null['det']) 

# Null occupancy probability 
backTransform(anam.null['state'])

#################################################
# CREATE THE DETECTION MODEL #
#################################################

## Pronghorn detection hypotheses are listed in Table S2.1.

## Table S2.3 lists which detection hypothesis groups proceeded (i.e. performed 
  # better than the null detection model p(.), did not have uninformative 
  # parameters, and supported our original hypothesis).

## Candidate detection models are found in Table S2.5.
  # Models were limited to 4 parameters total and were tested with no variables 
  # on the probability of habitat selection (ψ(.)).

##### Water Hypothesis Group ####

water.tank.dense <- occu( ~ water_tank_density_3_7_km
                          ~ 1, occu.anam, starts = c(-1, 0, 0))       
null.aicc - AICc(water.tank.dense) #worse than null 

water.dist.tank <- occu( ~ dist_water_tank 
                         ~ 1, occu.anam, starts = c(-1, 0, 0))       
null.aicc - AICc(water.dist.tank, k=2) #worse than null 

#none proceed

##### Precipitation Hypothesis Group####

# Daily precipitation 

precip.cat <- occu( ~ as.factor(precip.cat) 
                    ~ 1, occu.anam, starts = c(-1, 0, 0))
null.aicc - AICc(precip.cat) #worse than null 

precip.cm <- occu( ~ precip.cm ~ 1, occu.anam, starts = c(-1, 0, 0))
null.aicc - AICc(precip.cm) #worse than null 

# Long-term precipitation

day.since.rain <- occu( ~ days.since.rain 
                        ~ 1, occu.anam, starts = c(-1, 0, 0))    
null.aicc - AICc(day.since.rain, k=2) #better than null
confint(day.since.rain, level=0.85, type="det") #85% CI overlaps zero

rain.month <- occu( ~ rain.month ~ 1, occu.anam, starts = c(-1, 0, 0))    
null.aicc - AICc(rain.month, k=2) #worse than null 

rain.week <- occu( ~ rain.week ~ 1, occu.anam, starts = c(-1, 0, 0))    
null.aicc - AICc(rain.week, k=2) #worse than null 

# Water source and long-term precipitation interaction

precip.water.dense <-  occu( ~ water_tank_density_3_7_km * 
                               days.since.rain ~ 1, occu.anam)   
null.aicc - AICc(precip.water.dense, k=2) #worse than null 

#none proceed

##### Temperature Hypothesis Group ####

temp.max <- occu( ~ max.temp ~ 1, occu.anam, starts = c(-1, 0, 0))
null.aicc - AICc(temp.max, k=2) #worse than null 

# Temperature interactions

temp.water.dense.season <- occu( ~ max.temp * water_tank_density_3_7_km
                                 ~ 1, occu.anam)       
null.aicc - AICc(temp.water.dense.season, k=2) #worse than null 

temp.humid <- occu( ~ max.temp * humid ~ 1, occu.anam)      
null.aicc - AICc(temp.humid, k=2) #worse than null 

# none proceed 

#### Wind Hypothesis Group ####

wind <- occu( ~ wind ~ 1, occu.anam)       
null.aicc - AICc(wind , k=2) #worse than null

#none proceed

##### Predator Activity Hypothesis Group ####

coy.total <- occu( ~ coyote.count ~ 1, occu.anam, starts = c(-1, 0, 0))        
null.aicc - AICc(coy.total,k=2) #worse than null 

coy.hours <- occu( ~ coyote.active ~ 1, occu.anam, starts = c(-1, 0, 0))  
null.aicc - AICc(coy.hours,k=2) #better than null
confint(coy.hours, level=0.85, type="det") 
# 85% CI does not overlap zero but is positive (no support for hypothesis)

bob.hours <- occu( ~ bobcat.active ~ 1, occu.anam, starts = c(-1, 0, 0))
null.aicc - AICc(bob.hours,k=2) #worse than null 

bob.total <- occu( ~ bobcat.count ~ 1, occu.anam, starts = c(-1, 0, 0))
null.aicc - AICc(bob.total,k=2) #worse than null 

coy.bob.hours <- occu( ~ coy.bob.active ~ 1, occu.anam, starts = c(-1, 0, 0))
null.aicc - AICc(coy.bob.hours, k=2) #better than null
confint(coy.bob.hours, level=0.85, type="det") 
# 85% CI does not overlap zero but is positive (no support for hypothesis)

coy.bob.total <- occu( ~ coy.bob.count ~ 1, occu.anam, starts = c(-1, 0, 0)) 
null.aicc - AICc(coy.bob.total,k=2) #worse than null 

# none proceed

##### Human Activity Hypothesis Group#####

vehicle.total <- occu( ~ vehicle.count ~ 1, occu.anam, starts = c(-1, 0, 0))
null.aicc - AICc(vehicle.total,k=2) #worse than null

vehicle.hours <- occu( ~ vehicle.active ~ 1, occu.anam, starts = c(-1, -1, -1))     
null.aicc - AICc(vehicle.hours,k=2) #worse than null

people.hours <- occu( ~ people.active ~ 1, occu.anam, starts = c(-1, 0, 0))
null.aicc - AICc(people.hours,k=2) #better than null
confint(people.hours, level=0.85, type="det") # 85% CI overlaps zero

human.hours <- occu( ~ human.active ~ 1, occu.anam, starts = c(-1, 0, 0)) 
null.aicc - AICc(human.hours,k=2) #better than null
confint(human.hours, level=0.85, type="det") # 85% CI overlaps zero

# none proceed

##### Biotic Community Type Hypothesis Group ####

wood <- occu( ~ woodland_percent_3_7km ~ 1, occu.anam)          
null.aicc - AICc(wood, k=2)  #worse than null

ndvi <- occu( ~ NDVI_3_7km ~ 1, occu.anam)          
null.aicc - AICc(ndvi, k=2)  #worse than null

bio.com <- occu( ~ as.factor(biotic_com_2) ~ 1, occu.anam, starts = c(0, 1, 0))        
null.aicc - AICc(bio.com, k=2)  #better than null
confint(bio.com, level=0.85, type="det") # 85% CI overlaps zero

tree <- occu( ~ tree_density_5_70 ~ 1, occu.anam)
null.aicc - AICc(tree)  #better than null
confint(tree, level=0.85, type="det") # 85% CI does not overlap zero

#tree density proceeds

##### Vegetation Concealment Cover Hypothesis Group ####

obs.distance <- occu( ~ obs_dist ~ 1, occu.anam)       
null.aicc - AICc(obs.distance, k=2)  #worse than null

veg.cov.cam <- occu( ~ veg_cover_cam ~ 1, occu.anam)    
null.aicc - AICc(veg.cov.cam, k=2) #better than null
confint(veg.cov.cam, level=0.85, type="det") # 85% CI does not overlap zero

detect.angle <- occu( ~ detection_angle ~ 1, occu.anam, starts = c(-1, 0, 0))  
null.aicc - AICc(detect.angle, k=2) #better than null
confint(detect.angle, level=0.85, type="det") # 85% CI does not overlap zero

max.trig.dist <- occu( ~ max_trig_dist ~ 1, occu.anam, starts = c(-1, 0, 0))   
null.aicc - AICc(max.trig.dist, k=2) #worse than null

# veg cover proceeds

#### Sampling Period Hypothesis Group ####

### Sampling Period Length ###
samp.period <- occu( ~ sampling_period_length ~ 1, occu.anam)          
null.aicc - AICc(samp.period,k=2) #worse than null

### Camera Moved ###
cam.moved <- occu( ~ as.factor(cam_moved) ~ 1, occu.anam, starts = c(-1, 0, 0))
null.aicc - AICc(cam.moved, k=2) #worse than null

# none proceed

#################################################
# DETECTION MODEL SELECTION #
#################################################

# Variables that proceed (2)
  # tree density
  # veg cover

site.covs.cor <- site.covs %>% 
  select(tree_density_5_70, veg_cover_cam)

cor_site_covs <- cor(site.covs.cor, method='spearman') 
#correlated = 0.87 (use visual obstruction)

# Null model
mod_null <- occu( ~ 1 ~ 1, occu.anam)

# Single variable model
mod1 <- occu( ~ veg_cover_cam ~ 1, occu.anam)

# Model selection
top_mods <- model.sel(mod1, mod_null)

# Candidate detection models are found in Table S2.5.

# Run this code to see candidate detection models 
#write.xlsx(top_mods, file="Pronghorn Base Models.xlsx", 
         #  sheetName="Detection Models", append=T)  

# Top model diagnostics

summary(mod1)
AICc(mod_null)-AICc(mod1)

# Calculate 85% confidence intervals for variables 
confint(mod1, type = "det", level = 0.85)    
#CIs for detection variables - none overlap zero

confint(mod1, type = "state", level = 0.85)  
#CIs for detection variables - none overlap zero

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

occu.null.aicc <- AICc(mod1)

#### Water Hypothesis Group ####

water.tank.dense <- occu( ~ veg_cover_cam 
                          ~ water_tank_density_3_7_km, occu.anam)       
occu.null.aicc - AICc(water.tank.dense) #worse than null

water.dist.tank <- occu( ~ veg_cover_cam 
                         ~ dist_water_tank, occu.anam)       
occu.null.aicc - AICc(water.dist.tank, k=2) #worse than null

#none proceed 

#### Predator Activity Hypothesis Group####

bobcat.coyote <- occu( ~ veg_cover_cam ~ pred_count_avg, occu.anam)
occu.null.aicc - AICc(bobcat.coyote, k=2) #worse than null 

coyote <- occu( ~ veg_cover_cam ~ coy_count_avg, occu.anam)
occu.null.aicc - AICc(coyote, k=2) #worse than null 

# none proceed 

#### Biotic Community Type Hypothesis Group ####

woodland <- occu( ~ veg_cover_cam ~ tree_density_5_70, occu.anam)
occu.null.aicc - AICc(woodland) #worse than null 

shrubland <- occu( ~ veg_cover_cam ~ shrub_yucca_density, occu.anam)
occu.null.aicc - AICc(shrubland) #worse than null 

grassland <- occu( ~ veg_cover_cam ~ herbaceous_cov, occu.anam)
occu.null.aicc - AICc(grassland) #worse than null 

wood <- occu(  ~ veg_cover_cam ~ woodland_percent_3_7km, occu.anam)          
occu.null.aicc - AICc(wood, k=2)  #worse than null

ndvi <- occu(  ~ veg_cover_cam ~ NDVI_3_7km, occu.anam)          
occu.null.aicc - AICc(ndvi, k=2)  #worse than null

bio.com <- occu( ~ veg_cover_cam
                 ~ as.factor(biotic_com_2), occu.anam)          
occu.null.aicc - AICc(bio.com, k=2)  #better than null
confint(bio.com, type = "state", level = 0.85) #85% CI overlaps zero

# none proceed 

#### Vegetation Concealment Cover Hypothesis Group  ####

vert.cover <- occu( ~ veg_cover_cam ~ vertical_cover, occu.anam)
occu.null.aicc - AICc(vert.cover) #worse than null 

veg.cov <- occu( ~ veg_cover_cam ~ veg_cover, occu.anam)
occu.null.aicc - AICc(veg.cov) #worse than null 

# none proceed

##### Topography Hypothesis Group ####

topo.pos <- occu( ~ veg_cover_cam ~ topo_pos, occu.anam)
occu.null.aicc - AICc(topo.pos) #better than null
confint(topo.pos, type = "state", level = 0.85) 
# 85% CI does not overlap zero but is positive (no support for hypothesis)

slope <- occu( ~ veg_cover_cam ~ slope, occu.anam)
occu.null.aicc - AICc(slope) #better than null
confint(slope, type = "state", level = 0.85) #85% CI does not overlap zero

elev <- occu( ~ veg_cover_cam ~ elevation, occu.anam)
occu.null.aicc - AICc(elev) #better than null
confint(elev, type = "state", level = 0.85) #85% CI does not overlap zero

# slope proceeds

##### Coolness of Site Hypothesis Group ####

aspect <- occu( ~ veg_cover_cam ~ aspect, occu.anam)
occu.null.aicc - AICc(aspect) #worse than null 

heat.load <- occu( ~ veg_cover_cam ~  heat_load, occu.anam)
occu.null.aicc - AICc(heat.load) #worse than null 

#none proceed

#################################################
# OCCUPANCY MODEL SELECTION #
#################################################

# Variables that proceed (1)
  # slope

# Detection only model 
det.mod <- occu( ~ veg_cover_cam ~ 1, occu.anam)

# 1-variable occupancy models
mod1 <- occu( ~ veg_cover_cam ~  slope, occu.anam)

# Combine models in model selection table
top_mods <- model.sel(det.mod, mod1)

# Candidate occupancy models are found in Table S2.6.

# Run this code to see candidate occupancy models 
#write.xlsx(top_mods, file="Pronghorn Base Models.xlsx", 
        #   sheetName="Occupancy Models", append=T)  

# Top model diagnostics

summary(mod1)
AICc(det.mod) - AICc(mod1)

# The parameter estimates, standard errors, and 85% confidence intervals for
  #the best-supported model describing the probability of detection (p) and 
  #habitat selection (ψ) of pronghorns are listed in Table S2.8.

### Calculating 85% confidence intervals for variables 
confint(mod1, type = "det", level = 0.85)    
#CIs for detection variables do not overlap zero

confint(mod1, type = "state", level = 0.85)  
#CIs for occupancy variables do not overlap zero

# Calculate mean daily detection probability with the top detection model

backTransform(linearComb(mod1,                           
                         coefficients= c(1,0),  
                         type = 'det'))     
#0.017 or 1.7%

# Calculate the overall detection probability if K = 29
1-(1-0.017)^29

# 39% probability of detecting a pronghorn at least once given the minimum
  #sampling period length at a site

# Calculate the overall detection probability if K = 38 
1-(1-0.017)^38

# 48% probability of detecting a pronghorn at least once given the average
  #sampling period length at a site

# Calculate the occupancy probability

backTransform(linearComb(mod1,                
                         coefficients=c(1,0),   
                         type = 'state')) 
#0.122 OR 12%

#Predict detection probability for survey sites 

pred_prong <- predict(mod1,          
                      type = "det",                 
                      newdata = occu.anam@siteCovs)[c("Predicted",
                                                       "SE",
                                                       "lower",
                                                       "upper")]

pred_prong_df <- data.frame(Predicted = pred_prong$Predicted,
                                     StandardError = pred_prong$SE,
                                     lower = pred_prong$lower,
                                     upper = pred_prong$upper,
                                     site.covs)

# Predict occupancy probability for survey sites 

pred_occu <- predict(mod1,          
                      type = "state",                 
                      newdata = occu.anam@siteCovs)[c("Predicted",
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