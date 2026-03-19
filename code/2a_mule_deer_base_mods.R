###########################################################
########### Code for: base models for mule deer ###########
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
# SETUP CODE FOR MULE DEER OCCUPANCY MODELS #
###########################################################

# Read in edited .csv
detHist <- read.csv(file = "mule deer detection hist.csv", row.names = 1)

# Change from integer to numeric
detHist %>%
  mutate(across(1:51, as.numeric))

# Create detection history 
occu.odhe <- unmarkedFrameOccu(y=detHist, 
                               siteCovs = site.covs.scaled, 
                               obsCovs = obsCovs.scaled)

#### Null model ####
odhe.null <- occu( ~ 1 ~ 1, occu.odhe, linkPsi="logit", starts = c(2, -3), 
                   se = TRUE, control = list(maxit = 10000))

null.aicc <- AICc(odhe.null)

# Null detection probability 
backTransform(odhe.null['det']) 

# Null occupancy probability 
backTransform(odhe.null['state']) 

#################################################
# CREATE THE DETECTION MODEL #
#################################################

## Mule deer detection hypotheses are listed in Table S2.1.

## Table S2.3 lists which detection hypothesis groups proceeded (i.e. performed 
  # better than the null detection model p(.), did not have uninformative 
  # parameters, and supported our original hypothesis).

## Candidate detection models are found in Table S2.5.
  # Models were limited to 4 parameters total and were tested with no variables 
  # on the probability of habitat selection (ψ(.)).

##### Water Hypothesis Group ####

water.tank.dist <- occu( ~ dist_water_tank ~ 1, occu.odhe)
null.aicc - AICc(water.tank.dist, k=2) #worse than null 

water.tank.dense <- occu( ~ water_tank_density_1_9_km ~ 1, occu.odhe)       
null.aicc - AICc(water.tank.dense) #worse than null 

#none proceed

##### Precipitation Hypothesis Group####

# Daily precipitation 
precip.cat <- occu( ~ as.factor(precip.cat) ~ 1, occu.odhe, 
                      starts = c(2, 0, -3))
null.aicc - AICc(precip.cat)  #worse than null 

precip.cm <- occu( ~ precip.cm ~ 1, occu.odhe, starts = c(2, 0, -3))
null.aicc - AICc(precip.cm)  #worse than null 

# Long-term precipitation

days.since.rain <- occu( ~ days.since.rain ~ 1, occu.odhe, starts = c(2, 0, -3))    
null.aicc - AICc(days.since.rain, k=2)  #worse than null 

rain.week <- occu( ~ rain.week ~ 1, occu.odhe, starts = c(2, 0, -3))    
null.aicc - AICc(rain.week, k=2) #worse than null 

rain.month <- occu( ~ rain.month ~ 1, occu.odhe, starts = c(2, 0, -3))    
null.aicc - AICc(rain.month, k=2) #better than null 
confint(rain.month, level=0.85, type="det") # 85% CI does not overlap zero

# Water source and long-term precipitation interaction

rain.month.water.dense <- occu( ~ water_tank_density_1_9_km * rain.month 
                                ~ 1, occu.odhe)       
null.aicc - AICc(rain.month.water.dense, k=2) #worse than null 

# rain month proceeds

##### Temperature Hypothesis Group ####

temp.max <- occu( ~ max.temp ~ 1, occu.odhe, starts = c(-1, -1, 0))  
null.aicc - AICc(temp.max, k=2) #worse than null 

# Temperature interactions 

temp.water.dense <- occu( ~ max.temp * water_tank_density_1_9_km ~ 1, occu.odhe)   
null.aicc - AICc(temp.water.dense , k=2)  #worse than null 

temp.canopy <- occu( ~ max.temp * canopy_cov ~ 1, occu.odhe,
                       starts = c(-1, -1, 0, -1, -1))       
null.aicc - AICc(temp.canopy, k=2) #better than null 
confint(temp.canopy, level=0.85, type="det") # 85% CI overlaps zero

temp.humid <- occu( ~ max.temp * humid ~ 1, occu.odhe)      
null.aicc - AICc(temp.humid, k=2)  #worse than null

# none proceed

#### Wind Hypothesis Group ####

wind <- occu( ~ wind ~ 1, occu.odhe, starts = c(0, -1, 0))       
null.aicc - AICc(wind, k=2) #worse than null 

# none proceed

##### Predator Activity Hypothesis Group ####

coy.bob.hours <- occu( ~ coy.bob.active ~ 1, occu.odhe, starts = c(-1, 0, 0))
null.aicc - AICc(coy.bob.hours, k=2) #worse than null 

coy.bob.total <- occu( ~ coy.bob.count ~ 1, occu.odhe, starts = c(2, 0, -3)) 
null.aicc - AICc(coy.bob.total,k=2) #worse than null

coy.total <- occu( ~ coyote.count ~ 1, occu.odhe, starts = c(-1, -1, 0))    
null.aicc - AICc(coy.total,k=2) #worse than null 

coy.hours <- occu( ~ coyote.active ~ 1, occu.odhe, starts = c(-1, 0, 0))    
null.aicc - AICc(coy.hours,k=2) #worse than null 

bob.hours <- occu( ~ bobcat.active ~ 1, occu.odhe, starts = c(-1, -1, 0))
null.aicc - AICc(bob.hours,k=2) #worse than null 

bob.total <- occu( ~ bobcat.count ~ 1, occu.odhe, starts = c(-1, 0, 0))  
null.aicc - AICc(bob.total,k=2) #worse than null 

# none proceed

##### Livestock Activity Hypothesis Group####

cow.hours <- occu( ~ cow.active ~ 1, occu.odhe, starts = c(-1, -1, 0))       
null.aicc - AICc(cow.hours, k=2) #better than null 
confint(cow.hours, level=0.85, type="det") #85% CI does not overlap zero

cow.total <- occu( ~ cow.count ~ 1, occu.odhe, starts = c(2, 0, -3))       
null.aicc - AICc(cow.total,k=2) #better than null 
confint(cow.total, level=0.85, type="det") #85% CI does not overlap zero

stock.total <- occu( ~ livestock.count ~ 1, occu.odhe, starts = c(2, 0, -3))
null.aicc - AICc(stock.total,k=2) #better than null 
confint(stock.total, level=0.85, type="det") #85% CI does not overlap zero

stock.hours <- occu( ~ livestock.active ~ 1, occu.odhe, starts = c(2, 0, -3))     
null.aicc - AICc(stock.hours,k=2) #better than null 
confint(stock.hours, level=0.85, type="det") #85% CI does not overlap zero

# Livestock interactions

coy.cow <- occu( ~ coyote.count * livestock.count ~ 1, occu.odhe)     
null.aicc - AICc(coy.total,k=2) #worse than null

cow.water.dense <- occu( ~ livestock.count * water_tank_density_1_9_km  
                         ~ 1, occu.odhe)
null.aicc - AICc(cow.water.dense , k=2) #better than null 
confint(cow.water.dense, level=0.85, type="det") #85% CI overlaps zero

# livestock total proceeds 

##### Human Activity Hypothesis Group#####

vehicle.total <- occu( ~ vehicle.count ~ 1, occu.odhe, starts = c(2, 0, -3))
null.aicc - AICc(vehicle.total,k=2) #worse than null

vehicle.hours <- occu( ~ vehicle.active ~ 1, occu.odhe, starts = c(-1, 0, 0))     
null.aicc - AICc(vehicle.hours,k=2) #worse than null

people.hours <- occu( ~ people.active ~ 1, occu.odhe, starts = c(-1, 0, 0))
null.aicc - AICc(people.hours,k=2) #better than null 
confint(people.hours, level=0.85, type="det")  #85% CI does not overlap zero

human.hours <- occu( ~ human.active ~ 1, occu.odhe, starts = c(-1, -1, 0)) 
null.aicc - AICc(human.hours,k=2) #worse than null

# People hours proceeds 

##### Biotic Community Type Hypothesis Group ####

ndvi <- occu( ~ NDVI_1_9km ~ 1, occu.odhe, starts = c(2, 0, -3))   
null.aicc - AICc(ndvi,k=2) #better than null 
confint(ndvi, level=0.85, type="det") #85% CI does not overlap zero

wood <- occu( ~ woodland_percent_1_9km ~ 1, occu.odhe)          
null.aicc - AICc(wood, k=2) #better than null 
confint(wood, level=0.85, type="det") #85% CI does not overlap zero

bio.com <- occu( ~ as.factor(biotic_com_2) ~ 1, occu.odhe)          
null.aicc - AICc(bio.com, k=2) #better than null 
confint(bio.com, level=0.85, type="det") #85% CI does not overlap zero

tree <- occu( ~ tree_density_5_70 ~ 1, occu.odhe)
null.aicc - AICc(tree)  #better than null
confint(tree, level=0.85, type="det") #85% CI does not overlap zero

# NDVI proceeds

##### Vegetation Concealment Cover Hypothesis Group ####

obs.distance <- occu( ~ obs_dist ~ 1, occu.odhe)      
null.aicc - AICc(obs.distance, k=2) #better than null
confint(obs.distance, level=0.85, type="det") # 85% CI does not overlap zero

veg.cover.cam <- occu( ~ veg_cover_cam ~ 1, occu.odhe) 
null.aicc - AICc(veg.cover.cam, k=2) #better than null
confint(veg.cover.cam, level=0.85, type="det") # 85% CI does not overlap zero

detect.angle <- occu( ~ detection_angle ~ 1, occu.odhe, starts = c(-1, -1, 0))    
null.aicc - AICc(detect.angle, k=2) #worse than null

max.trig.dist <- occu( ~ max_trig_dist ~ 1, occu.odhe, starts = c(-1, 0, 0))   
null.aicc - AICc(max.trig.dist, k=2)#better than null
confint(max.trig.dist, level=0.85, type="det") # 85% CI does not overlap zero

# max trigger distance proceeds

#### Sampling Period Hypothesis Group ####

### Sampling Period Length ###
samp.period <- occu( ~ sampling_period_length 
                     ~ 1, occu.odhe, starts = c(2, 0, -3))          
null.aicc - AICc(samp.period,k=2)  #worse than null

### Camera Moved ###
cam.moved <- occu( ~ as.factor(cam_moved) ~ 1, occu.odhe, starts = c(2, 0, -3))  
null.aicc - AICc(cam.moved, k=2) #worse than null

# none proceed

#################################################
# DETECTION MODEL SELECTION #
#################################################

# Variables that proceed (5)
  # NDVI
  # max trigger dist
  # livestock count
  # people active
  # rain month

# Check correlations between detection variables 

# site-level variables 
site.covs.cor <- site.covs %>% 
  select(max_trig_dist, NDVI_1_9km)

cor_site_covs <- cor(site.covs.cor, method='spearman')  
# none correlated above |0.7|

#  observation-level variables

obs_cor <-cor(data.frame(
  livestock = as.vector(obsCovs$livestock.count),
  human = as.vector(obsCovs$people.active),
  rain = as.vector(obsCovs$rain.month)
))
# none correlated above |0.7|

# Run this code to get correlation matrix in Excel 
# write.xlsx(cor_obs, file="Correlations mule deer.xlsx") 

# Null model
mod_null <- occu( ~ 1 ~ 1, occu.odhe)

# Single variable models 
mod1 <- occu( ~ max_trig_dist ~ 1, occu.odhe)
mod2 <- occu( ~ NDVI_1_9km ~ 1, occu.odhe)
mod3 <- occu( ~ people.active ~ 1, occu.odhe)
mod4 <- occu( ~ livestock.count ~ 1, occu.odhe)
mod5 <- occu( ~ rain.month ~ 1, occu.odhe)

# 2 variable models 
mod6 <- occu( ~ max_trig_dist + NDVI_1_9km 
              ~ 1, occu.odhe, starts = c(2, 5, 0, -3))
mod7 <- occu( ~ max_trig_dist + people.active ~ 1, occu.odhe)
mod8 <- occu( ~ max_trig_dist + livestock.count ~ 1, occu.odhe)
mod9 <- occu( ~ max_trig_dist + rain.month ~ 1, occu.odhe)
mod10 <- occu( ~ NDVI_1_9km + people.active ~ 1, occu.odhe)
mod11 <- occu( ~ NDVI_1_9km + livestock.count ~ 1, occu.odhe)
mod12 <- occu( ~ NDVI_1_9km + rain.month ~ 1, occu.odhe)
mod13 <- occu( ~ people.active + livestock.count ~ 1, occu.odhe)
mod14 <- occu( ~ people.active + rain.month ~ 1, occu.odhe)
mod15 <- occu( ~ livestock.count + rain.month ~ 1, occu.odhe)

# Model selection
top_mods <- model.sel(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod10, 
                      mod11, mod12, mod13, mod14, mod15, mod_null)

# Candidate detection models are found in Table S2.5.

# Run this code to see candidate detection models 
 #write.xlsx(top_mods, file="Mule Deer Base Models.xlsx", 
         #  sheetName="Detection Models", append=T)  

# Top model diagnostics

summary(mod6)
null.aicc - AICc(mod6)

# Calculate the 85% confidence intervals for variables 
confint(mod6, type = "det", level = 0.85)    
#CIs for detection variables - none overlap zero

confint(mod6, type = "state", level = 0.85)  
#CIs for occupancy intercept - overlaps zero

# Calculate multicollinearity

unmarked::vif(mod6, type = "det") # no colliniarity - all below 2

#################################################
# CREATE THE HABITAT SELECTION (OCCUPANCY) MODEL#
#################################################

## Mule deer habitat selection hypotheses are listed in Table S2.2.

## Table S2.4 lists which occupancy hypothesis groups proceeded (i.e. performed 
  # better than the null occupancy model ψ(.), did not have uninformative 
  # parameters, and supported our original hypothesis).

## Candidate habitat selection (occupancy) models are found in Table S2.6.
  # Models were limited to 6 parameters and were derived from the best-supported
  # detection model (mod6 above).

occu.null.aicc <- AICc(mod6)

#### Water Hypothesis Group ####

water.dist <- occu( ~ max_trig_dist + NDVI_1_9km 
                         ~ dist_water_tank, occu.odhe)
occu.null.aicc - AICc(water.dist) #worse than null 

water.tank.dense <- occu( ~ max_trig_dist + NDVI_1_9km 
                          ~ water_tank_density_1_9_km, occu.odhe)       
occu.null.aicc - AICc(water.tank) #worse than null 

#none proceed

#### Predator Activity Hypothesis Group####

bobcat.coyote <- occu( ~ max_trig_dist + NDVI_1_9km
                       ~ pred_count_avg, occu.odhe)
occu.null.aicc - AICc(bobcat.coyote, k=2) #worse than null 

coyote <- occu( ~ max_trig_dist + NDVI_1_9km 
                ~ coy_count_avg, occu.odhe)
occu.null.aicc - AICc(coyote, k=2) #worse than null 

# none proceed 

#### Livestock Activity Hypothesis Group####

stock <- occu( ~ max_trig_dist + NDVI_1_9km 
               ~ livestock_count_avg, occu.odhe)
occu.null.aicc - AICc(stock, k=2) #better than null
confint(stock, type = "state", level = 0.85) #85% CI does not overlap zero

cow <- occu( ~ max_trig_dist + NDVI_1_9km 
             ~ cow_count_avg, occu.odhe)
occu.null.aicc - AICc(cow, k=2) #better than null
confint(cow, type = "state", level = 0.85) #85% CI does not overlap zero

# livestock count proceeds

#### Biotic Community Type Hypothesis Group ####

trees <- occu( ~  max_trig_dist + NDVI_1_9km 
               ~  tree_density_5_70, occu.odhe)
occu.null.aicc - AICc(trees) #better than null
confint(trees, level=0.85, type="state") #85% CI does not overlap zero

ndvi <- occu( ~ max_trig_dist + NDVI_1_9km 
              ~ NDVI_1_9km, occu.odhe)
occu.null.aicc - AICc(ndvi) #better than null
confint(ndvi, level=0.85, type="state") #85% CI does not overlap zero

shrubland <- occu( ~ max_trig_dist + NDVI_1_9km 
                   ~ shrub_yucca_density, occu.odhe)
occu.null.aicc - AICc(shrubland) #better than null
confint(shrubland, level=0.85, type="state") #85% CI does not overlap zero

grassland <- occu( ~ max_trig_dist + NDVI_1_9km 
                   ~ herbaceous_cov, occu.odhe)
occu.null.aicc - AICc(grassland) #better than null
confint(grassland, level=0.85, type="state") #85% CI does not overlap zero

wood <- occu( ~ max_trig_dist + NDVI_1_9km 
              ~ woodland_percent_1_9km, occu.odhe)          
occu.null.aicc - AICc(wood, k=2) #better than null
confint(wood, level=0.85, type="state") #85% CI does not overlap zero

bio.com <- occu( ~ max_trig_dist + NDVI_1_9km
                 ~ as.factor(biotic_com_2), occu.odhe)          
occu.null.aicc - AICc(bio.com, k=2) #better than null
confint(bio.com, level=0.85, type="state") #85% CI does not overlap zero

# Shrubland proceeds

#### Vegetation Concealment Cover Hypothesis Group  ####

vert.cover <- occu( ~ max_trig_dist + NDVI_1_9km 
                    ~ vertical_cover, occu.odhe)
occu.null.aicc - AICc(vert.cover)  #worse than null 

veg.cover <- occu( ~ max_trig_dist + NDVI_1_9km 
                   ~ veg_cover, occu.odhe)
occu.null.aicc - AICc(veg.cover) #better than null
confint(veg.cover, type = "state", level = 0.85) #85% CI does not overlap zero

# veg cov proceeds 

##### Topography Hypothesis Group ####

topo.pos <- occu( ~ max_trig_dist + NDVI_1_9km 
                  ~ topo_pos, occu.odhe)
occu.null.aicc - AICc(topo.pos) #worse than null 

slope <- occu( ~ max_trig_dist + NDVI_1_9km 
                  ~  slope, occu.odhe)
occu.null.aicc - AICc(slope) #better than null
confint(slope, level=0.85, type="state") #85% CI does not overlap zero

elev <- occu( ~ max_trig_dist + NDVI_1_9km 
              ~ elevation, occu.odhe)
occu.null.aicc - AICc(elev)  #worse than null 

# Slope proceeds

##### Coolness of Site Hypothesis Group ####

aspect <- occu( ~ max_trig_dist + NDVI_1_9km 
                ~ aspect, occu.odhe)
occu.null.aicc - AICc(aspect) #worse than null 

heat.load <- occu( ~ max_trig_dist + NDVI_1_9km 
                   ~  heat_load, occu.odhe)
occu.null.aicc - AICc(heat.load) #worse than null 

#none proceed

#################################################
# OCCUPANCY MODEL SELECTION #
#################################################

# Variables that proceed (4)
  # livestock_count_avg
  # shrub_yucca_density
  # veg_cover
  # slope

# Check correlations between occupancy variables 
site.covs.cor <- site.covs %>% 
  select(livestock_count_avg, shrub_yucca_density, 
         veg_cover, slope)

cor_site_covs <- cor(site.covs.cor, method='spearman') 
# none correlated above |0.7|

# Detection only model 
det.mod <- occu( ~ max_trig_dist + NDVI_1_9km  ~ 1, occu.odhe)

# 1-variable occupancy models
mod1 <- occu( ~ max_trig_dist + NDVI_1_9km ~ livestock_count_avg, occu.odhe)
mod2 <- occu( ~ max_trig_dist + NDVI_1_9km ~ slope, occu.odhe)
mod3 <- occu( ~ max_trig_dist + NDVI_1_9km ~ shrub_yucca_density, occu.odhe)
mod4 <- occu( ~ max_trig_dist + NDVI_1_9km ~ veg_cover, occu.odhe)

# 2-variable occupancy models
mod5 <- occu( ~ max_trig_dist + NDVI_1_9km 
              ~ livestock_count_avg + slope, occu.odhe)
mod6 <- occu( ~ max_trig_dist + NDVI_1_9km 
              ~ livestock_count_avg + shrub_yucca_density, occu.odhe)
mod7 <- occu( ~ max_trig_dist + NDVI_1_9km 
              ~ livestock_count_avg + veg_cover, occu.odhe)
mod8 <- occu( ~ max_trig_dist + NDVI_1_9km 
              ~ slope + shrub_yucca_density, occu.odhe)
mod9 <- occu( ~ max_trig_dist + NDVI_1_9km 
              ~ slope + veg_cover, occu.odhe)
mod10 <- occu( ~ max_trig_dist + NDVI_1_9km 
               ~ shrub_yucca_density + veg_cover, occu.odhe)

# Combine models in model selection table
top_mods <- model.sel(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9,
                      mod10, det.mod)

# Candidate occupancy models are found in Table S2.6.

# Run this code to see candidate occupancy models 
#write.xlsx(top_mods, file="Mule Deer Base Models.xlsx", 
        #   sheetName="Occupancy Models", append=T)  

# Top Model Diagnostics 

AICc(det.mod) - AICc(mod10)
summary(mod10)

# The parameter estimates, standard errors, and 85% confidence intervals for
  #the best-supported model describing the probability of detection (p) and 
  #habitat selection (ψ) of mule deer are listed in Table S2.7.

# Calculate the 85% confidence intervals for variables 
confint(mod10, type = "det", level = 0.85)   
#CIs for detection variables do not overlap zero

confint(mod10, type = "state", level = 0.85)
#CIs for occupancy variables do not overlap zero

# Calculate multicollinearity

unmarked::vif(mod10, type = "det") # no colliniarity - all below 2
unmarked::vif(mod10, type = "state") # no colliniarity - all below 2

# Calculate mean daily detection probability with the top detection model

backTransform(linearComb(mod10,                           
                         coefficients= c(1,0,0),  
                         type = 'det'))     
#0.0642 or 6.4%

# Calculate the overall detection probability if K = 29
1-(1-0.0642)^29

# 85% probability of detecting a mule deer at least once given the minimum 
  #sampling period length at a site

# Calculate the overall detection probability if K = 38 
1-(1-0.0642)^38

# 92% probability of detecting a mule deer at least once given the average
  #sampling period length at a site

# Calculate the occupancy probability

backTransform(linearComb(mod10,                
                         coefficients=c(1,0,0),   
                         type = 'state')) 
#0.748 OR 74.8%


# Predict daily detection probability for survey sites 

pred_det <- predict(mod10,          
                      type = "det",                 
                      newdata = occu.odhe@siteCovs)[c("Predicted",
                                                      "SE",
                                                      "lower",
                                                      "upper")]

pred_det_df <- data.frame(Predicted = pred_det$Predicted,
                            StandardError = pred_det$SE,
                            lower = pred_det$lower,
                            upper = pred_det$upper,
                            site.covs)

# Predict occupancy probability for survey sites 

pred_occu <- predict(mod10,          
                     type = "state",                 
                     newdata = occu.odhe@siteCovs)[c("Predicted",
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