###########################################################
######## Code for: base models for American badgers #######
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
# SETUP CODE FOR AMERICAN BADGER OCCUPANCY MODELS #
###########################################################

# Read in edited .csv
detHist <- read.csv(file = "badger detection hist.csv", 
                    row.names = 1)

# Change from integer to numeric

detHist <- detHist %>%
  mutate(across(1:51, as.numeric))

# Create detection history 
occu.tata <- unmarkedFrameOccu(y=detHist, 
                               siteCovs = site.covs.scaled, 
                               obsCovs = obsCovs.scaled)

# Null model 
tata.null <- occu( ~ 1 ~ 1, occu.tata, linkPsi="logit", se=TRUE,
                  control = list(maxit = 10000))
null.aicc <- AICc(tata.null)

# Null detection probability 
backTransform(tata.null['det'])

# Null occupancy probability 
backTransform(tata.null['state'])

#################################################
# CREATE THE DETECTION MODEL #
#################################################

## American badger detection hypotheses are listed in Table S2.1.

## Table S2.3 lists which detection hypothesis groups proceeded (i.e. performed 
  # better than the null detection model p(.), did not have uninformative 
  # parameters, and supported our original hypothesis).

## Candidate detection models are found in Table S2.5.
  # Models were limited to 4 parameters total and were tested with no variables 
  # on the probability of habitat selection (ψ(.)).

##### Water Hypothesis Group ####

water.tank <- occu( ~ water_tank_density_1_6_km 
                    ~ 1, occu.tata, starts = c(-1, 0, 0))       
null.aicc - AICc(water.tank) #worse than null

water.dist.tank <- occu( ~ dist_water_tank 
                         ~ 1, occu.tata, starts = c(-1, 0, 0))      
null.aicc - AICc(water.dist.tank, k=2)  #worse than null

#none proceed 

##### Precipitation Hypothesis Group####

# Daily precipitation 

precip.cat <- occu( ~ as.factor(precip.cat) ~ 1, occu.tata, starts = c(-1, 0, 0))
null.aicc - AICc(precip.cat) #better than null
confint(precip.cat, level=0.85, type="det") # 85% CI does not overlap zero 

precip.cm <- occu( ~ precip.cm ~ 1, occu.tata, starts = c(-1, 0, 0))
null.aicc - AICc(precip.cm) #better than null
confint(precip.cm, level=0.85, type="det") # 85% CI does not overlap zero 

# Long-term precipitation

days.since.rain <- occu( ~ days.since.rain ~ 1 , occu.tata)    
null.aicc - AICc(days.since.rain, k=2) #better than null
confint(days.since.rain, level=0.85, type="det") # 85% CI does not overlap zero 

rain.month<- occu( ~ rain.month ~ 1, occu.tata)    
null.aicc - AICc(rain.month, k=2) #better than null
confint(rain.month, level=0.85, type="det") # 85% CI does not overlap zero

rain.week <- occu( ~ rain.week ~ 1, occu.tata)    
null.aicc - AICc(rain.week, k=2) #worse than null 

# Precipitation interaction
precip.water.dense <- occu( ~ water_tank_density_1_6_km * 
                              days.since.rain ~ 1, occu.tata)       
null.aicc - AICc(precip.water.dense , k=2)  #better than null
confint(precip.water.dense, level=0.85, type="det") # 85% CI overlaps zero

# precip cat proceeds

#### Humidity Hypothesis Group ####

humidity <- occu( ~ humid ~ 1, occu.tata, starts = c(-1, 0, 0))       
null.aicc - AICc(humidity, k=2) #better than null
confint(humidity, level=0.85, type="det") #85 CI does not overlap zero

# humidity proceeds

##### Temperature Hypothesis Group ####

temp.max <- occu( ~ max.temp ~ 1, occu.tata, starts = c(-1, -1, 0))    
null.aicc - AICc(temp.max, k=2)  #worse than null 

# Temperature interactions

temp.water.dense <- occu( ~ max.temp * water_tank_density_1_6_km ~ 1, occu.tata)
null.aicc - AICc(temp.water.dense, k=2) #worse than null 

temp.canopy <- occu( ~ max.temp * canopy_cov ~ 1, occu.tata)      
null.aicc - AICc(temp.canopy, k=2) #worse than null 

temp.humid <- occu( ~ max.temp * humid ~ 1, occu.tata)      
null.aicc - AICc(temp.humid, k=2) #better than null
confint(temp.humid, level=0.85, type="det") 
# 85% CI does not overlap zero but is positive (no support for hypothesis)

# none proceed

#### Wind Hypothesis Group ####

wind <- occu( ~ wind ~ 1, occu.tata)       
null.aicc - AICc(wind, k=2) #worse than null 

# none proceed

##### Predator Activity Hypothesis Group ####

coy.hours <- occu( ~ coyote.active ~ 1, occu.tata, starts = c(-1, -1, 0)) 
null.aicc - AICc(coy.hours,k=2) #worse than null

coy.total <- occu( ~ coyote.count ~ 1, occu.tata, starts = c(-1, -1, 0)) 
null.aicc - AICc(coy.total,k=2) #worse than null

#none proceed

#### Prey Activity Hypothesis Group ####

lagomorph.total <- occu( ~ lagomorph.count 
                         ~ 1, occu.tata, starts = c(-1, -1, 0))
null.aicc - AICc(lagomorph.total,k=2) #worse than null

lagomorph.hours <- occu( ~ lagomorph.active 
                         ~ 1, occu.tata, starts = c(-1, -1, -1))     
null.aicc - AICc(lagomorph.hours,k=2) #worse than null

rodent.total <- occu( ~ rodent.count ~ 1, occu.tata, starts = c(-1, 0, 0))
null.aicc - AICc(rodent.total,k=2) #worse than null

rodent.hours <- occu( ~ rodent.active ~ 1, occu.tata, starts = c(-1, -1, -1))     
null.aicc - AICc(rodent.hours,k=2) #better than null
confint(rodent.hours, level=0.85, type="det") #85 CI does not overlap zero

lago.rodent.total <- occu( ~ lago.rodent.count 
                           ~ 1, occu.tata, starts = c(-1, -1, 0))
null.aicc - AICc(lago.rodent.total,k=2) #worse than null

lago.rodent.hours <- occu( ~ lago.rodent.active 
                           ~ 1, occu.tata, starts = c(-1, -1, -1))     
null.aicc - AICc(lago.rodent.hours,k=2) #worse than null

cottontail.total <- occu( ~ cottontail.count 
                          ~ 1, occu.tata, starts = c(-1, -1, 0))
null.aicc - AICc(cottontail.total,k=2) #worse than null

cottontail.hours <- occu( ~ cottontail.active 
                          ~ 1, occu.tata, starts = c(-1, -1, -1))     
null.aicc - AICc(cottontail.hours,k=2) #worse than null

jackrabbit.total <- occu( ~ jackrabbit.count 
                          ~ 1, occu.tata, starts = c(-1, -1, 0))
null.aicc - AICc(jackrabbit.total,k=2)#worse than null

jackrabbit.hours <- occu( ~ jackrabbit.active 
                          ~ 1, occu.tata, starts = c(-1, -1, -1))     
null.aicc - AICc(jackrabbit.hours,k=2) #worse than null

# rodent hours proceeds

##### Human Activity Hypothesis Group#####

vehicle.total <- occu( ~ vehicle.count ~ 1, occu.tata, starts = c(-1, -1, 0))
null.aicc - AICc(vehicle.total,k=2) #worse than null

vehicle.hours <- occu( ~ vehicle.active ~ 1, occu.tata, starts = c(-1, -1, -1))     
null.aicc - AICc(vehicle.hours,k=2) #worse than null

people.hours <- occu( ~ people.active ~ 1, occu.tata, starts = c(-1, 0, 0))
null.aicc - AICc(people.hours,k=2) #worse than null

human.hours <- occu( ~ human.active ~ 1, occu.tata, starts = c(-1, 0, 0)) 
null.aicc - AICc(human.hours,k=2) #worse than null

# none proceed

##### Biotic Community Type Hypothesis Group ####

ndvi.site <- occu( ~ NDVI_1_6km ~ 1, occu.tata)          
null.aicc - AICc(ndvi.site,k=2) #better than null
confint(ndvi.site, level=0.85, type="det") #85 CI does not overlap zero

wood <- occu( ~ woodland_percent_1_6km ~ 1, occu.tata)          
null.aicc - AICc(wood, k=2) #better than null
confint(wood, level=0.85, type="det") #85 CI does not overlap zero

bio.com <- occu( ~ as.factor(biotic_com_2) ~ 1, occu.tata, 
                   starts = c(0, -4, 0))        
null.aicc - AICc(bio.com)  #worse than null

tree <- occu( ~ tree_density_5_70 ~ 1, occu.tata)
null.aicc - AICc(tree)  #worse than null

#NDVI proceeds

##### Vegetation Concealment Cover Hypothesis Group ####

obs.distance <- occu( ~ obs_dist ~ 1, occu.tata)       
null.aicc - AICc(obs.distance, k=2) #better than null
confint(obs.distance, level=0.85, type="det")  #85 CI does not overlap zero

veg.cov.cam <- occu( ~ veg_cover_cam ~ 1, occu.tata)       
null.aicc - AICc(veg.cov.cam, k=2) #better than null
confint(veg.cov.cam, level=0.85, type="det") #85 CI does not overlap zero

detect.angle <- occu( ~ detection_angle ~ 1, occu.tata, starts = c(-1, 0, 0))
null.aicc - AICc(detect.angle, k=2) #better than null
confint(detect.angle, level=0.85, type="det") #85 CI does not overlap zero

max.trig.dist <- occu( ~ max_trig_dist ~ 1, occu.tata, starts = c(-1, 0, 0)) 
null.aicc - AICc(max.trig.dist, k=2) #worse than null

# detection angle proceeds

#### Sampling Period Hypothesis Group ####

### Sampling Period Length ###
samp.period <- occu( ~ sampling_period_length ~ 1, occu.tata)          
null.aicc - AICc(samp.period,k=2) #worse than null

### Camera Moved ###
cam.moved <- occu( ~ as.factor(cam_moved) ~ 1, occu.tata)   
null.aicc - AICc(cam.moved, k=2)  #worse than null

## Camera Type ##
cam.type <- occu( ~ as.factor(badger.cam) ~ 1, occu.tata) 
null.aicc - AICc(cam.type, k=2) #worse than null

# none proceed

#################################################
# DETECTION MODEL SELECTION #
#################################################

# Variables that proceed (5)
  #rodent activity
  #detection_angle
  #NDVI_1_6km
  #humidity 
  #precip cat

# Null model
tata.null <- occu( ~ 1 ~ 1, occu.tata)

#single variable models 
mod1 <- occu( ~ rodent.active ~ 1, occu.tata)
mod2 <- occu( ~ NDVI_1_6km ~ 1, occu.tata)
mod3 <- occu( ~ detection_angle ~ 1, occu.tata)
mod4 <- occu( ~ humid ~ 1, occu.tata)
mod5 <- occu( ~ as.factor(precip.cat) ~ 1, occu.tata)

#2- variable models 
mod6 <- occu( ~ rodent.active + NDVI_1_6km ~ 1, occu.tata)
mod7 <- occu( ~ rodent.active + detection_angle ~ 1, occu.tata)
mod8 <- occu( ~ rodent.active + humid ~ 1, occu.tata)
mod9 <- occu( ~ rodent.active + as.factor(precip.cat) ~ 1, occu.tata)
mod10 <- occu( ~ NDVI_1_6km + detection_angle ~ 1, occu.tata)
mod11 <- occu( ~ NDVI_1_6km + humid ~ 1, occu.tata)
mod12 <- occu( ~ NDVI_1_6km + as.factor(precip.cat) ~ 1, occu.tata)
mod13 <- occu( ~ detection_angle + humid ~ 1, occu.tata)
mod14 <- occu( ~ detection_angle + as.factor(precip.cat) ~ 1, occu.tata)
mod15 <- occu( ~ humid + as.factor(precip.cat) ~ 1, occu.tata)

# Model selection
top_mods <- model.sel(
  mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod12, 
  mod13, mod14, mod15, tata.null
)

# Candidate detection models are found in Table S2.5.

# Run this code to see candidate detection models 
# write.xlsx(top_mods, file="Badger Base Models.xlsx", 
           #sheetName="Detection Models", append=T)

# Top model diagnostics

summary(mod14)
AICc(tata.null) - AICc(mod14) 

# Calculate multicollinearity
unmarked::vif(mod14, type = "det") # all below 2

### Calculating 85% confidence intervals for variables 
confint(mod14, type = "det", level = 0.85)    
#CIs for detection variables - none overlap zero

confint(mod14, type = "state", level = 0.85)  
#CIs for occupancy intercept - does not overlap zero

#################################################
# CREATE THE HABITAT SELECTION (OCCUPANCY) MODEL#
#################################################

## American badger habitat selection hypotheses are listed in Table S2.2.

## Table S2.4 lists which occupancy hypothesis groups proceeded (i.e. performed 
  # better than the null occupancy model ψ(.), did not have uninformative 
  # parameters, and supported our original hypothesis).

## Candidate habitat selection (occupancy) models are found in Table S2.6.
  # Models were limited to 6 parameters and were derived from the best-supported
  # detection model (mod14 above).

occu.null.aicc <- AICc(mod14)

####Water Hypothesis Group#####

water.tank <- occu( ~ detection_angle + as.factor(precip.cat)
                    ~ water_tank_density_1_6_km, occu.tata)       
occu.null.aicc - AICc(water.tank)  #worse than null

water.dist.tank <- occu(~ detection_angle + as.factor(precip.cat)
                        ~ dist_water_tank, occu.tata)      
occu.null.aicc - AICc(water.dist.tank, k=2)  #worse than null

#none proceed

#### Prey Activity Hypothesis Group####

rodent <- occu( ~ detection_angle + as.factor(precip.cat)
                ~ rodent_count_avg, occu.tata)
occu.null.aicc - AICc(rodent, k=2) #worse than null

#none proceed

#### Vegetation Concealment Cover Hypothesis Group  ####

vert.cover <- occu( ~ detection_angle + as.factor(precip.cat)
                    ~ vertical_cover, occu.tata)
occu.null.aicc - AICc(vert.cover) #better than null
confint(vert.cover, level=0.85, type="state") #85% CI does not overlap zero

veg.cov <- occu( ~ detection_angle + as.factor(precip.cat)
                 ~ veg_cover, occu.tata)
occu.null.aicc - AICc(veg.cov) #worse than null 

# vertical cover proceeds

#### Biotic Community Type Hypothesis Group ####

woodland <- occu( ~ detection_angle + as.factor(precip.cat)
                  ~ tree_density_5_70, occu.tata)
occu.null.aicc - AICc(woodland) #worse than null 

ndvi.site <- occu( ~ detection_angle + as.factor(precip.cat)  
                   ~ NDVI_1_6km, occu.tata)          
occu.null.aicc - AICc(ndvi.site,k=2) #worse than null 

shrubland <- occu( ~ detection_angle + as.factor(precip.cat)
                   ~ shrub_yucca_density, occu.tata)
occu.null.aicc - AICc(shrubland) #worse than null 

grassland <- occu( ~ detection_angle + as.factor(precip.cat)
                   ~ herbaceous_cov, occu.tata)
occu.null.aicc - AICc(grassland) #worse than null 

wood <- occu( ~ detection_angle + as.factor(precip.cat)
              ~ woodland_percent_1_6km, occu.tata)          
occu.null.aicc - AICc(wood, k=2) #worse than null

bio.com <- occu( ~ detection_angle + as.factor(precip.cat)
                 ~ as.factor(biotic_com_2), occu.tata, 
                   starts = c(0, 1, -4, -0.25, 0.5))        
occu.null.aicc - AICc(bio.com, k=2)  #worse than null

#none proceed 

##### Topography Hypothesis Group ####

topo.pos <- occu( ~ detection_angle + as.factor(precip.cat)
                  ~ topo_pos, occu.tata)
occu.null.aicc - AICc(topo.pos) #worse than null 

slope <- occu( ~ detection_angle + as.factor(precip.cat)
               ~ slope, occu.tata)
occu.null.aicc - AICc(slope) #better than null
confint(slope, level=0.85, type="state") #85% CI does not overlap zero

elev <- occu( ~ detection_angle + as.factor(precip.cat)
              ~ elevation, occu.tata)
occu.null.aicc - AICc(elev) #worse than null

#slope proceeds 

##### Coolness of Site Hypothesis Group ####

aspect <- occu( ~ detection_angle + as.factor(precip.cat)
                ~ aspect, occu.tata)
occu.null.aicc - AICc(aspect)   #worse than null 

heat.load <- occu( ~ detection_angle + as.factor(precip.cat)
                   ~ heat_load, occu.tata)
occu.null.aicc - AICc(heat.load) #worse than null

#none proceeds

#################################################
# OCCUPANCY MODEL SELECTION #
#################################################

# Variables that proceed (2)
  # vertical_cover
  # slope

# Check correlations between occupancy variables 
site.covs.cor <- site.covs %>% 
  select(slope,  vertical_cover)

cor_site_covs <- cor(site.covs.cor, method='spearman')
# none correlated above |0.7|

# Detection only model 
det.mod <- occu( ~ detection_angle + as.factor(precip.cat) ~ 1, occu.tata)

# single variable mods
mod1 <- occu( ~ detection_angle + as.factor(precip.cat) ~  slope, occu.tata)
mod2 <- occu( ~ detection_angle + as.factor(precip.cat)
              ~  vertical_cover, occu.tata)

#2-variable mods
mod3 <- occu( ~ detection_angle + as.factor(precip.cat)
              ~  vertical_cover + slope, occu.tata)

# Combine models in model selection table
top_mods <- model.sel(det.mod, mod1, mod2, mod3)

# Candidate occupancy models are found in Table S2.6.

# Run this code to see candidate occupancy models 
#write.xlsx(top_mods, file="Badger Base Models.xlsx", 
         # sheetName="Occupancy Models", append=T)

# Top model diagnostics

summary(mod3)
AICc(det.mod) - AICc(mod3)

# The parameter estimates, standard errors, and 85% confidence intervals for
  #the best-supported model describing the probability of detection (p) and 
  #habitat selection (ψ) of American badgers are listed in Table S2.12.

### Calculating 85% confidence intervals for variables 
confint(mod3, type = "det", level = 0.85)    
#CIs for detection variables do not overlap zero

confint(mod3, type = "state", level = 0.85)  
#CIs for occupancy intercept and vertical cover overlap zero 

# Calculate multicollinearity
unmarked::vif(mod3, type = "det") # all below 2
unmarked::vif(mod3, type = "state") # all below 2

# Calculate mean daily detection probability with the top detection model

backTransform(linearComb(mod3,                           
                         coefficients= c(1,0,0),  
                         type = 'det'))     
#0.0264  or 2.6%

### Occupancy probability

backTransform(linearComb(mod3,                
                         coefficients=c(1,0,0),   
                         type = 'state')) 
#0.379 OR 37.9%

### Calculating overall detection probability if K = 29
1-(1-0.0264)^29

#54% probability of detecting an American badger at least once given the minimum
  #sampling period length at a site

### Calculating overall detection probability if K = 38
1-(1-0.0264)^38

# 64% probability of detecting an American badger at least once given the average
  #sampling period length at a site

#Predict detection probability for survey sites 

p_full <- predict(mod3, type = "det")

# Get number of sites and occasions
n_sites <- numSites(occu.tata)
n_occasions <- ncol(getY(occu.tata))
site_index <- rep(1:n_sites, each = n_occasions)

# Add site info and summarize
p_site <- p_full %>%
  mutate(site = site_index) %>%
  group_by(site) %>%
  summarise(
    mean.p = mean(Predicted, na.rm = TRUE),
    mean.SE = sqrt(sum(SE^2, na.rm = TRUE)) / n(),
    lower85 = mean.p - qnorm(0.925) * mean.SE,
    upper85 = mean.p + qnorm(0.925) * mean.SE
  )

#Predict occupancy probability for survey sites 

pred_occu <- predict(mod3,          
                     type = "state",                 
                     newdata = occu.tata@siteCovs)[c("Predicted",
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
