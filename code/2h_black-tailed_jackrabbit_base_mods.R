###########################################################
## Code for: Turbine visibility is a strong predictor of ##
## altered habitat selection by terrestrial mammals at a ##
######## wind energy facility in central New Mexico #######
###########################################################
### This script is for the creation of the black-tailed ###
################ jackrabbit base models ###################
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
# SETUP CODE FOR BLACK-TAILED JACKRABBITS OCCUPANCY MODELS #
###########################################################

# Read in edited .csv
detHist <- read.csv(file = "jackrabbit detection hist.csv", row.names = 1)
# Read in site.covs.scaled
site.covs.scaled <- readRDS("site_covs_scaled.RData")
# Read in obsCovs.scaled
obsCovs.scaled <- readRDS("obsCovs_scaled.RData")

# Change from integer to numeric
detHist %>%
  mutate(across(1:51, as.numeric))

# Create detection history 
occu.leca <- unmarkedFrameOccu(y=detHist, 
                               siteCovs = site.covs.scaled, 
                               obsCovs = obsCovs.scaled)

#### Null model ####
leca.null <- occu( ~ 1 ~ 1, occu.leca , linkPsi="logit", se=TRUE,
                  control = list(maxit = 10000))
null.aicc <- AICc(leca.null)

# Null detection probability (Table 2)
backTransform(leca.null['det'])

### Calculating overall detection probability if K = 38
1-(1-0.196)^38

# 99% probability of detecting a gray fox at least once given the average
#sampling period length at a site and the null detection model (see Table 2)

# Null occupancy probability 
backTransform(leca.null['state'])

#################################################
# CREATE THE DETECTION MODEL #
#################################################

## Black-tailed jackrabbit detection hypotheses are listed in Table S2.1.

## Table S2.3 lists which detection hypothesis groups proceeded (i.e. performed 
  # better than the null detection model p(.), did not have uninformative 
  # parameters, and supported our original hypothesis).

## Candidate detection models are found in Table S2.5.
  # Models were limited to 4 parameters total and were tested with no variables 
  # on the probability of habitat selection (ψ(.)).

##### Water Hypothesis Group####

water.dist.tank <- occu( ~ dist_water_tank ~ 1, occu.leca, starts = c(-1, 0, 0))      
null.aicc - AICc(water.dist.tank, k=2) #better than null
confint(water.dist.tank, level=0.85, type="det") #85% CI does not overlap zero 

water.tank.dense <- occu( ~ water_tank_density_0_9km
                    ~ 1, occu.leca, starts = c(-1, 0, 0))       
null.aicc - AICc(water.tank.dense) #worse than null

# dist water tank proceeds

##### Precipitation Hypothesis Group####

# Daily precipitation 

precip.cat <- occu( ~ as.factor(precip.cat) ~ 1, occu.leca, starts = c(-1, 0, 0))
null.aicc - AICc(precip.cat) #better than null
confint(precip.cat, level=0.85, type="det") #85% CI does not overlap zero 

precip.cm <- occu( ~ precip.cm ~ 1, occu.leca, starts = c(-1, 0, 0))
null.aicc - AICc(precip.cm) #better than null
confint(precip.cm, level=0.85, type="det") #85% CI does not overlap zero 

# Long-term precipitation

days.since.rain <- occu( ~ days.since.rain ~ 1, occu.leca, starts = c(-1, 0, 0))    
null.aicc - AICc(days.since.rain, k=2) #better than null
confint(days.since.rain, level=0.85, type="det") #85% CI does not overlap zero 

rain.month <- occu( ~ rain.month ~ 1, occu.leca, starts = c(-1, 0, 0))    
null.aicc - AICc(rain.month, k=2) #better than null
confint(rain.month, level=0.85, type="det") #85% CI does not overlap zero 

rain.week <- occu( ~ rain.week ~ 1, occu.leca, starts = c(-1, 0, 0))    
null.aicc - AICc(rain.week, k=2) #better than null
confint(rain.week, level=0.85, type="det") #85% CI does not overlap zero 

#  Water source and long-term precipitation interaction

precip.water.dense <- occu( ~ water_tank_density_0_9km * rain.month 
                            ~ 1, occu.leca)       
null.aicc - AICc(precip.water.dense , k=2)  #better than null
confint(precip.water.dense, level=0.85, type="det") # 85% CI overlaps zero

# precip cat proceeds 

##### Temperature Hypothesis Group ####

temp.max <- occu( ~ max.temp ~ 1, occu.leca, starts = c(-1, 0, 0))
null.aicc - AICc(temp.max, k=2) #better than null
confint(temp.max, level=0.85, type="det") 
# 85% CI does not overlaps zero but is positive (no support for hypothesis)

# Temperature interactions

temp.water.dense.season <- occu( ~ max.temp * water_tank_density_0_9km
                                 ~ 1, occu.leca)       
null.aicc - AICc(temp.water.dense.season, k=2)  #worse than null

temp.canopy <- occu( ~ max.temp * canopy_cov ~ 1, occu.leca)    
null.aicc - AICc(temp.canopy, k=2) #worse than null

temp.humid <- occu( ~ max.temp * humid ~ 1, occu.leca)      
null.aicc - AICc(temp.humid, k=2) #worse than null

# none proceed

#### Wind Hypothesis Group ####

wind <- occu( ~ wind ~ 1, occu.leca, starts = c(-1, 0, 0))       
null.aicc - AICc(wind, k=2) #better than null
confint(wind, level=0.85, type="det") #85% CI does not overlap zero

# wind proceeds

##### Predator Activity Hypothesis Group ####

coy.hours <- occu( ~ coyote.active ~ 1, occu.leca, starts = c(-1, 0, 0)) 
null.aicc - AICc(coy.hours,k=2) #worse than null

coy.total <- occu( ~ coyote.count ~ 1, occu.leca, starts = c(-1, -1, 0)) 
null.aicc - AICc(coy.total,k=2) #worse than null

coy.bob.hours <- occu( ~ coy.bob.active ~ 1, occu.leca, starts = c(-1, 0, 0)) 
null.aicc - AICc(coy.bob.hours,k=2)  #worse than null

coy.bob.total <- occu( ~ coy.bob.count ~ 1, occu.leca, starts = c(-1, -1, 0)) 
null.aicc - AICc(coy.bob.total,k=2) #worse than null

all.pred.hours <- occu( ~ all.pred.active ~ 1, occu.leca, starts = c(-1, 0, 0)) 
null.aicc - AICc(all.pred.hours,k=2) #worse than null

all.pred.total <- occu( ~ all.pred.count ~ 1, occu.leca, starts = c(-1, -1, 0)) 
null.aicc - AICc(all.pred.total,k=2) #worse than null

bob.hours <- occu( ~ bobcat.active ~ 1, occu.leca, starts = c(-1, 0, 0)) 
null.aicc - AICc(bob.hours,k=2)  #worse than null

bob.total <- occu( ~ bobcat.count ~ 1, occu.leca, starts = c(-1, -1, 0)) 
null.aicc - AICc(bob.total,k=2) #worse than null

meso.hours <- occu( ~ meso.active ~ 1, occu.leca, starts = c(-1, -1, 0))        
null.aicc - AICc(meso.hours,k=2) #worse than null

meso.total <- occu( ~ meso.count ~ 1, occu.leca, starts = c(-1, -1, 0))        
null.aicc - AICc(meso.total,k=2) #worse than null

# none proceed

##### Livestock Activity Hypothesis Group####

cow.hours <- occu( ~ cow.active ~ 1, occu.leca, starts = c(-1, 0, 0))       
null.aicc - AICc(cow.hours, k=2) #worse than null

cow.total <- occu( ~ cow.count ~ 1, occu.leca, starts = c(-1, 0, 0))       
null.aicc - AICc(cow.total,k=2) #worse than null

stock.total <- occu( ~ livestock.count ~ 1, occu.leca, starts = c(-1, 0, 0))
null.aicc - AICc(stock.total,k=2) #worse than null

stock.hours <- occu( ~ livestock.active ~ 1, occu.leca, starts = c(-1, 0, 0))     
null.aicc - AICc(stock.hours,k=2) #worse than null

# livestock and water tank interactions 

stock.water.dense.season <- occu( ~ cow.active * water_tank_density_0_9km
                                  ~ 1, occu.leca)          
null.aicc - AICc(stock.water.dense.season,k=2) #worse than null

# none proceed

##### Biotic Community Type Hypothesis Group ####

ndvi <- occu( ~ NDVI_0_9km ~  1, occu.leca)
null.aicc - AICc(ndvi)  #better than null
confint(ndvi, level=0.85, type="det") # 85% CI does not overlap zero 

wood <- occu( ~ woodland_percent_0_9km ~ 1, occu.leca)          
null.aicc - AICc(wood, k=2) #better than null
confint(wood, level=0.85, type="det") # 85% CI does not overlap zero 

bio.com <- occu( ~ as.factor(biotic_com_2) ~ 1, occu.leca, starts = c(0, 1, 0))
null.aicc - AICc(bio.com, k=2)  #better than null
confint(bio.com, level=0.85, type="det") 
# 85% CI does not overlaps zero but no support for hypothesis

tree <- occu( ~  tree_density_5_70 ~ 1, occu.leca)
null.aicc - AICc(tree)  #better than null
confint(tree, level=0.85, type="det") # 85% CI does not overlap zero 

# NDVI proceeds

##### Vegetation Concealment Cover Hypothesis Group ####

obs.distance <- occu( ~ obs_dist ~ 1, occu.leca)       
null.aicc - AICc(obs.distance, k=2) #better than null
confint(obs.distance, level=0.85, type="det") # 85% CI does not overlap zero 

veg.cov.cam <- occu( ~ veg_cover_cam_under_1m 
                     ~ 1, occu.leca, starts = c(-1, -1, 0))       
null.aicc - AICc(veg.cov.cam, k=2)  #better than null
confint(veg.cov.cam, level=0.85, type="det") # 85% CI does not overlap zero 

detect.angle <- occu( ~ detection_angle ~ 1, occu.leca, starts = c(-1, 0, 0))  
null.aicc - AICc(detect.angle, k=2) #better than null
confint(detect.angle, level=0.85, type="det") 
# 85% CI does not overlaps zero but negative (no support for hypothesis)

max.trig.dist <- occu( ~ max_trig_dist ~ 1, occu.leca, starts = c(-1, 0, 0))   
null.aicc - AICc(max.trig.dist, k=2) #better than null
confint(max.trig.dist, level=0.85, type="det") # 85% CI does not overlap zero 

# veg cover proceeds

#### Sampling Period Hypothesis Group ####

### Sampling Period Length ###
samp.period <- occu( ~ sampling_period_length ~ 1, occu.leca)          
null.aicc - AICc(samp.period,k=2)  #better than null
confint(samp.period, level=0.85, type="det")  # 85% CI does not overlap zero  

### Camera Moved ###
cam.moved <- occu( ~ as.factor(cam_moved) ~ 1, occu.leca,starts = c(-1, 0, 0))   
null.aicc - AICc(cam.moved, k=2) #worse than null

## Camera Type ##
cam.type <- occu( ~ as.factor(jackrabbit.cam) ~ 1, occu.leca) 
null.aicc - AICc(cam.type, k=2) #worse than null

# sampling period proceeds

#################################################
# DETECTION MODEL SELECTION #
#################################################

# Variables that proceed (6)
  # precipitation 
  # NDVI
  # sampling period
  # veg cover
  # water tank dist 
  # wind

# Check correlations between detection variables 

# Read in site-level covariates
site.covs <- read.csv("site_covs.csv", nrows = 102, header = TRUE)

# site-level variables 
site.covs.cor <- site.covs %>% 
  select(sampling_period_length, veg_cover_cam_under_1m, NDVI_0_9km,
         dist_water_tank)

cor_site_covs <- cor(site.covs.cor, method='spearman') 
# none correlated above |0.7|

# Read in observation-level covariates
obs.covs <- readRDS("obsCovs.RData")

#  observation-level variables
obs_cor <-cor(data.frame(
  precip = as.vector(obs.covs$precip.cat),
  wind = as.vector(obs.covs$wind)
))
# none correlated above |0.7|

# Run this code to get correlation matrix in Excel 
#write.xlsx(obs_cor, file="Correlations jackrabbit.xlsx") 

# Null model
mod_null <- occu( ~ 1 ~ 1, occu.leca)

# single variable models 
mod1 <- occu( ~ NDVI_0_9km ~ 1, occu.leca)
mod2 <- occu( ~ sampling_period_length ~ 1, occu.leca)
mod3 <- occu( ~ veg_cover_cam_under_1m ~ 1, occu.leca)
mod4 <- occu( ~ dist_water_tank ~ 1, occu.leca)
mod5 <- occu( ~ as.factor(precip.cat) ~ 1, occu.leca)
mod6 <- occu( ~ wind ~ 1, occu.leca)

# 2-variable models 
mod7 <- occu( ~ NDVI_0_9km + veg_cover_cam_under_1m ~ 1, occu.leca)
mod8 <- occu( ~ NDVI_0_9km + dist_water_tank ~ 1, occu.leca)
mod9 <- occu( ~ NDVI_0_9km + as.factor(precip.cat) ~ 1, occu.leca)
mod10 <- occu( ~ sampling_period_length + veg_cover_cam_under_1m
               ~ 1, occu.leca)
mod11 <- occu( ~ sampling_period_length + dist_water_tank ~ 1, occu.leca)
mod12 <- occu( ~ sampling_period_length + as.factor(precip.cat) ~ 1, occu.leca)
mod13 <- occu( ~ veg_cover_cam_under_1m + dist_water_tank
               ~ 1, occu.leca)
mod14 <- occu( ~ veg_cover_cam_under_1m + as.factor(precip.cat)
               ~ 1, occu.leca)
mod15 <- occu( ~ dist_water_tank + as.factor(precip.cat) ~ 1, occu.leca)
mod16 <- occu( ~ NDVI_0_9km + sampling_period_length ~ 1, occu.leca)
mod17 <- occu( ~ NDVI_0_9km + wind ~ 1,  occu.leca)
mod18 <- occu( ~ wind + veg_cover_cam_under_1m ~ 1, occu.leca)
mod19 <- occu( ~ wind + dist_water_tank ~ 1, occu.leca)
mod20 <- occu( ~ wind + as.factor(precip.cat) ~ 1, occu.leca)
mod21 <- occu( ~ wind + sampling_period_length ~ 1, occu.leca)

# Model selection
cand.models <- list(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9, 
                    mod10, mod11, mod12, mod13, mod14, mod15, mod16, mod17,
                    mod18, mod19, mod20, mod21, mod_null)

modnames <- c("mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7", "mod8", 
              "mod9", "mod10", "mod11", "mod12", "mod13", "mod14", "mod15", 
              "mod16", "mod17", "mod18", "mod19","mod20", "mod21", "mod_null")

aicc_table <- aictab(cand.set = cand.models, modnames = modnames, sort = TRUE)
print(aicc_table)

# Candidate detection models are found in Table S2.5.

# Top model diagnostics 

summary(mod14)
AICc(mod_null) - AICc(mod14)

# Calculate the 85% confidence intervals for variables 
confint(mod14, type = "det", level = 0.85)    
#CIs for detection variables - none overlap zero

confint(mod14, type = "state", level = 0.85)  
#CIs for occupancy intercept - does not overlap zero

# Calculate multicollinearity

unmarked::vif(mod14, type = "det") # no colliniarity - all below 2

#################################################
# CREATE THE HABITAT SELECTION (OCCUPANCY) MODEL#
#################################################

## Black-tailed jackrabbit habitat selection hypotheses are listed in Table S2.2.

## Table S2.4 lists which occupancy hypothesis groups proceeded (i.e. performed 
  # better than the null occupancy model ψ(.), did not have uninformative 
  # parameters, and supported our original hypothesis).

## Candidate habitat selection (occupancy) models are found in Table S2.6.
  # Models were limited to 6 parameters and were derived from the best-supported
  # detection model (mod1 above).

occu.null.aicc <- AICc(mod14)

####Water Hypothesis Group####

water.dist.tank <- occu( ~ as.factor(precip.cat) + veg_cover_cam_under_1m
                         ~ dist_water_tank, occu.leca)      
occu.null.aicc - AICc(water.dist.tank, k=2) #worse than null 

water.tank.dense <- occu( ~ as.factor(precip.cat) + veg_cover_cam_under_1m
                          ~ water_tank_density_0_9km, occu.leca)       
occu.null.aicc - AICc(water.tank.dense) #worse than null

#none proceed

#### Predator Activity Hypothesis Group####

bobcat.coyote <- occu( ~ as.factor(precip.cat) + veg_cover_cam_under_1m
                       ~ pred_count_avg, occu.leca)
occu.null.aicc - AICc(bobcat.coyote, k=2) #worse than null 

coyote <- occu( ~ as.factor(precip.cat) + veg_cover_cam_under_1m
                ~ coy_count_avg, occu.leca)
occu.null.aicc - AICc(coyote, k=2) #worse than null 

# none proceed 

#### Livestock Activity Hypothesis Group####

cow <- occu( ~ as.factor(precip.cat) + veg_cover_cam_under_1m
             ~ cow_count_avg, occu.leca)
occu.null.aicc - AICc(cow, k=2) #better than null
confint(cow, level=0.85, type="state")
# 85% CI does not overlaps zero but positive (no support for hypothesis)

# none proceed 

#### Biotic Community Type Hypothesis Group ####

tree <- occu( ~  as.factor(precip.cat) + veg_cover_cam_under_1m
              ~  tree_density_5_70, occu.leca)
occu.null.aicc - AICc(tree) #better than null
confint(tree, level=0.85, type="state") #85% CI does not overlap zero

ndvi <- occu( ~ as.factor(precip.cat) + veg_cover_cam_under_1m
              ~ NDVI_0_9km, occu.leca)          
occu.null.aicc - AICc(ndvi, k=2) #better than null
confint(ndvi, level=0.85, type="state") #85% CI does not overlap zero

shrubland <- occu( ~  as.factor(precip.cat) + veg_cover_cam_under_1m
                   ~  shrub_yucca_density, occu.leca)
occu.null.aicc - AICc(shrubland)  #worse than null

grassland <- occu( ~  as.factor(precip.cat) + veg_cover_cam_under_1m
                   ~  herbaceous_cov, occu.leca)
occu.null.aicc - AICc(grassland) #worse than null 

wood <- occu(~  as.factor(precip.cat) + veg_cover_cam_under_1m
               ~ woodland_percent_0_9km, occu.leca)          
occu.null.aicc - AICc(wood, k=2) #better than null
confint(wood, level=0.85, type="state") #85% CI does not overlap zero

bio.com <- occu( ~ as.factor(precip.cat) + veg_cover_cam_under_1m
                 ~ as.factor(biotic_com_2), occu.leca)        
occu.null.aicc - AICc(bio.com, k=2)  #better than null
confint(bio.com, level=0.85, type="state") #85% CI does not overlap zero

# woodland percent proceeds

#### Vegetation Concealment Cover Hypothesis Group  ####

vert.cover <- occu( ~  as.factor(precip.cat) + veg_cover_cam_under_1m
                    ~  vertical_cover, occu.leca)
occu.null.aicc - AICc(vert.cover) #worse than null 

veg.cov.diff <- occu( ~  as.factor(precip.cat) + veg_cover_cam_under_1m
                      ~  veg_cover_diff, occu.leca)
occu.null.aicc - AICc(veg.cov.diff) #worse than null 

# none proceed 

##### Topography Hypothesis Group ####

topo.pos <- occu( ~  as.factor(precip.cat) + veg_cover_cam_under_1m
                  ~  topo_pos, occu.leca)
occu.null.aicc - AICc(topo.pos) #worse than null 

slope <- occu( ~  as.factor(precip.cat) + veg_cover_cam_under_1m
               ~  slope, occu.leca)
occu.null.aicc - AICc(slope)  #better than null
confint(slope, level=0.85, type="state") #85% CI does not overlap zero

elev <- occu( ~  as.factor(precip.cat) + veg_cover_cam_under_1m
              ~  elevation, occu.leca)
occu.null.aicc - AICc(elev)  #better than null
confint(elev, level=0.85, type="state") #85% CI does not overlap zero

# slope proceeds

##### Coolness of Site Hypothesis Group ####

aspect <- occu( ~  as.factor(precip.cat) + veg_cover_cam_under_1m
                ~  aspect, occu.leca)
occu.null.aicc - AICc(aspect) #worse than null 

heat.load <- occu( ~  as.factor(precip.cat) + veg_cover_cam_under_1m
                   ~  heat_load, occu.leca)
occu.null.aicc - AICc(heat.load) #worse than null

#none proceeds

#################################################
# OCCUPANCY MODEL SELECTION #
#################################################

# Variables that proceed (2)
  # slope
  # woodland percent

# Check correlations between occupancy variables 
site.covs.cor <- site.covs %>% 
  select(slope, woodland_percent_0_9km)

cor_site_covs <- cor(site.covs.cor, method='spearman') 
# none correlated above |0.7|

# Detection only model 
det.mod <- occu(  ~  as.factor(precip.cat) + 
                     veg_cover_cam_under_1m ~ 1, occu.leca)

# 1-variable occupancy models
mod1 <- occu( ~  as.factor(precip.cat) +
                 veg_cover_cam_under_1m
              ~  woodland_percent_0_9km, occu.leca)
mod2 <- occu( ~  as.factor(precip.cat) +
                 veg_cover_cam_under_1m
              ~  slope, occu.leca)

# 2-variable occupancy models
mod3 <- occu( ~  as.factor(precip.cat) +
                 veg_cover_cam_under_1m
              ~  woodland_percent_0_9km + slope, occu.leca)

# Model selection
cand.models <- list(mod1, mod2, mod3, det.mod)

modnames <- c("mod1", "mod2", "mod3", "det.mod")

aicc_table <- aictab(cand.set = cand.models, modnames = modnames, sort = TRUE)
print(aicc_table)

# Candidate occupancy models are found in Table S2.6.

# Top model diagnostics

summary(mod3)
AICc(det.mod) - AICc(mod3)

# The parameter estimates, standard errors, and 85% confidence intervals for
  #the best-supported model describing the probability of detection (p) and 
  #habitat selection (ψ) of black-tailed jackrabbits are listed in Table S2.14.

# Calculate the 85% confidence intervals for variables 
confint(mod3, type = "det", level = 0.85)    
#CIs for detection variables do not overlap zero

confint(mod3, type = "state", level = 0.85)  
#CIs for occupancy variables do not overlap zero

# Calculate multicollinearity
unmarked::vif(mod3, type = "det") # all below 2
unmarked::vif(mod3, type = "state") # all below 2

# Calculate mean daily detection probability with the top detection model

backTransform(linearComb(mod3,                           
                         coefficients= c(1,0,0),  
                         type = 'det'))     
#0.179 or 18%

### Occupancy probability

backTransform(linearComb(mod3,                
                         coefficients=c(1,0,0),   
                         type = 'state')) 
#0.87	OR 87%

### Calculating overall detection probability if K = 29
1-(1-0.179)^29

# 99% probability of detecting a black-tailed jackrabbit at least once given 
  #the minimum sampling period length at a site

### Calculating overall detection probability if K = 38
1-(1-0.179)^38

# 99% probability of detecting a black-tailed jackrabbit at least once given 
  #the average sampling period length at a site

#Predict detection probability for survey sites 

newdata <- data.frame(
  precip.cat = 0,
  veg_cover_cam_under_1m = seq(
    min(siteCovs(occu.leca)$veg_cover_cam_under_1m),
    max(siteCovs(occu.leca)$veg_cover_cam_under_1m),
    length = 102
  )
)

mod3 <- occu( ~ precip.cat + veg_cover_cam_under_1m
              ~ woodland_percent_0_9km + slope, occu.leca)

pred_det <- predict(mod3, type = "det", newdata = newdata)[
  ,c("Predicted","SE","lower","upper")
]

pred_det_df <- data.frame(Predicted = pred_det$Predicted,
                          StandardError = pred_det$SE,
                          lower = pred_det$lower,
                          upper = pred_det$upper,
                          site.covs)

#Predict occupancy probability for survey sites 

pred_occu <- predict(mod3,          
                     type = "state",                 
                     newdata = occu.leca@siteCovs)[c("Predicted",
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