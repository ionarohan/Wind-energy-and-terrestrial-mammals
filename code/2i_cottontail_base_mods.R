###########################################################
####### Code for: base models for desert cottontails ######
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
library(AICcmodavg)

###########################################################
# SETUP CODE FOR DESERT COTTONTAIL OCCUPANCY MODELS #
###########################################################

# Read in edited .csv
detHist <- read.csv(file = "cottontail_detection_hist.csv", row.names = 1)
# Read in site.covs.scaled
site.covs.scaled <- readRDS("site_covs_scaled.RData")
# Read in obsCovs.scaled
obsCovs.scaled <- readRDS("obsCovs_scaled.RData")

# Change from integer to numeric
detHist %>%
  mutate(across(1:51, as.numeric))

# Create detection history 
occu.syau <- unmarkedFrameOccu(y=detHist, 
                               siteCovs = site.covs.scaled, 
                               obsCovs = obsCovs.scaled)

#### Null model ####
syau.null <- occu( ~1 ~ 1, occu.syau , linkPsi="logit", se=TRUE,
                  control = list(maxit = 10000))
null.aicc <- AICc(syau.null)

# Null detection probability 
backTransform(syau.null['det'])

# Null occupancy probability 
backTransform(syau.null['state'])

#################################################
# CREATE THE DETECTION MODEL #
#################################################

## Desert cottontail detection hypotheses are listed in Table S2.1.

## Table S2.3 lists which detection hypothesis groups proceeded (i.e. performed 
  # better than the null detection model p(.), did not have uninformative 
  # parameters, and supported our original hypothesis).

## Candidate detection models are found in Table S2.5.
  # Models were limited to 4 parameters total and were tested with no variables 
  # on the probability of habitat selection (ψ(.)).

##### Water Hypothesis Group ####

water.dist.tank <- occu( ~ dist_water_tank ~ 1, occu.syau, starts = c(-1, 0, 0))
null.aicc - AICc(water.dist.tank, k=2) #worse than null 

# none proceed 

##### Precipitation Hypothesis Group####

# Daily precipitation 

precip.cat <- occu( ~ as.factor(precip.cat) 
                    ~ 1, occu.syau, starts = c(-1, -1, 0))
null.aicc - AICc(precip.cat) #worse than null

precip.cm <- occu( ~ precip.cm ~ 1, occu.syau, starts = c(-1, -1, 0))
null.aicc - AICc(precip.cm) #worse than null

# Long-term precipitation

days.since.rain <- occu( ~ days.since.rain ~ 1, occu.syau, starts = c(-1, 0, 0))    
null.aicc - AICc(days.since.rain, k=2) #worse than null

rain.month <- occu( ~ rain.month ~ 1, occu.syau, starts = c(-1, 0, 0))    
null.aicc - AICc(rain.month, k=2) #worse than null

rain.week <- occu( ~ rain.week ~ 1, occu.syau, starts = c(-1, 0, 0))    
null.aicc - AICc(rain.week, k=2) #worse than null

# Water source and long-term precipitation interaction

precip.water.dist <- occu( ~ dist_water_tank * rain.week ~ 1, occu.syau)
null.aicc - AICc(precip.water.dist, k=2) #worse than null

# none proceed

##### Temperature Hypothesis Group ####

test.temp.max <- occu( ~ max.temp ~ 1, occu.syau, starts = c(-1, 0, 0))
null.aicc - AICc(test.temp.max, k=2) #worse than null

# Temperature interactions

temp.water.dist <- occu( ~ max.temp * dist_water_tank ~ 1, occu.syau)
null.aicc - AICc(temp.water.dist, k=2) #worse than null

temp.canopy <- occu( ~ max.temp * canopy_cov ~ 1, occu.syau)      
null.aicc - AICc(temp.canopy, k=2) #better than null
confint(temp.canopy, level=0.85, type="det")  #85% CI does not overlap zero

temp.humid <- occu( ~ max.temp * humid ~ 1, occu.syau)      
null.aicc - AICc(temp.humid, k=2) #worse than null

# canopy cover interaction proceeds

#### Wind Hypothesis Group ####

wind <- occu( ~ wind ~ 1, occu.syau, starts = c(-1, 0, 0))    
null.aicc - AICc(wind, k=2) #worse than null

# none proceeds

##### Predator Activity Hypothesis Group ####

coy.hours <- occu( ~ coyote.active ~ 1, occu.syau, starts = c(-1, 0, 0)) 
null.aicc - AICc(coy.hours,k=2) #better than null
confint(coy.hours, level=0.85, type="det") 
# 85% CI does not overlap zero but is positive (no support for hypothesis)

coy.total <- occu( ~ coyote.count ~ 1, occu.syau, starts = c(-1, -1, 0)) 
null.aicc - AICc(coy.total,k=2) #better than null
confint(coy.total, level=0.85, type="det") 
# 85% CI does not overlap zero but is positive (no support for hypothesis)

coy.bob.hours <- occu( ~ coy.bob.active ~ 1, occu.syau, starts = c(-1, -1, 0)) 
null.aicc - AICc(coy.bob.hours,k=2) #better than null
confint(coy.bob.hours, level=0.85, type="det") 
# 85% CI does not overlap zero but is positive (no support for hypothesis)

coy.bob.total <- occu( ~ coy.bob.count ~ 1, occu.syau, starts = c(-1, -1, 0)) 
null.aicc - AICc(coy.bob.total,k=2) #better than null
confint(coy.bob.total, level=0.85, type="det") 
# 85% CI does not overlap zero but is positive (no support for hypothesis)

all.pred.hours <- occu( ~ all.pred.active ~ 1, occu.syau, starts = c(-1, -1, 0))
null.aicc - AICc(all.pred.hours,k=2) #better than null
confint(all.pred.hours, level=0.85, type="det") 
# 85% CI does not overlap zero but is positive (no support for hypothesis)

all.pred.total <- occu( ~ all.pred.count ~ 1, occu.syau, starts = c(-1, -1, 0)) 
null.aicc - AICc(all.pred.total,k=2) #better than null
confint(all.pred.total, level=0.85, type="det")
# 85% CI does not overlap zero but is positive (no support for hypothesis)

bob.hours <- occu( ~ bobcat.active ~ 1, occu.syau, starts = c(-1, -1, 0)) 
null.aicc - AICc(bob.hours,k=2) #worse than null

bob.total <- occu( ~ bobcat.count ~ 1, occu.syau, starts = c(-1, -1, 0)) 
null.aicc - AICc(bob.total,k=2) #worse than null

meso.hours <- occu( ~ meso.active ~ 1, occu.syau,  starts = c(-1, -1, 0))        
null.aicc - AICc(meso.hours,k=2) #worse than null

meso.total <- occu( ~ meso.count ~ 1, occu.syau,  starts = c(-1, -1, 0))        
null.aicc - AICc(meso.total,k=2) #worse than null

# none proceed

##### Livestock Activity Hypothesis Group####

cow.hours <- occu( ~ cow.active ~ 1, occu.syau, starts = c(-1, 0, 0))       
null.aicc - AICc(cow.hours, k=2)  #worse than null

cow.total <- occu( ~ cow.count ~ 1, occu.syau, starts = c(-1, 0, 0))       
null.aicc - AICc(cow.total,k=2) #better than null
confint(cow.total, level=0.85, type="det") # 85% CI does not overlap zero

stock.total <- occu( ~ livestock.count ~ 1, occu.syau, starts = c(-1, 0, 0))
null.aicc - AICc(stock.total,k=2) #better than null
confint(stock.total, level=0.85, type="det") # 85% CI does not overlap zero

stock.hours <- occu( ~ livestock.active ~ 1, occu.syau,  starts = c(-1, 0, 0))     
null.aicc - AICc(stock.hours,k=2) #worse than null

# Potential interactions

stock.pred <- occu( ~ cow.count * coyote.active ~ 1, occu.syau)          
null.aicc - AICc(stock.pred,k=2) #better than null
confint(stock.pred, level=0.85, type="det") # 85% CI overlaps zero

stock.water.dist <- occu( ~ cow.count * dist_water_tank ~ 1, occu.syau)          
null.aicc - AICc(stock.water.dist,k=2) #worse than null

# cow total proceeds

##### Biotic Community Type Hypothesis Group ####

ndvi <- occu( ~ NDVI_0_1km ~1 , occu.syau)          
null.aicc - AICc(ndvi, k=2) #worse than null

wood <- occu( ~ woodland_percent_0_1km ~ 1, occu.syau)          
null.aicc - AICc(wood, k=2) #worse than null

bio.com <- occu( ~ as.factor(biotic_com_2) ~ 1, occu.syau, starts = c(0, -3, 1))
null.aicc - AICc(bio.com, k=2)  #better than null
confint(bio.com, level=0.85, type="det") # 85% CI does not overlap zero

tree <- occu( ~  tree_density_5_70 ~ 1, occu.syau)
null.aicc - AICc(tree)  #worse than null

#bio com proceeds

##### Vegetation Concealment Cover Hypothesis Group ####

obs.distance <- occu( ~ obs_dist ~ 1, occu.syau)       
null.aicc - AICc(obs.distance, k=2) #worse than null

veg.cover.cam <- occu( ~ veg_cover_cam_under_1m ~ 1, occu.syau)   
null.aicc - AICc(veg.cover.cam, k=2) #better than null
confint(veg.cover.cam, level=0.85, type="det") # 85% CI does not overlap zero

detect.angle <- occu( ~ detection_angle ~ 1, occu.syau, starts = c(-1, 0, 0))  
null.aicc - AICc(detect.angle, k=2) #worse than null

max.trig.dist <- occu( ~ max_trig_dist ~ 1, occu.syau, starts = c(-1, -1, 0))   
null.aicc - AICc(max.trig.dist, k=2) #worse than null

# veg cover under 1 M proceeds

#### Sampling Period Hypothesis Group ####

### Sampling Period Length ###
samp.period <- occu( ~ sampling_period_length
                     ~ 1, occu.syau, starts = c(-1, -1, 0))         
null.aicc - AICc(samp.period,k=2) #worse than null

### Camera Moved ###
cam.moved <- occu( ~ as.factor(cam_moved) ~ 1, occu.syau)   
null.aicc - AICc(cam.moved, k=2) #better than null
confint(cam.moved, level=0.85, type="det") 
# 85% CI does not overlap zero but is positive (no support for hypothesis)

# none proceed

#################################################
# DETECTION MODEL SELECTION #
#################################################

# Parameters that proceed
  # veg cover under 1 m
  # canopy cover * max.temp
  # biotic community 
  # cow count

# Check correlations between detection variables 

# Read in site-level covariates
site.covs <- read.csv("site_covs.csv", nrows = 102, header = TRUE)

# site-level variable correlations
site.covs.cor <- site.covs %>%
  select(veg_cover_cam_under_1m, canopy_cov, biotic_com_2)

cor_site_covs <- cor(site.covs.cor, method='spearman') 

# Read in observation-level covariates
obs.covs <- readRDS("obsCovs.RData")

#  observation-level variables
obs_cor <-cor(data.frame(
  temp = as.vector(obs.covs$max.temp),
  cow = as.vector(obs.covs$cow.count)
))
# none correlated above |0.7|

# Run this code to get correlation matrix in Excel 
  # write.xlsx(obs_cor, file="Correlations cottontail.xlsx") 

#Null model
syau.null <- occu( ~ 1 ~ 1, occu.syau)

# Single effect model
mod1 <- occu( ~ veg_cover_cam_under_1m ~ 1, occu.syau)
mod2 <- occu( ~ cow.count ~ 1, occu.syau)
mod3 <- occu( ~ as.factor(biotic_com_2) ~ 1, occu.syau)
mod4 <- occu( ~ canopy_cov * max.temp ~ 1, occu.syau)

# multiple effect mods
mod5 <- occu( ~ veg_cover_cam_under_1m + cow.count
              ~ 1, occu.syau)
mod6 <- occu( ~ veg_cover_cam_under_1m + as.factor(biotic_com_2)
              ~ 1, occu.syau)
mod7 <- occu( ~ cow.count + as.factor(biotic_com_2)
              ~ 1, occu.syau)

# Model selection
cand.models <- list(mod1, mod2, mod3, mod4, mod5, mod6, mod7, syau.null)

modnames <- c("mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7", 
              "syau.null")

aicc_table <- aictab(cand.set = cand.models, modnames = modnames, sort = TRUE)
print(aicc_table)

# Candidate detection models are found in Table S2.5.

# Top model diagnostics 
AICc(syau.null) - AICc(mod6)
summary(mod6)

# Calculate multicollinearity
unmarked::vif(mod6, type = "det") # all below 2

# Calculate the 85% confidence intervals for variables 
confint(mod6, type = "det", level = 0.85)    
#CIs for detection variables - none overlap zero

confint(mod6, type = "state", level = 0.85)  
#CIs for occupancy intercept - does not overlap zero


#################################################
# CREATE THE HABITAT SELECTION (OCCUPANCY) MODEL#
#################################################

## Desert cottontail habitat selection hypotheses are listed in Table S2.2.

## Table S2.4 lists which occupancy hypothesis groups proceeded (i.e. performed 
  # better than the null occupancy model ψ(.), did not have uninformative 
  # parameters, and supported our original hypothesis).

## Candidate habitat selection (occupancy) models are found in Table S2.6.
  # Models were limited to 6 parameters and were derived from the best-supported
  # detection model (mod1 above).

occu.null.aicc <- AICc(mod6)

####Water Hypothesis Group####

water.dist.tank <- occu(  ~ veg_cover_cam_under_1m + as.factor(biotic_com_2) 
                          ~ dist_water_tank, occu.syau)      
occu.null.aicc - AICc(water.dist.tank, k=2) #worse than null 

#none proceed

#### Predator Activity Hypothesis Group####

bobcat.coyote <- occu( ~ veg_cover_cam_under_1m + as.factor(biotic_com_2)     
                       ~ pred_count_avg, occu.syau)
occu.null.aicc - AICc(bobcat.coyote, k=2) #worse than null 

coyote <- occu( ~ veg_cover_cam_under_1m + as.factor(biotic_com_2) 
                ~ coy_count_avg, occu.syau)
occu.null.aicc - AICc(coyote, k=2) #worse than null 

# none proceed 

#### Livestock Activity Hypothesis Group####

cow <- occu( ~ veg_cover_cam_under_1m + as.factor(biotic_com_2)             
             ~ cow_count_avg, occu.syau)
occu.null.aicc - AICc(cow, k=2) #worse than null 

# none proceed 

#### Biotic Community Type Hypothesis Group ####

tree <- occu( ~ veg_cover_cam_under_1m + as.factor(biotic_com_2)       
              ~ tree_density_5_70, occu.syau)
occu.null.aicc - AICc(tree) #better than null
confint(tree, level=0.85, type="state") #85% CI does not overlap zero

ndvi <- occu( ~ veg_cover_cam_under_1m + as.factor(biotic_com_2)      
              ~ NDVI_0_1km, occu.syau)
occu.null.aicc - AICc(ndvi) #better than null
confint(ndvi, level=0.85, type="state") #85% CI does not overlap zero

shrubland <- occu( ~ veg_cover_cam_under_1m + as.factor(biotic_com_2)      
                   ~ shrub_yucca_density, occu.syau)
occu.null.aicc - AICc(shrubland)  #worse than null 

grassland <- occu( ~ veg_cover_cam_under_1m + as.factor(biotic_com_2)    
                   ~ herbaceous_cov, occu.syau)
occu.null.aicc - AICc(grassland) #better than null
confint(grassland, level=0.85, type="state") #85% CI does not overlap zero

wood <- occu( ~ veg_cover_cam_under_1m + as.factor(biotic_com_2)  
              ~ woodland_percent_0_1km, occu.syau)          
occu.null.aicc - AICc(wood, k=2) #better than null
confint(wood, level=0.85, type="state") #85% CI does not overlap zero

bio.com <- occu( ~ veg_cover_cam_under_1m + as.factor(biotic_com_2)
                 ~ as.factor(biotic_com_2), occu.syau)      
occu.null.aicc - AICc(bio.com, k=2)  #better than null
confint(bio.com, level=0.85, type="state") #85% CI does not overlap zero

# NDVI proceeds

#### Vegetation Concealment Cover Hypothesis Group  ####

vert.cover <- occu( ~ veg_cover_cam_under_1m + as.factor(biotic_com_2)     
                    ~ vertical_cover, occu.syau)
occu.null.aicc - AICc(vert.cover) #worse than null 

veg.cover.diff <- occu( ~ veg_cover_cam_under_1m + as.factor(biotic_com_2)      
                        ~ veg_cover_diff, occu.syau)
occu.null.aicc - AICc(veg.cover.diff) #worse than null 

# none proceed 

##### Topography Hypothesis Group ####

topo.pos <- occu( ~ veg_cover_cam_under_1m + as.factor(biotic_com_2)       
                  ~ topo_pos, occu.syau)
occu.null.aicc - AICc(topo.pos) #better than null
confint(topo.pos, level=0.85, type="state") 
# 85% CI does not overlap zero but is positive (no support for hypothesis)

slope <- occu( ~ veg_cover_cam_under_1m + as.factor(biotic_com_2)        
               ~  slope, occu.syau)
occu.null.aicc - AICc(slope) #worse than null 

elev <- occu( ~ veg_cover_cam_under_1m + as.factor(biotic_com_2)           
              ~ elevation, occu.syau)
occu.null.aicc - AICc(elev) #worse than null 

# none proceeds

##### Coolness of Site Hypothesis Group ####

aspect <- occu( ~ veg_cover_cam_under_1m + as.factor(biotic_com_2)       
                ~ aspect, occu.syau)
occu.null.aicc - AICc(aspect) #worse than null 

heat.load <- occu( ~ veg_cover_cam_under_1m + as.factor(biotic_com_2)      
                   ~  heat_load, occu.syau)
occu.null.aicc - AICc(heat.load) #worse than null

#none proceeds

#################################################
# OCCUPANCY MODEL SELECTION #
#################################################

# Variables that proceed (1)
  # NDVI

# Detection only model 
det.mod <- occu( ~ veg_cover_cam_under_1m + as.factor(biotic_com_2) 
                 ~ 1, occu.syau)

#single variable mods
mod1 <- occu( ~  veg_cover_cam_under_1m + as.factor(biotic_com_2)        
              ~  NDVI_0_1km, occu.syau)

# Model selection
cand.models <- list(mod1, det.mod)

modnames <- c("mod1", "det.mod")

aicc_table <- aictab(cand.set = cand.models, modnames = modnames, sort = TRUE)
print(aicc_table)

# Candidate occupancy models are found in Table S2.6.

# Top model diagnostics

summary(mod1)
AICc(det.mod)-AICc(mod1)

# The parameter estimates, standard errors, and 85% confidence intervals for
  #the best-supported model describing the probability of detection (p) and 
  #habitat selection (ψ) of desert cottontail are listed in Table S2.15.

### Calculating 85% confidence intervals for variables 
confint(mod1, type = "det", level = 0.85)    
#CIs for detection variables do not overlap zero

confint(mod1, type = "state", level = 0.85)  
#CIs for occupancy variables do not overlap zero

# Calculate mean daily detection probability with the top detection model

backTransform(linearComb(mod1,                           
                         coefficients= c(1,0,0),  
                         type = 'det'))     
# 0.122 OR 12%


### Occupancy probability

backTransform(linearComb(mod1,                
                         coefficients=c(1,0),   
                         type = 'state')) 
#0.381 OR 38%

### Calculating overall detection probability if K = 29
1-(1-0.122)^29

# 98% probability of detecting a desert cottontail at least once given the 
  #minimum sampling period length at a site

### Calculating overall detection probability if K = 38
1-(1-0.122)^38

# 99% probability of detecting a desert cottontail at least once given the 
  #minimum sampling period length at a site

# predict detection 
p_full <- predict(mod1, type = "det")

# Get number of sites and occasions
n_sites <- numSites(occu.syau)
n_occasions <- ncol(getY(occu.syau))
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

pred_occu <- predict(mod1,          
                     type = "state",                 
                     newdata = occu.syau@siteCovs)[c("Predicted",
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