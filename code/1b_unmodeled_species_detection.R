###########################################################
## Code for: Turbine visibility is a strong predictor of ##
## altered habitat selection by terrestrial mammals at a ##
######## wind energy facility in central New Mexico #######
###########################################################
### This script is used to calculate the null detection ###
### probabilities for the species that were not modeled ###
############## for the creation Table 2 ###################
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
  #is the homewd directory on line 25.

#Set home working directory
#e.g. homewd = "C:/Users/ionar/Desktop/R Repository/Wind-energy-and-terrestrial-mammals/"
homewd = "<insert your folder here and end with a forward slash>"

#Set wd to data folder on your local computer 
setwd(paste0(homewd, "data/"))

# Load packages 
library(unmarked)
library(tidyverse)

###########################################################
# SETUP CODE FOR BARBARY SHEEP #
###########################################################

# Read in edited .csv
detHist <- read.csv(file = "barbary_sheep_detection_hist.csv", row.names = 1)
# Read in site.covs.scaled
site.covs.scaled <- readRDS("site_covs_scaled.RData")
# Read in obsCovs.scaled
obsCovs.scaled <- readRDS("obsCovs_scaled.RData")

# Change from integer to numeric
detHist <- detHist %>%
  mutate(across(1:51, as.numeric))

# Create detection history 
occu.barb <- unmarkedFrameOccu(y=detHist, 
                               siteCovs = site.covs.scaled, 
                               obsCovs = obsCovs.scaled)

barb.null <- occu(~1~1, occu.barb, linkPsi="logit", starts = c(-1, -1), se=TRUE,
                 control = list(maxit = 10000))

# Null detection probability 
backTransform(barb.null['det'])

### Calculate the overall detection probability if K = 38
1-(1-0.000263)^38

# 1% probability of detecting a barbary sheep at least once given the 
# average sampling period length at a site and the null detection model 
# (see Table 2)

# Null occupancy probability 
backTransform(barb.null['state'])

###########################################################
# SETUP CODE FOR ELK #
###########################################################

# Read in edited .csv
detHist <- read.csv(file = "elk_detection_hist.csv", row.names = 1)

# Change from integer to numeric
detHist %>%
  mutate(across(1:51, as.numeric))

# Create detection history 
occu.ceca <- unmarkedFrameOccu(y=detHist, 
                               siteCovs = site.covs.scaled, 
                               obsCovs = obsCovs.scaled)

ceca.null <- occu(~1~1, occu.ceca , linkPsi="logit", se=TRUE,
                  control = list(maxit = 20000))

# Null detection probability 
backTransform(ceca.null['det'])

### Calculate the overall detection probability if K = 38
1-(1-0.0349)^38

# 74% probability of detecting an elk at least once given the 
# average sampling period length at a site and the null detection model 
# (see Table 2)

# Null occupancy probability 
backTransform(ceca.null['state'])

###########################################################
# SETUP CODE FOR COLLARED PECARY #
###########################################################

# Read in edited .csv
detHist <- read.csv(file = "javelina_detection_hist.csv", row.names = 1)

# Change from integer to numeric
detHist <- detHist %>%
  mutate(across(1:51, as.numeric))

# Create detection history 
occu.jav <- unmarkedFrameOccu(y=detHist, 
                               siteCovs = site.covs.scaled, 
                               obsCovs = obsCovs.scaled)

jav.null <- occu(~1~1, occu.jav, linkPsi="logit", starts = c(-1, -1), se=TRUE,
                  control = list(maxit = 10000))

# Null detection probability 
backTransform(jav.null['det'])

### Calculate the overall detection probability if K = 38
1-(1-0.0795)^38

# 96% probability of detecting a collared pecary at least once given the 
# average sampling period length at a site and the null detection model 
# (see Table 2)

# Null occupancy probability 
backTransform(jav.null['state'])

###########################################################
# SETUP CODE FOR BOBCAT #
###########################################################

# Read in edited .csv
detHist <- read.csv(file = "bobcat_detection_hist.csv", row.names = 1)

# Change from integer to numeric
detHist <- detHist %>%
  mutate(across(1:51, as.numeric))

# Create detection history 
occu.lyru <- unmarkedFrameOccu(y=detHist, 
                               siteCovs = site.covs.scaled, 
                               obsCovs = obsCovs.scaled)

#### Null model ####
lyru.null <- occu(~1~1, occu.lyru, linkPsi="logit", se=TRUE)

# Null detection probability 
backTransform(lyru.null['det'])

### Calculate the overall detection probability if K = 38
1-(1-0.0207)^38

# 55% probability of detecting a bobcat at least once given the 
# average sampling period length at a site and the null detection model 
# (see Table 2)

# Null occupancy probability 
backTransform(lyru.null['state'])

###########################################################
# SETUP CODE FOR RACCOON #
###########################################################

# Read in edited .csv
detHist <- read.csv(file = "raccoon_detection_hist.csv", row.names = 1)

# Change from integer to numeric
detHist <- detHist %>%
  mutate(across(1:51, as.numeric))

# Create detection history 
occu.rac <- unmarkedFrameOccu(y=detHist, 
                              siteCovs = site.covs.scaled, 
                              obsCovs = obsCovs.scaled)

rac.null <- occu(~1~1, occu.rac, linkPsi="logit", starts = c(-1, -1), se=TRUE,
                 control = list(maxit = 10000))

# Null detection probability 
backTransform(rac.null['det'])

### Calculate the overall detection probability if K = 38
1-(1-0.000287)^38

# 1% probability of detecting a raccoon at least once given the 
  # average sampling period length at a site and the null detection model 
  # (see Table 2)

# Null occupancy probability 
backTransform(rac.null['state'])

###########################################################
# SETUP CODE FOR AMERICAN HOG-NOSED SKUNK #
###########################################################

# Read in edited .csv
detHist <- read.csv(file = "hog_skunk_detection_hist.csv", row.names = 1)

# Change from integer to numeric
detHist <- detHist %>%
  mutate(across(1:51, as.numeric))

# Create detection history 
occu.hog <- unmarkedFrameOccu(y=detHist, 
                              siteCovs = site.covs.scaled, 
                              obsCovs = obsCovs.scaled)

hog.null <- occu(~1~1, occu.hog, linkPsi="logit", starts = c(-1, -1), se=TRUE,
                 control = list(maxit = 10000))

# Null detection probability 
backTransform(hog.null['det'])

### Calculate the overall detection probability if K = 38
1-(1-0.000256)^38

# 1% probability of detecting a hog-nosed skunk at least once given the 
# average sampling period length at a site and the null detection model 
# (see Table 2)

# Null occupancy probability 
backTransform(hog.null['state'])

###########################################################
# SETUP CODE FOR NA PORCUPINE #
###########################################################

# Read in edited .csv
detHist <- read.csv(file = "porcupine_detection_hist.csv", row.names = 1)

# Change from integer to numeric
detHist <- detHist %>%
  mutate(across(1:51, as.numeric))

# Create detection history 
occu.pine <- unmarkedFrameOccu(y=detHist, 
                              siteCovs = site.covs.scaled, 
                              obsCovs = obsCovs.scaled)

pine.null <- occu(~1~1, occu.pine, linkPsi="logit", starts = c(-1, -1), se=TRUE,
                 control = list(maxit = 10000))

# Null detection probability 
backTransform(pine.null['det']).

### Calculate the overall detection probability if K = 38
1-(1-0.0135)^38

# 40% probability of detecting a porcupine at least once given the 
# average sampling period length at a site and the null detection model 
# (see Table 2)

# Null occupancy probability 
backTransform(pine.null['state'])
