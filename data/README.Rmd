Data files used in analyzing and plotting the effects of wind energy facilities on the habitat selection of terrestrial mammals in central New Mexico rangelands.

# Covariates

Site-level variables are contained in the file "site_covs.csv", these variables do not change over the duration the site was surveyed. Observation-level variables are in separate .csv files and are unique to every 24-hour sampling period. The "pre_model_code.R" code must be run before all species base models, wind models, and plot creation, as it sets up the variables for use in the occupancy model framework. 

# Base Models

Species' base models describe their habitat selection without accounting for the effects of the wind energy facilities. These models were created before adding the wind variables to the base model to determine whether accounting for the wind energy variables resulted in more explanatory power than the base model alone.

# Wind Models 

Each species' wind model was created by adding the 6 wind energy variables (wind farm interior, turbine visibility, distance to nearest wind turbine, density of wind turbines surrounding the site, distance to nearest turbine access road, and density of turbine access roads surrounding the site) to the base model to determine how the wind energy facilities affected the habitat selection of each terrestrial mammal species.
