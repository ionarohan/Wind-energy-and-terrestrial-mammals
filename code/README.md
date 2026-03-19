Code used for analyzing and plotting the effects of wind energy facilities on the habitat selection of terrestrial mammals in central New Mexico rangelands. Note that the "1_pre_model_code.R" script must be run before all species base models, wind models, and plot creation, as it sets up the variables for use in the occupancy model framework. 

# Base Models

Species' base models are numbered "2" and describe their habitat selection without accounting for the effects of the wind energy facilities. These models were created before adding the wind variables to the base model to determine whether accounting for the wind energy variables resulted in more explanatory power than the base model alone.

# Wind Models 

Species' wind models are numbered "3" with the letter corresponding to their base model. Each species' wind model was created by adding the 6 wind energy variables (wind farm interior, turbine visibility, distance to nearest wind turbine, density of wind turbines surrounding the site, distance to nearest turbine access road, and density of turbine access roads surrounding the site) to the base model to determine whether accounting for the wind variables resulted in more explanatory power than the base model alone.

# Plots 

The plots for the variables describing each species' probability of detection and habitat selection within their base models can be created using the script "2j_all species base plots.R". The plots showing the effects of the wind variables on the species' probability of habitat selection can be created using the script "3j_all species wind plots.R".
