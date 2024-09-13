#######################################################################
#######################################################################
##     This script was created by Dr. Jen Cruz as part of            ##
##            the Applied Population Ecology Class                  ###
##                                                                   ##  
## Here we import our cleaned data for season 1 of our occurrence    ##
#  observations for Piute ground squirrels at the NCA and run a      ##
## closed population occupancy analysis. See Mackenzie et al. 2002   ##
## for details of the model. The occupancy model is hierarchical with #
# two components: (1) an ecological submodel linking occupancy to    ##
## environmental predictors at the site. (2) an observation submodel ##
## linking our detection probability to relevant predictors.         ##
##                                                                   ##
#######################################################################
##### Set up your workspace and load relevant packages -----------
# Clean your workspace to reset your R environment. #
rm( list = ls() )
# Check that you are in the right project folder
getwd()

# Install new packages from "CRAN" repository. # 
install.packages( "unmarked" ) #package for estimating occupancy, N-mixtures, 
#and some multinomial approaches for capture data
install.packages( "MuMIn") # package for model selection and evaluation
# load packages:

library( tidyverse )#includes dplyr, tidyr and ggplot2
library( unmarked ) 
library( MuMIn )
## end of package load ###############
###################################################################
#### Load or create data -----------------------------------------
# set directory where your data are:

datadir <- paste( getwd(), "/Data/", sep = "" )
# load our cleaned data

closeddf <- read.csv( file = paste( datadir, "closedf.csv", sep = ""),
                      header = TRUE )
#view
head( closeddf ); dim( closeddf ) 
#### End of data load -------------
####################################################################
##### Ready data for analysis --------------
# the unmarked function has several functions to make data inport #
# easy
# We need to define which predictors we will link to which responses #
# We expect detection to be influenced by observer effects, but it could also #
# be affected by amount of cover obstructing visibility (so potentially a #
# negative relationship with sagebrush). #
# We expect occupancy to be influenced by habitat (sagebrush and cheatgrass) #
# Why don't we include temperature in this model for one season?
# Answer: The Feb.minT and AprMay.maxT are constant throughout the season. Without any change to max and min temps, we would not expect that there would be a significant effect from these weather variables on detection.
#
# Let's define our unmarked dataframe:
# Start by defining which columns represent the response (observed occurrences)
umf <- unmarkedFrameOccu( y = as.matrix( closeddf[ ,c("pres.j1", "pres.j2", "pres.j3")]),
                          # Define predictors at the site level:
                          siteCovs = closeddf[ ,c("sagebrush", "cheatgrass")],
                          # Define predictors at the survey level as a list:
                          obsCovs = list( obsv = closeddf[ ,c("observer.j1", "observer.j2", "observer.j3")] ) ) 
#characters coverted to factors 
summary(umf)
#now scale ecological predictors:
sc <- apply( siteCovs(umf), MARGIN = 2, FUN = scale )
# We replace the predictors in our unmarked dataframe with the scaled values:
siteCovs( umf ) <- sc
# Why do we scale predictors?
# Answer: Scaling predictors allows them to be compared to each other as well as improve the efficiency of the model. Without scaling the predictors, comparisons cannot be made effectively. Scaling predictors can also make interpretation of the data easier. Convergence of data is more efficient when there isn't a wide range of values and loads faster. Scaling factors can also ensure that certain predictors are not over represented in the model. 
#
# View summary of unmarked dataframe:
summary( umf )
# What does it tell us?
#observations, covariates, observation-level covariates by technician.
### end data prep -----------
### Analyze data ------------------------------------------
# We are now ready to perform our analysis. Since the number of predictors #
# is reasonable for the sample size, and there were no issues with #
# correlation, we focus on a single full, additive model:
fm.closed <- occu( ~ 1 + obsv + sagebrush
                   ~ 1 + sagebrush + cheatgrass, data = umf )
# Note that we start with the observation submodel, linking it to the intercept # 
# and observer effect, obsv. We then define the ecological submodel as related #
# to sagebrush and cheatgrass. We end by defining the data to be used.

# View model results:
fm.closed

# We can also estimate confidence intervals for coefficients in #
# ecological submodel:
confint( fm.closed, type = "state" )
# Why do we call them coefficients and not predictors?
# Answer: These terms refer to different components of the model. Predictors are the independent variables that are used to explain/predict the variation in the response variable. Coefficients are the parameters of the model that represent the strength/direction of the relationship between each predictor and the response variable.  
#
# coefficients for detection submodel:
confint( fm.closed, type = 'det' )
#
# Based on the overlap of the 95% CIs for your predictor coefficients, #
# can you suggest which may be important to each of your responses? #
# Answer: Based on the CIs generated, tech2, tech3, and tech4 have values that do not include zero, suggesting that they are significant to the model. The intercept and sagebrush CIs include zero, suggesting that these predictors might not have a meaningful effect on the response variable. 
# 
#############end full model ###########
###### Model selection ---------------------------------------
# Indiscriminate model selection has become popular in recent years. #
# Although we do not believe this is a suitable approach here, we #
# demonstrate two approaches for running various reduced, additive models: #

# We start by manually running alternative models:
( fm.2 <- occu( ~ 1 + obsv + sagebrush  ~ 1 + sagebrush, data = umf ) )
( fm.3 <- occu( ~ 1 + obsv + sagebrush ~ 1 + cheatgrass, data = umf ) )
( fm.4 <- occu( ~ 1 + obsv + sagebrush ~ 1, data = umf ) )
( fm.5 <- occu( ~ 1 + obsv ~ 1 + sagebrush + cheatgrass, data = umf ) )
( fm.6 <- occu( ~ 1 + obsv ~ 1 + sagebrush , data = umf ) )
( fm.7 <- occu( ~ 1 + obsv ~ 1 + cheatgrass, data = umf ) )
( fm.8 <- occu( ~ 1 + obsv ~ 1, data = umf ) )
( fm.9 <- occu( ~ 1 + sagebrush ~ 1 + sagebrush + cheatgrass, data = umf ) )
( fm.10 <- occu( ~ 1 + sagebrush ~ 1 + sagebrush , data = umf ) )
( fm.11 <- occu( ~ 1 + sagebrush ~ 1 + cheatgrass, data = umf ) )
( fm.12 <- occu( ~ 1 + sagebrush ~ 1, data = umf ) )
( fm.13 <- occu( ~ 1 ~ 1 + sagebrush + cheatgrass, data = umf ) ) 
( fm.14 <- occu( ~ 1 ~ 1 + sagebrush , data = umf ) )
( fm.15 <- occu( ~ 1 ~ 1 + cheatgrass, data = umf ) )
( fm.16 <- occu( ~ 1 ~ 1, data = umf ) )
# Use unmarked function we create a list of model options:
fms <- fitList( 'psi(sagebrush + cheatgrass)p(obsv+sagebrush)' = fm.closed,
                'psi(sagebrush)p(obsv+sagebrush)' = fm.2,
                'psi(cheatgrass)p(obsv+sagebrush)' = fm.3,
                'psi(.)p(obsv+sagebrush)' = fm.4,
                'psi(sagebrush + cheatgrass)p(obsv)' = fm.5,
                'psi(sagebrush)p(obsv)' = fm.6,
                'psi(cheatgrass)p(obsv)' = fm.7,
                'psi(.)p(obsv)' = fm.8,
                'psi(sagebrush + cheatgrass)p(sagebrush)' = fm.9,
                'psi(sagebrush)p(sagebrush)' = fm.10,
                'psi(cheatgrass)p(sagebrush)' = fm.11,
                'psi(.)p(sagebrush)' = fm.12,
                'psi(sagebrush + cheatgrass)p(.)' = fm.13,
                'psi(sagebrush)p(.)' = fm.14,
                'psi(cheatgrass)p(.)' = fm.15,
                'psi(.)p(.)' = fm.16 )
#Note this uses the traditional (.) format to signify an intercept only model.
# We use unmarked function modSel() to compare models using AIC:
unmarked::modSel(fms )

# Alternatively, to run all possible model combinations automatically we can #
# use the dredge() function in the MuMIn package. This package allows you to #
# select alternative Information Criterion metrics including AIC, AICc, QAIC, BIC # 
modelList <- MuMIn::dredge( fm.closed, rank = 'AIC' )
#fixed terms are psi(Int) and p(Int)
#view model selection table:
modelList
# Which model(s) was/were the most supported? 
# Answer: Model 8 was the most supported based on the AIC score being the lowest (154.0) with the weight of 0.740 suggesting that this model has a 74% chance of being the best model of the group. Model 14 was the second had the second lowest AIC (158.0) but had a weight of 0.099, suggesting that it had less than a 10% chance of being the best model.  
#
# Does this change the inference from running the full model alone? How?
# Answer: Yes, model selection improves the level of confidence in the significant predictors and provides a more reliable inference. Using the full model alone could cause the model to become overfitted and lead to less reliable inferences. 
# 
# When is model selection a suitable approach?
# Answer: Model selection is a suitable approach when you want to compare/choose among several models to find the one that best explains the data or make predictions while avoiding overfitting. Situations where you have multiple candidate models, if you are unsure of the significance of the predictors, or when you are worried about the complexity of the model. 
#
# What would our estimates of occupancy be if we had not done any modeling?
# calculate naive occupancy by assigning a site as occupied if occurrence was #
# detected in any of the surveys, and as empty if occurrence was not detected #
# in any of the surveys:
y.naive <- ifelse( rowSums( closeddf[ ,c("pres.j1", "pres.j2", "pres.j3")])>0,1,0)

# What are the estimates of occupancy from our models:
# Calculate Best Unbiased Predictors of site occupancy from each model:
# Estimate conditional occupancy at each site:
re <- ranef( fm.closed )
# the use those to estimate occupancy with the bup() function:
y.est.fm.closed <-round( bup(re, stat="mean" )) # Posterior mean
# Repeat this process for other top model and the null:
y.est.fm.14 <-round( bup(ranef(fm.14), stat="mean" ) ) # Posterior mean
y.est.fm.16 <-round( bup(ranef(fm.16), stat="mean" ) ) # Posterior mean
# Compare results among them:
y.est.fm.closed - y.naive
y.est.fm.closed - y.est.fm.14
y.est.fm.closed - y.est.fm.16
#view together
data.frame( y.naive, y.est.fm.closed, y.est.fm.14, y.est.fm.16 )
# What do these results tell us about the importance of accounting for effects that impact detection?
# Answer: The results highlight the critical importance of accounting for effects that impact detection. The Naive detection models can underestimate true occupancy because the fail to account for factors like observer skills/training, environmental conditions, or site characteristics that influence detection probability. The full model (model 8) takes into account these effects and provides a more accurate estimate of detection/occupancy, generating a more reliable inference.
#
# What was the estimated mean occupancy while keeping #
# sagebrush and cheatgrass at their mean values:
backTransform( linearComb( fm.closed, coefficients = c(1,0,0) , 
                           type = "state" ) )
# Note we transform the occupancy response (defined as state by unmarked) back #
# from the logit scale. The ecological model has 1 intercept and two predictors.#
# The predictors are scaled so their mean is 0, the intercept is 1, thus: c(1,0,0).#
# What was our estimated occupancy?
# Answer: The estimated occupancy generated is 0.382, which aligns with the estimate provided in the table. Suggesting that based on the logistic transformation of the linear combination, the predicted occupancy is around 38.2% when sagebrush and cheatgrass are at their mean (zero after scaling).
#
# What about our mean probability of detection for each observer?
# We start with observer 1:
backTransform( linearComb( fm.closed, coefficients = c(1,0,0,0,0), type = "det" ) )
#observer 2:
backTransform( linearComb( fm.closed, coefficients = c(1,1,0,0,0), type = "det" ) )
#observer 3:
backTransform( linearComb( fm.closed, coefficients = c(1,0,1,0,0), type = "det" ) )
#observer 4:
backTransform( linearComb( fm.closed, coefficients = c(1,0,0,1,0), type = "det" ) )
#mean occupancy for obsv 1 at mean % sagebrush:
backTransform( linearComb( fm.closed, coefficients = c(1,0,0,0,1), type = "det" ) )

# What do these results tell us about what drives occupancy and detection of #
# Piute ground squirrels in 2007?
# Answer:These results highlight the need to account for both observer effects and environmental variables (sagebrush and cheatgrass) when assessing occupancy/detection. Different observers have varying levels of effectiveness and environmental factors can significantly influence the detection probabilities. Observer 1 had the highest detection probability (51.2%) and observer 4 had the lowest detection probability (3.2%). 
#

# end of analysis ######

############################################################################
################## Save your data and workspace ###################

# This time we want to save our workspace so that we have access to all #
# the objects that we created during our analyses. #
save.image( "OccAnalysisWorkspace.RData" )

# Why don't we want to re-run the analyses every time instead?
# Answer: Re-running the analyses can be important for validating results and answering different questions but it is more practical to run an initial comprehensive analysis and then make targeted updates as necessary. Computational costs of re-running models can become expensive. Having a baseline model helps in understanding the impact of new additional variables or significant changes and can act as a reference point. 
#

########## End of saving section ##################################

############# END OF SCRIPT #####################################