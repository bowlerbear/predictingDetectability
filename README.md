# Predicting species' detectabilities

Scripts are:

"01_DetectionProb_analysis.R" 
The script runs the distance models on the available data, including: 
- null models to estimate constant mean detection probability 
- covariates models to test the effect of forest cover and road/paths
- annual models to estimate detection probabilities for each year separately.

"02_Traits_lm_models.R" 
This script includes:
- multiple regression to explore the relationship between species traits and mean detection probability
- uses a leave-one-out approach with the model to test predictive ability
- estimates the error from using these predicted estimate over the original estimate
- same as above except for the effects of covariates on detection probability

"03_Traits_brt_models.R"
- boosted regression models to explore the relationship between species traits and mean detection probability
- uses a leave-one-out approach with the model to test predictive ability
- estimates the error from using these predicted estimate over the original estimate
- same as above except for the effects of covariates on detection probability

"04_Annual_analysis.R"
This script includes:
- Comparison of estimates of detection probabilities across years
- Analysis of annual trait associations
- Effects of using trait-based estimates on population size estimates

"05_Subsampled_models.R" 
This script subsampled the original dataset to compare the effect on the detection probabilty estimates. 

"HPC_Abundance_models.R"
Script for HPC to run brt models to estimate site-level relative abundance estimates

"HPC_Detection_boostraps.R"
Script for HPC to run models to bootstrap the detection probability estimates

