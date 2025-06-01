# README:
# The script runs the distance models on the available data, including:
# null models to estimate constant mean detection probability
# covariates modelt to test the effect of forest cover and road/paths
# annual models to estimate detection probabilities for each year separately

### libraries ####

library(tidyverse)
library(lubridate)
library(Distance)
library(ggthemr)
ggthemr("fresh") 

# read in objects -----------------------------------------------------------------------------

data <- readRDS("data/data.rds") # data set of species observations (has to be requested from DOF)
info <- readRDS("data/info.rds") # data set with survey information (has to be requested from DOF)
environData <- readRDS("data/environData.rds") # dataset of environmental covariates for each survey location

# distance model functions --------------------------------------------------------------------

distanceModel <- function(myspecies, data, model="null"){
  
  #species data
  dataS <- data %>% 
    filter(Species==myspecies) %>%
    filter(type!="vinter") # remove winter data
  
  #survey info data
  infoS <- info %>%
    filter(type != "vinter") # remove winter data
  
  #get all transects for the obs
  allData <- inner_join(infoS, dataS)
  
  #make long format
  detectionData <- allData %>%
    dplyr::select(kvadratnr,Year, X.0, X.1, X.2) %>%
    pivot_longer(starts_with("X."), names_to = "distance", values_to = "nu") %>%
    dplyr::filter(!is.na(nu)) %>%
    dplyr::mutate(distance = recode(distance, X.0 = "12.5", X.1 = "37.5", X.2 = "75")) %>%
    dplyr::mutate(distance = as.numeric(as.character(distance))) %>%
    rename(size = nu) %>%
    filter(size>0) %>%
    filter(!is.na(size))
  
  #add on environmental data
  detectionData <- inner_join(detectionData, environData)
  
  if(model=="null"){ # standard distance model with no covariates
    
    hr.model <- ds(detectionData, truncation=100,
                   transect="line",
                   formula = ~1,
                   cutpoints = c(0,25,50,100),
                   key = "hn")
    
    temp <- summary(hr.model)
    
    return(data.frame(Species = myspecies,
                      P = predict(hr.model)$fitted[1],
                      P_se = predict(hr.model,  se.fit=TRUE)$se.fit[1],
                      ESW = predict(hr.model,  esw=TRUE)$fitted[1],
                      ESW_se = predict(hr.model, esw=TRUE, se.fit=TRUE)$se.fit[1],
                      AIC = AIC(hr.model)$AIC,
                      N = nrow(detectionData)))
    
  }else if(model=="covariates"){ #distance model including the effects of path/road cover and forest cover
    
    hr.model <- ds(detectionData, truncation=100,
                   transect="line",
                   formula = ~lines_pathroad + lines_forest,
                   cutpoints = c(0,25,50,100),
                   key = "hn")
    
    temp <- summary(hr.model)
    
    return(data.frame(Species = myspecies,
                      P = temp$ds$average.p,
                      P_se = temp$ds$average.p.se,
                      ESW = temp$ds$average.p*100,
                      ESW_se = temp$ds$average.p.se*100,
                      AIC = AIC(hr.model)$AIC,
                      N = nrow(detectionData),
                      Coef_pathroad = temp$ds$coeff$key.scale$estimate[2],
                      Coef_forest = temp$ds$coeff$key.scale$estimate[3],
                      Z_pathroad =  temp$ds$coeff$key.scale$estimate[2]/temp$ds$coeff$key.scale$se[2],
                      Z_forest =  temp$ds$coeff$key.scale$estimate[3]/temp$ds$coeff$key.scale$se[3]))
    
  }
  
}

# fit distance models -------------------------------------------------------------------------

## constant -----------------------------------------------------------------------------------

# null model with no covariates

distanceSp <- lapply(species, function(x){
  distanceModel(x, data, model="null")
}) %>%
  reduce(rbind)

saveRDS(distanceSp, file="outputs/distance_passerines_constant.rds")

#plot the detection probabilities (ESW = effectvive strip width)
distanceSp %>%
ggplot() + 
  geom_point(aes(x = fct_reorder(Species, ESW), y = ESW/100)) +
  xlab("") + ylab("Detection probability")+
  coord_flip()

ggsave("plots/Figure_S1.png", height = 8.5, width = 5)

## covariates ------------------------------------------------------------------------------

#model including the effects of paths/roads and forests

distanceSp_covs <- lapply_with_error(species, function(x){
  distanceModel(x, data, model="covariates")
}) %>%
  reduce(rbind)

saveRDS(distanceSp, file="outputs/distance_passerines_covariates.rds")

## annual models --------------------------------------------------------------------------

# fit null models to each year separately

distanceSp <- lapply(species, function(x){
  print(x)
  distanceModel(x, subset(data,Year==2014), model="null")
}) %>%
  reduce(rbind)

saveRDS(distanceSp, file="outputs/distance_passerines_constant_2014.rds")

distanceSp <- lapply(species, function(x){
  print(x)
  distanceModel(x, subset(data,Year==2015), model="null")
}) %>%
  reduce(rbind)

saveRDS(distanceSp, file="outputs/distance_passerines_constant_2015.rds")


distanceSp <- lapply(species, function(x){
  print(x)
  distanceModel(x, subset(data,Year==2016), model="null")
}) %>%
  reduce(rbind)

saveRDS(distanceSp, file="outputs/distance_passerines_constant_2016.rds")


distanceSp <- lapply(species, function(x){
  print(x)
  distanceModel(x, subset(data,Year==2017), model="null")
}) %>%
  reduce(rbind)

saveRDS(distanceSp, file="outputs/distance_passerines_constant_2017.rds")

### end ####
