# README:
# This script subsampled the original dataset to compare the effect on the detection probabilty estimates. 

# Libraries --------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(patchwork)
library(ggthemr)
ggthemr("fresh") 

# Effect of sample size ------------------------------------------------------------------

# subsample the data to test the effect on the estimated detection probability

#identify most sampled species
topSpecies <- data %>%
  group_by(Species) %>%
  summarise(nuRecs = length(Species)) %>%
  ungroup() %>%
  arrange(desc(nuRecs)) %>%
  slice_head(n=20) %>% # all above 1000 records
  pull("Species")

#draw up data frame to define the subsampling space (N is size of subsample)
sim_grid <- expand_grid(Species = topSpecies, N = seq(10,1000, by=10))

#function - modified version of function used in script 01
subsetsampleDistanceModel <- function(myspecies, data, samplesize = 50){
  
  #and just one species
  dataS <- data %>% 
    filter(Species==myspecies) %>%
    filter(type!="vinter") 
  
  infoS <- info %>%
    filter(type != "vinter")
  
  #get all transects for the obs
  allData <- inner_join(infoS, dataS)
  
  #make long
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
  
  #sample
  d <-   detectionData[sample(nrow(detectionData), samplesize, replace=FALSE),]
  
  hr.model <- ds(d, 
                 truncation=100,
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
                    N = nrow(d)))
  
}

subsampledESW<- list()

for(i in length(holder):nrow(sim_grid)){
  
  holder[[i]] <- subsetsampleDistanceModel(myspecies = sim_grid$Species[i], 
                                           data, 
                                           samplesize = sim_grid$N[i])
}

saveRDS(subsampledESW, file="outputs/subsampledESWs_top20_10.rds")

# Compare with original estimate using whole dataset and trait-based estimate
leave_one_out_results <- readRDS("outputs/leave_one_out_results.rds") # made in script 01
subsampledESW <- subsampledESW %>% bind_rows() %>% inner_join(.,leave_one_out_results, by="Species")

#compare the difference between
subsampledESW <- subsampledESW %>%
  mutate(obsDiff = abs(ESW - observed_ESW)/100,
         predDiff = abs(predicted_ESW - observed_ESW) /100) 

ggplot(subsampledESW) +
  geom_point(aes(x = N, y = obsDiff), alpha=0.2) +
  geom_smooth(aes(x = N, y = obsDiff),method="gam") +
  geom_line(aes(x=N, y=predDiff), linetype="dashed", color="red") +
  facet_wrap(~Species) + 
  xlab("Sample size") + ylab("Difference from best estimate") +
  theme(strip.text = element_text(size = 8))+
  ylim(0, 0.3)

ggsave(file = "plots/Figure_6.png", width = 6.75, height = 5)

# end --------------------------------------------------------------------------