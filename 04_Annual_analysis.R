# README:
# This script includes:
# (1) Comparison of estimates of detection probabilities across years
# (2) Analysis of annual trait associations
# (3) Effects of using trait-based estimates on population size estimates

# libraries --------------------------------------------------------------------

library(tidyverse)
library(broom)
library(lubridate)
library(patchwork)
library(ggExtra)
library(corrplot)
library(DescTools)
library(ggrepel)
library(modelsummary)
library(ggthemr)
ggthemr("fresh") 

# Get annual detection probabilities  --------------------------------------------------

#see script 01 for where these are made
nullDF2014 <- readRDS("outputs/distance_passerines_constant_2014.rds") %>% add_column(Year=2014)
nullDF2015 <- readRDS("outputs/distance_passerines_constant_2015.rds") %>% add_column(Year=2015)
nullDF2016 <- readRDS("outputs/distance_passerines_constant_2016.rds") %>% add_column(Year=2016)
nullDF2017 <- readRDS("outputs/distance_passerines_constant_2017.rds") %>% add_column(Year=2017)

# Combine all -------------------------------------------------------------------

nullDFannual <- bind_rows(nullDF2014,
                          nullDF2015,
                          nullDF2016,
                          nullDF2017) %>%
                dplyr::select(ESW,Species,Year) %>%
                pivot_wider(names_from = Year,
                            values_from = ESW) %>%
                janitor::clean_names() %>%
                filter(complete.cases(.))

# Correlation ------------------------------------------------------------------------

#how correlated are species detecton probabilities across years?
pairs(nullDFannual[,-1])
cor(nullDFannual[,-1])

#the concordance correlation coefficients
CCC(nullDFannual$x2014, nullDFannual$x2015)$rho.c
#est    lwr.ci    upr.ci
#1 0.8003876 0.6986817 0.8703836
CCC(nullDFannual$x2015, nullDFannual$x2016)$rho.c
#est    lwr.ci    upr.ci
#1 0.8401798 0.7570643 0.8965335
CCC(nullDFannual$x2016, nullDFannual$x2017)$rho.c
#est    lwr.ci    upr.ci
#1 0.8401798 0.7570643 0.8965335

# Plotting --------------------------------------------------------------------------

#swatch()
panela <- ggplot(nullDFannual) +
  geom_point(aes(x = x2014/100, y = x2015/100, colour = "2014 vs 2015")) +
  geom_point(aes(x = x2015/100, y = x2016/100, colour = "2015 vs 2016")) +
  geom_point(aes(x = x2016/100, y = x2017/100, colour = "2016 vs 2017")) +
  xlab("Previous detection probabilty") +
  ylab("Detection probability") +
  geom_abline(linetype = "dashed") +
  scale_colour_manual(
    values = c(
      "2014 vs 2015" = swatch()[1],
      "2015 vs 2016" = swatch()[5],
      "2016 vs 2017" = swatch()[9]
    )
  ) +
  labs(colour = "Year Comparison") + 
  theme(
    legend.position.inside = c(0.1, 0.9))  
panela
# this eventually is Fig. 5a - combined with other plots below

# Trait associations --------------------------------------------------

traits <- readRDS("traits/traits.rds")

#some tweaks
table(traits$Habitat)
traits$Habitat[traits$Habitat=="Woodland"] <- "Forest"
#for black redstart, change to:
traits$Habitat[traits$Habitat=="Rock"] <- "Human Modified"

traitsDF <- nullDFannual %>%
  rename(Species = species) %>%
  inner_join(., traits, by="Species") 

lm2014 <- lm(x2014 ~ log(Mass) + ForStrat.ground + Habitat + 
            Trophic.Level, data=traitsDF)
summary(lm2014)
anova(lm2014)
#Response: x2014
#Df  Sum Sq Mean Sq F value    Pr(>F)    
#log(Mass)        1  6375.9  6375.9 33.4843 2.547e-07 ***
#  ForStrat.ground  1  1171.8  1171.8  6.1542   0.01584 *  
#  Habitat          4   902.9   225.7  1.1855   0.32605    
#  Trophic.Level    2   708.4   354.2  1.8601   0.16423


lm2015 <- lm(x2015 ~ log(Mass) + ForStrat.ground + Habitat + 
             Trophic.Level, data=traitsDF)
summary(lm2015)
anova(lm2015)
#Response: x2015
#Df Sum Sq Mean Sq F value    Pr(>F)    
#log(Mass)        1 8326.9  8326.9 52.0076 9.244e-10 ***
#  ForStrat.ground  1 1388.1  1388.1  8.6697   0.00455 ** 
#  Habitat          4  875.8   219.0  1.3676   0.25550    
#Trophic.Level    2 1830.2   915.1  5.7156   0.00527 **
  

lm2016 <- lm(x2016 ~ log(Mass) + ForStrat.ground + Habitat + 
             Trophic.Level, data=traitsDF)
summary(lm2016)
anova(lm2016)
#Response: x2016
#Df  Sum Sq Mean Sq F value    Pr(>F)    
#log(Mass)        1 10843.9 10843.9 64.8523 3.201e-11 ***
#  ForStrat.ground  1   685.9   685.9  4.1023  0.047136 *  
#  Habitat          4  1033.9   258.5  1.5458  0.200080    
#Trophic.Level    2  1810.0   905.0  5.4125  0.006815 ** 

lm2017 <- lm(x2017 ~ log(Mass) + ForStrat.ground + Habitat + 
             Trophic.Level, data=traitsDF)
summary(lm2017)
anova(lm2017)
#Response: x2017
#Df Sum Sq Mean Sq F value    Pr(>F)    
#Df Sum Sq Mean Sq F value    Pr(>F)    
#log(Mass)        1 6785.5  6785.5 42.5795 1.416e-08 ***
#  ForStrat.ground  1 1123.1  1123.1  7.0475   0.01007 *  
#  Habitat          4 1875.4   468.9  2.9421   0.02719 *  
#  Trophic.Level    2 1521.8   760.9  4.7747   0.01179 *  

modelsummary(lm2014, statistic = c("conf.low", "conf.high"), 
             shape = term ~ model + statistic,
             output="plots/Linear_Model_Ouput_2014.docx")

modelsummary(lm2015, statistic = c("conf.low", "conf.high"), 
             shape = term ~ model + statistic,
             output="plots/Linear_Model_Ouput_2015.docx")

modelsummary(lm2016, statistic = c("conf.low", "conf.high"), 
             shape = term ~ model + statistic,
             output="plots/Linear_Model_Ouput_2016.docx")

modelsummary(lm2017, statistic = c("conf.low", "conf.high"), 
             shape = term ~ model + statistic,
             output="plots/Linear_Model_Ouput_2017.docx")

#outputs included in Table S4

# Prediction -----------------------------------------------------------------------

#lets test how well a trait-based estimate using last years data can be used to predict this year
# of detection probability

## 2014 to 2015 ------------------------------------------------------------

#we'll use the 2014 data to make trait predictions 
lm1 <- lm(x2014 ~ log(Mass) + ForStrat.ground + Habitat + 
            Trophic.Level, data=traitsDF)
traitsDF$preds2014 <- predict(lm1)

#we'll then compare it to what was actually observed in 2015
summary(abs(traitsDF$x2015 - traitsDF$x2014)) #compare to past estimate
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.000   1.781   4.843   7.447   9.089  38.505
summary(abs(traitsDF$x2015 - traitsDF$preds2014))#compare to trait-based estimate
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.2056  3.7681  6.4069  9.2844 13.1165 38.4673

CCC(traitsDF$x2015, traitsDF$x2014)$rho.c # stronger correlation with past estimate
#est    lwr.ci    upr.ci
#1 0.8003876 0.6986817 0.8703836
CCC(traitsDF$x2015, traitsDF$preds2014)$rho.c
#est    lwr.ci    upr.ci
#1 0.6620723 0.5381027 0.7580043

#package
errors1415 <- data.frame(error = c(CCC(traitsDF$x2015, traitsDF$x2014)$blalt$delta,
                               CCC(traitsDF$x2015, traitsDF$preds2014)$blalt$delta),
                     type = c(rep("past estimate",71), rep("trait-based estimate",71)),
                     species = rep(traitsDF$Species, 2))

## 2015 to 2016 ------------------------------------------------------------

#we'll use the 2015 data to make trait predictions 
lm1 <- lm(x2015 ~ log(Mass) + ForStrat.ground + Habitat + 
            Trophic.Level, data=traitsDF)
traitsDF$preds2015 <- predict(lm1)

#we'll then compare it to what was actually observed in 2015
summary(abs(traitsDF$x2016 - traitsDF$x2015)) #compare to past estimate
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.000   1.643   4.370   6.631   7.829  42.790 

summary(abs(traitsDF$x2016 - traitsDF$preds2015)) #compare to trait-based estimate
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.3249  3.6650  7.2086  9.7670 14.0649 31.3305 

CCC(traitsDF$x2016, traitsDF$x2015)$rho.c # stronger correlation with previous estimate
#est    lwr.ci    upr.ci
#1 0.8467414 0.7656606 0.9013361
CCC(traitsDF$x2016, traitsDF$preds2015)$rho.c
#est    lwr.ci    upr.ci
#1 0.7037067 0.5823608 0.7943768

#package
errors1516 <- data.frame(error = c(CCC(traitsDF$x2016, traitsDF$x2015)$blalt$delta,
                               CCC(traitsDF$x2016, traitsDF$preds2015)$blalt$delta),
                     type = c(rep("past estimate",71), rep("trait-based estimate",71)),
                     species = rep(traitsDF$Species, 2))

## 2016 to 2017 ------------------------------------------------------------

#we'll use the 2016 data to make trait predictions 
lm1 <- lm(x2016 ~ log(Mass) + ForStrat.ground + Habitat + 
            Trophic.Level, data=traitsDF)
traitsDF$preds2016 <- predict(lm1)

#we'll then compare it to what was actually observed in 2015
summary(abs(traitsDF$x2017 - traitsDF$x2016)) #compare to past estimate
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0306  1.1366  2.9019  6.0088  7.4847 41.4089

summary(abs(traitsDF$x2017 - traitsDF$preds2016)) #compare to trait-based estimate
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.6034  3.4449  7.0159  9.8701 14.9256 33.5758

CCC(traitsDF$x2017, traitsDF$x2016)$rho.c # corr greater with trait-based estimate
#est    lwr.ci    upr.ci
#1 0.8401798 0.7570643 0.8965335

CCC(traitsDF$x2017, traitsDF$preds2016)$rho.c
#est    lwr.ci    upr.ci
#1 0.6790767 0.5377321 0.7832938

#package
errors1617 <- data.frame(error = c(CCC(traitsDF$x2017, traitsDF$x2016)$blalt$delta,
                               CCC(traitsDF$x2017, traitsDF$preds2016)$blalt$delta),
                     type = c(rep("past estimate",71), rep("trait-based estimate",71)),
                     species = rep(traitsDF$Species, 2))

# final set

lm1 <- lm(x2017 ~ log(Mass) + ForStrat.ground + Habitat + 
            Trophic.Level, data=traitsDF)
traitsDF$preds2017 <- predict(lm1)

# Tests --------------------------------------------------------------

test1415 <- errors1415 %>% 
  dplyr::select(species, error, type) %>%
  pivot_wider(names_from = "type", values_from = "error") %>%
  janitor::clean_names()

mean(abs(test1415$past_estimate) > abs(test1415$trait_based_estimate)) #41%
species1415_diff <- test1415 %>% filter(abs(past_estimate) > abs(trait_based_estimate)) # species for which past estimate is larger
wilcox.test(abs(test1415$past_estimate), abs(test1415$trait_based_estimate), Paired=TRUE, exact=FALSE)
#W = 2059, p-value = 0.05999

test1516 <- errors1516 %>% 
  dplyr::select(species, error, type) %>%
  pivot_wider(names_from = "type", values_from = "error") %>%
  janitor::clean_names()

mean(abs(test1516$past_estimate) > abs(test1516$trait_based_estimate)) #34%
species1516_diff <- test1516 %>% filter(abs(past_estimate) > abs(trait_based_estimate)) # species past errorlarger
wilcox.test(abs(test1516$past_estimate), abs(test1516$trait_based_estimate), Paired=TRUE, exact=FALSE)
#W = 1771, p-value = 0.002243

test1617 <- errors1617 %>% 
  dplyr::select(species, error, type) %>%
  pivot_wider(names_from = "type", values_from = "error") %>%
  janitor::clean_names()

mean(abs(test1617$past_estimate) > abs(test1617$trait_based_estimate)) #31%
species1617_diff <- test1617 %>% filter(abs(past_estimate) > abs(trait_based_estimate))
wilcox.test(abs(test1617$past_estimate), abs(test1617$trait_based_estimate), Paired=TRUE, exact=FALSE)
#W = 1504, p-value = 3.393e-05

#spcecies for which past estimate error is larger for all years
intersect(species1415_diff$species[species1415_diff$species %in% species1516_diff$species],
          species1516_diff$species[species1516_diff$species %in% species1617_diff$species])

# trait-based prediction better for all years:
#[1] "Buteo buteo"        "Corvus corone"      "Hippolais icterina" "Oenanthe oenanthe" 
#[5] "Periparus ater"

# Combined plot ------------------------------------------------------

# combine data
errors1415 <- errors1415 %>%
  mutate(years="2014 to 2015")

errors1516 <- errors1516 %>%
  mutate(years="2015 to 2016")

errors1617 <- errors1617 %>%
  mutate(years="2016 to 2017")

errors_all <- rbind(errors1415, errors1516, errors1617)

# plot the data (eventually forming Fig 5b, combined with other plots below)
panelb <- ggplot(errors_all, aes(x = type, y = abs(error/100), group = species)) +
  geom_point(aes(color=years)) +
  geom_line(aes(color=years), alpha = 0.3) +  
  xlab("Approach") +
  ylab("|Detection error|") +
  facet_wrap(~ years) +
  scale_colour_manual(
    values = c(
      "2014 to 2015" = swatch()[1],
      "2015 to 2016" = swatch()[5],
      "2016 to 2017" = swatch()[9]
    )) +
  scale_x_discrete(labels = c("past estimate" = "Past estimate", "trait-based estimate" = "Trait-based estimate")) +
  theme(legend.position = "none",
        strip.text = element_text(size = 12))
panelb

# Population size estimate --------------------------------------------

# lets see the implications of using the predicted/previous det prob instead of observed det prob for a year

## 2014/2015 ---------------------------------------------------------------

# lets first continue with the scenario that was want to predict for 2015 but only have data for 2014

#correct all the site level predictions by detection and sum abundances across all sites (= total pop size)
popOutput1415 <- lapply(sort(unique(traitsDF$Species)), function(myspecies){
  
  #get raw site level predictions of raw abundance
  sitePreds_sp <- list.files("../data/models", full.names = TRUE) %>%
    str_subset("xggridpreds") %>%
    str_subset(strsplit(myspecies,"/")[[1]][1]) %>%
    readRDS() %>%
    dplyr::select(species:preds) %>%
    mutate(terr_cover = 1- buffer_seawater) %>%
    mutate(Rel_Abund = (preds * terr_cover))
  
  #correct the predictions by the two estimates of det prob (previous estimate or trait-based estimate)
  sitePreds_sp <- sitePreds_sp %>%
    inner_join(., traitsDF, by = join_by(species==Species)) %>%
    mutate(True_Abs_Abund = Rel_Abund * 100/x2015, # using truth
           Prev_Abs_Abund = Rel_Abund * 100/x2014, # using previous estimate
           Trait_Abs_Abund = Rel_Abund * 100/preds2014) # using trait-based estimate
  
  #calculate the total population size using the observed and predicted det prob
  sitePreds_sp %>%
    summarise(True_total_pop = sum(True_Abs_Abund),
              Prev_total_pop = sum(Prev_Abs_Abund),
              Trait_total_pop = sum(Trait_Abs_Abund)) %>%
    mutate(Pop_error_w_previous = (True_total_pop - Prev_total_pop)/True_total_pop,
           Pop_error_w_traits = (True_total_pop - Trait_total_pop)/True_total_pop) %>%
    add_column(Species = myspecies)
  
}) %>% bind_rows()

saveRDS(popOutput1415, file="outputs/popOutput1415.rds")

## 2015/2016 ---------------------------------------------------------------

#repeat 

popOutput1516 <- lapply(sort(unique(traitsDF$Species)), function(myspecies){
  
  sitePreds_sp <- list.files("../data/models", full.names = TRUE) %>%
    str_subset("xggridpreds") %>%
    str_subset(strsplit(myspecies,"/")[[1]][1]) %>%
    readRDS() %>%
    dplyr::select(species:preds) %>%
    mutate(terr_cover = 1- buffer_seawater) %>%
    mutate(Rel_Abund = (preds * terr_cover))
  
  #correct the predictions by the two estimates of det prob (previous estimate or trait-based estimate)
  sitePreds_sp <- sitePreds_sp %>%
    inner_join(., traitsDF, by = join_by(species==Species)) %>%
    mutate(True_Abs_Abund = Rel_Abund * 100/x2016, # using truth
           Prev_Abs_Abund = Rel_Abund * 100/x2015, # using previous estimate
           Trait_Abs_Abund = Rel_Abund * 100/preds2015) # using trait-based estimate
  
  #calculate the total population size using the observed and predicted det prob
  sitePreds_sp %>%
    summarise(True_total_pop = sum(True_Abs_Abund),
              Prev_total_pop = sum(Prev_Abs_Abund),
              Trait_total_pop = sum(Trait_Abs_Abund)) %>%
    mutate(Pop_error_w_previous = (True_total_pop - Prev_total_pop)/True_total_pop,
           Pop_error_w_traits = (True_total_pop - Trait_total_pop)/True_total_pop) %>%
    add_column(Species = myspecies)
  
}) %>% bind_rows()

saveRDS(popOutput1516, file="outputs/popOutput1516.rds")

## 2016/2017 ---------------------------------------------------------------

#repeat

popOutput1617 <- lapply(sort(unique(traitsDF$Species)), function(myspecies){
  
  sitePreds_sp <- list.files("../data/models", full.names = TRUE) %>%
    str_subset("xggridpreds") %>%
    str_subset(strsplit(myspecies,"/")[[1]][1]) %>%
    readRDS() %>%
    dplyr::select(species:preds) %>%
    mutate(terr_cover = 1- buffer_seawater) %>%
    mutate(Rel_Abund = (preds * terr_cover))
  
  #correct the predictions by the two estimates of det prob (previous estimate or trait-based estimate)
  sitePreds_sp <- sitePreds_sp %>%
    inner_join(., traitsDF, by = join_by(species==Species)) %>%
    mutate(True_Abs_Abund = Rel_Abund * 100/x2017, # using truth
           Prev_Abs_Abund = Rel_Abund * 100/x2016, # using previous estimate
           Trait_Abs_Abund = Rel_Abund * 100/preds2016) # using trait-based estimate
  
  #calculate the total population size using the observed and predicted det prob
  sitePreds_sp %>%
    summarise(True_total_pop = sum(True_Abs_Abund),
              Prev_total_pop = sum(Prev_Abs_Abund),
              Trait_total_pop = sum(Trait_Abs_Abund)) %>%
    mutate(Pop_error_w_previous = (True_total_pop - Prev_total_pop)/True_total_pop,
           Pop_error_w_traits = (True_total_pop - Trait_total_pop)/True_total_pop) %>%
    add_column(Species = myspecies)
  
}) %>% bind_rows()

saveRDS(popOutput1617, file="outputs/popOutput1617.rds")

## Combine plots ----------------------------------------------------------------------

popOutput1415 <- readRDS("outputs/popOutput1415.rds")
popOutput1516 <- readRDS("outputs/popOutput1516.rds")
popOutput1617 <- readRDS("outputs/popOutput1617.rds")

# clean and combine the data 
popOutput1415 <- popOutput1415 %>%
  mutate(year="2014 to 2015") %>%
  dplyr::select(Species, Pop_error_w_previous, Pop_error_w_traits, year) %>%
  pivot_longer(cols=2:3, names_to ="Type", values_to ="error")

popOutput1516 <- popOutput1516 %>%
  mutate(year="2015 to 2016") %>%
  dplyr::select(Species, Pop_error_w_previous, Pop_error_w_traits, year) %>%
  pivot_longer(cols=2:3, names_to ="Type", values_to ="error")

popOutput1617 <- popOutput1617 %>%
  mutate(year="2016 to 2017") %>%
  dplyr::select(Species, Pop_error_w_previous, Pop_error_w_traits, year) %>%
  pivot_longer(cols=2:3, names_to ="Type", values_to ="error")

popOutput_all <- rbind(popOutput1415, popOutput1516, popOutput1617)

#fix labels
popOutput_all$Type <- ifelse(popOutput_all$Type == "Pop_error_w_previous", "Past estimate", "Trait-based estimate")

# Create the faceted plot
panelc <- ggplot(popOutput_all, aes(x = Type, y = abs(error*100))) +
  geom_boxplot(aes(fill = year), alpha = 0.5) +  # Boxplots colored by year
  geom_jitter(aes(colour = year), alpha = 0.3, width = 0.2) +  # Points colored by year
  xlab("Prediction approach") + 
  ylab("|Population error (%)|") +
  scale_fill_manual(
    values = c(
      "2014 to 2015" = swatch()[1],
      "2015 to 2016" = swatch()[5],
      "2016 to 2017" = swatch()[9]
    )) +
  scale_colour_manual(
    values = c(
      "2014 to 2015" = swatch()[1],
      "2015 to 2016" = swatch()[5],
      "2016 to 2017" = swatch()[9]
    )) +
  labs(fill = "Year Comparison", colour = "Year Comparison") +  # Legend titles
  facet_wrap(~ year, scales = "free_y") +
  theme(legend.position = "none",
        strip.text = element_text(size = 12))  # Remove the legend
panelc

# Combine all three plots ------------------------------------------------------

combined_fig <- (patchwork::wrap_elements(panela) / 
                   patchwork::wrap_elements(panelb) /
                   patchwork::wrap_elements(panelc)) +
  plot_annotation(tag_levels = list(c("a", "b", "c")))
combined_fig

ggsave("plots/Figure_5.png", height=9, width=9, units="in")

# Tests --------------------------------------------------------------

test1415 <- popOutput1415 %>% 
  dplyr::select(Species, error, Type) %>%
  pivot_wider(names_from = "Type", values_from = "error") %>%
  janitor::clean_names()

summary(abs(test1415$pop_error_w_previous))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.03569 0.10466 0.14467 0.19837 0.64634 
summary(abs(test1415$pop_error_w_traits))
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.003579 0.064958 0.123070 0.179672 0.262046 0.684516
mean(abs(test1415$pop_error_w_previous) > abs(test1415$pop_error_w_traits)) #44%
species1415_diff <- test1415 %>% filter(abs(pop_error_w_previous) > abs(pop_error_w_traits))
wilcox.test(abs(test1415$pop_error_w_previous), abs(test1415$pop_error_w_traits), Paired=TRUE, exact=FALSE)
#W = 2119, p-value = 0.1018

test1516 <- popOutput1516 %>% 
  dplyr::select(Species, error, Type) %>%
  pivot_wider(names_from = "Type", values_from = "error") %>%
  janitor::clean_names()

summary(abs(test1516$pop_error_w_previous))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.03525 0.08898 0.12323 0.17486 0.59499
summary(abs(test1516$pop_error_w_traits))
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.006175 0.075810 0.142288 0.181300 0.283129 0.593797 
mean(abs(test1516$pop_error_w_previous) > abs(test1516$pop_error_w_traits)) #35%
species1516_diff <- test1516 %>% filter(abs(pop_error_w_previous) > abs(pop_error_w_traits))
wilcox.test(abs(test1516$pop_error_w_previous), abs(test1516$pop_error_w_traits), Paired=TRUE, exact=FALSE)
#W = 1799, p-value = 0.003264

test1617 <- popOutput1617 %>% 
  dplyr::select(Species, error, Type) %>%
  pivot_wider(names_from = "Type", values_from = "error") %>%
  janitor::clean_names()

summary(abs(test1617$pop_error_w_previous))
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.0007729 0.0192812 0.0594643 0.1235438 0.1192659 1.3827086
summary(abs(test1617$pop_error_w_traits))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.01393 0.07163 0.14368 0.18733 0.28974 0.62034
mean(abs(test1617$pop_error_w_previous) > abs(test1617$pop_error_w_traits)) #31%
species1617_diff <- test1617 %>% filter(abs(pop_error_w_previous) > abs(pop_error_w_traits))
wilcox.test(abs(test1617$pop_error_w_previous), abs(test1617$pop_error_w_traits), Paired=TRUE, exact=FALSE)
#W = 1540, p-value = 6.376e-05

intersect(species1415_diff$species[species1415_diff$species %in% species1516_diff$species],
          species1516_diff$species[species1516_diff$species %in% species1617_diff$species])
# "Buteo buteo"        "Corvus corone"      "Hippolais icterina" "Periparus ater"

#overall difference
medianPerSpecies <- bind_rows(popOutput1415, popOutput1516, popOutput1617) %>%
  group_by(Species, Type) %>%
  summarise(med = median(abs(error))) %>%
  ungroup()

summary(medianPerSpecies$med[medianPerSpecies$Type=="Pop_error_w_previous"])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.04082 0.08294 0.12359 0.17531 0.59499
summary(medianPerSpecies$med[medianPerSpecies$Type=="Pop_error_w_traits"])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.03264 0.08045 0.12720 0.17328 0.24504 0.59380

# end --------------------------------------------------------------------------
