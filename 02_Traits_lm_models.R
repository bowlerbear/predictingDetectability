# README:
# This script includes:
# (1) multiple regression to explore the relationship between species traits and mean detection probability
# (2) uses a leave-one-out approach with the model to test predictive ability
# (3) estimates the error from using these predicted estimate over the original estimate
# (4) same as (1) except for the effects of covariates on detection probability

# Libraries --------------------------------------------------------------------

library(tidyverse)
library(broom)
library(lubridate)
library(patchwork)
library(ggExtra)
library(ape)
library(phytools)
library(modelsummary)
library(ggthemr)
ggthemr("fresh") 

# Get estimated detection probabilities ---------------------------------------

#these are calculated in the script 01
detDF <- readRDS("outputs/distance_passerines_constant.rds")

# Get tree --------------------------------------------------------------------

#get phylo tree to test phylogenetic signal later

tree <- read.nexus("traits/output.nex")

#get tip labels
labels <- read.csv("traits/BLIOCPhyloMasterTax.csv", sep=";") %>% 
  filter(My.name!="") %>%
  dplyr::select(TipLabel, My.name) 

#add missing row for the hooded crow
hooded.crow <- data.frame(TipLabel = "Corvus_corone", My.name = "Corvus cornix")
labels <- bind_rows(labels, hooded.crow)

#order data frame in same as tree
labels <- labels[order(match(labels$TipLabel, tree[[1]]$tip.label)),]

# Get traits file --------------------------------------------------------------

#already compiled ecological and morphological traits for each species
traits <- readRDS("traits/traits.rds") # see Table S2 for information
head(traits)

#some small tweaks
traits$Habitat[traits$Habitat=="Woodland"] <- "Forest"
#for black redstart, change since only one value in rock, and species also highly urban
traits$Habitat[traits$Habitat=="Rock"] <- "Human Modified"

#combine
traitsDF <- detDF %>% inner_join(., traits, by="Species") 

# Explore det probs -------------------------------------------------------------

summary(traitsDF$P)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.2498  0.4353  0.5436  0.5512  0.6070  1.0000

#Species with minimum
traitsDF[which(traitsDF$P == min(traitsDF$P)),]
#Species with maximum
traitsDF[which(traitsDF$P == max(traitsDF$P)),]

# Plot relationships -----------------------------------------------------------

#We want to look at relationships between species traits and their detection probabilites
#Note: ESW is effective strip width

#Figure 2 of the paper
g1 <- qplot(Mass, ESW/100, data=traitsDF) + 
  stat_smooth(method="lm") +
  ylab("Detection probability") +
  xlab("Mass (g)") +
  scale_x_log10()

g2 <- qplot(ForStrat.ground, ESW/100, data=traitsDF) + 
  stat_smooth(method="lm")+
  ylab("Detection probability") +
  xlab("Ground foraging  %")

g3 <- qplot(Habitat, ESW/100, data=traitsDF, geom="boxplot") + 
  geom_jitter(width = 0.2,colour="black", alpha = 0.5)+
  ylab("Detection probability") +
  coord_flip() 

g4 <- qplot(Trophic.Level, ESW/100, data=traitsDF, geom="boxplot") +
  geom_jitter(width = 0.2,colour="black", alpha = 0.5)+
  ylab("Detection probability") +
  xlab("Trophic niche") +
  coord_flip()

g5 <- qplot(as.factor(flockSize), ESW/100, data=traitsDF, geom="boxplot") +
  geom_jitter(width = 0.2,colour="black", alpha = 0.5)+
  ylab("Detection probability") +
  xlab("Flocking propensity") +
  coord_flip()

# Arrange g4, g5, g2, and g3 
gRest <- (g2 + g3) /(g5 + g4) 

# Combine g1 with gRest, with g1 on top and gRest below
final_plot <- g1 / gRest +   
  plot_annotation(tag_levels = list(c("a", "b", "c", "d", "e"))) +
  plot_layout(guides = 'collect') 

# Print the final plot
final_plot

ggsave("plots/Figure_2.png", width=7, height=9, units="in")

# Multiple regression ------------------------------------------------------------------

lm1 <- lm(ESW ~ log(Mass) + ForStrat.ground + Habitat + 
             Trophic.Level + flockSize, data=traitsDF)
summary(lm1)

#pull out coefficients and plot them
coefDF <- data.frame(Param = names(coefficients(lm1)),
                     estimate = as.numeric(coefficients(lm1)),
                     lower = as.numeric(confint(lm1)[,1]),
                     upper = as.numeric(confint(lm1)[,2]))

# drop to significant variables
lm1 <- lm(ESW ~ log(Mass) + ForStrat.ground + Habitat + Trophic.Level, data=traitsDF)
summary(lm1)

#organise this table (included in Table S4)
modelsummary(lm1, statistic = c("conf.low", "conf.high", "p.value"), 
             shape = term ~ model + statistic,
             output="plots/Linear_Model_Ouput_all.docx")

# Test phylo signal --------------------------------------------------------------

#test directly on detection probability
labels$trait <-detDF$ESW[match(labels$My.name, detDF$Species)]
my.trait <-setNames(labels$trait,labels$TipLabel)

#loop through each tree
pVals <- sapply(1:length(tree), function(i){
  phylosig(tree[[i]], my.trait, method="lambda", test=TRUE, nsim=999)$P
})
summary(pVals)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#3.106e-06 6.129e-06 7.734e-06 9.220e-06 1.172e-05 2.437e-05 

#test on model residuals after accounting for trait effects
detDF$resids <- resid(lm1)
labels$trait <-detDF$resids[match(labels$My.name, detDF$Species)]
my.trait <-setNames(labels$trait,labels$TipLabel)

#loop through each tree
pVals <- sapply(1:length(tree), function(i){
   phylosig(tree[[i]], my.trait, method="lambda", test=TRUE, nsim=999)$P
})
summary(pVals)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1       1       1       1       1       1

## Leave-one-out predictions -------------------------------------------------------------------

leave_one_out_function <- function(species_name){

  dat_filtered <- traitsDF %>%
    dplyr::filter(Species != species_name)

  mod <- lm(ESW ~ log(Mass) + ForStrat.ground + Habitat +
              Trophic.Level,
            data=dat_filtered)

  new_dat <- traitsDF %>%
    dplyr::filter(Species==species_name)

  predicted_ESW <- predict(mod, new_dat)

  loo_summary <- data.frame(Species=species_name,
                            observed_ESW=new_dat$ESW,
                            predicted_ESW=predicted_ESW)

  return(loo_summary)

}

leave_one_out_results <- bind_rows(lapply(unique(traitsDF$Species),
                                          leave_one_out_function))

saveRDS(leave_one_out_results, file="outputs/leave_one_out_results.rds")

## Compare trait predictions and original estimates -----------------------------------------------

leave_one_out_results <- readRDS("outputs/leave_one_out_results.rds")

#concordance correlation coefficient
DescTools::CCC(leave_one_out_results$observed_ESW,
               leave_one_out_results$predicted_ESW)$rho.c
#est    lwr.ci    upr.ci
#1 0.6656559 0.5239742 0.7714912

DescTools::CCC(leave_one_out_results$observed_ESW,
               leave_one_out_results$predicted_ESW)$C.b
#0.969114

#mean dfference between observed and predicted effective strip width in m
summary(abs(leave_one_out_results$observed_ESW-leave_one_out_results$predicted_ESW))
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.03852  3.81933  6.72528  9.55518 13.74976 39.58919

cor.test(leave_one_out_results$observed_ESW,
         leave_one_out_results$predicted_ESW)

# Pearson's product-moment correlation
# 
# data:  leave_one_out_results$observed_ESW and leave_one_out_results$predicted_ESW
# t = 7.9072, df = 70, p-value = 2.708e-11
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.5413442 0.7924410
# sample estimates:
#       cor 
# 0.6868706


#Figure 3
Fig3a <- ggplot(leave_one_out_results, 
                aes(x=observed_ESW/100, y=predicted_ESW/100))+
  geom_point()+
  theme_minimal()+
  theme(axis.text=element_text(color="black"),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5))+
  xlab("Detection probability")+
  ylab("Trait-based detection probability")+
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  xlim(0.15,1) + ylim(0.15,1)+
  annotate("text", x=0.4, y=0.95, label= "rho = 0.69, MAD = 9.56 m") 
  #geom_smooth(method="lm")
Fig3a

saveRDS(Fig3a, file="plots/Figure_3a.rds")

# Body size regression ------------------------------------------------------------------------

lm_size <- lm(ESW ~ log(Mass), data=traitsDF)
summary(lm_size)

## Leave-one-out predictions -------------------------------------------------------------------

# model only use size

leave_one_out_function <- function(species_name){

  dat_filtered <- traitsDF %>%
    dplyr::filter(Species != species_name)

  mod <- lm(ESW ~ log(Mass), data=dat_filtered)

  new_dat <- traitsDF %>%
    dplyr::filter(Species==species_name)

  predicted_ESW <- predict(mod, new_dat)

  loo_summary <- data.frame(Species=species_name,
                            observed_ESW=new_dat$ESW,
                            predicted_ESW=predicted_ESW)

  return(loo_summary)

}

leave_one_out_results_size <- bind_rows(lapply(unique(traitsDF$Species),
                                          leave_one_out_function))
saveRDS(leave_one_out_results_size, file="outputs/leave_one_out_results_size.rds")

## Compare trait predictions and original estimates -----------------------------------------------

leave_one_out_results_size <- readRDS("outputs/leave_one_out_results_size.rds")

DescTools::CCC(leave_one_out_results_size$observed_ESW,
               leave_one_out_results_size$predicted_ESW)$rho.c
#est    lwr.ci    upr.ci
#1 0.5506402 0.3951272 0.6754714

DescTools::CCC(leave_one_out_results_size$observed_ESW,
               leave_one_out_results_size$predicted_ESW)$C.b
#0.9046672

summary(abs(leave_one_out_results_size$observed_ESW-leave_one_out_results_size$predicted_ESW))
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.05897  4.73465  9.43455 10.97224 16.36088 40.90627

cor.test(leave_one_out_results_size$observed_ESW,
         leave_one_out_results_size$predicted_ESW)

#Pearson's product-moment correlation

# data:  leave_one_out_results_size$observed_ESW and leave_one_out_results_size$predicted_ESW
# t = 6.4183, df = 70, p-value = 1.41e-08
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.4388841 0.7364839
# sample estimates:
#       cor 
# 0.6086661 

#Figure 3
Fig3b <- ggplot(leave_one_out_results_size, aes(x=observed_ESW/100, y=predicted_ESW/100))+
  geom_point()+
  theme_minimal()+
  theme(axis.text=element_text(color="black"),
        axis.line.x = element_line(color = "black", size = 0.5),
        axis.line.y = element_line(color = "black", size = 0.5))+
  xlab("Detection probability")+
  ylab("Trait-based detection probability")+
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  xlim(0.15,1) + ylim(0.15,1)+
  annotate("text", x=0.4, y=0.95, label= "rho = 0.61, MAD = 10.97 m")
  #geom_smooth(method="lm")
Fig3b

saveRDS(Fig3b, file="plots/Figure_3b.rds")

# Population size estimate --------------------------------------------

#We now want to understand the implications of using the predicted detection probability 
#instead of direct estimates based on observed data 

constantESW <- readRDS("outputs/distance_passerines_constant.rds") # see script 01
bootsESW <- readRDS("outputs/distance_passerines_constant_boot.rds") #see HPC_Detection_bootstraps (used to calc CIs)
leave_one_out_results <- readRDS("outputs/leave_one_out_results.rds") # made above

#correct all the site level predictions by detection and sum abundances across all sites (= total pop size)
popOutput <- data.frame()
allspecies <- sort(unique(constantESW$Species))

for(s in 1:length(allspecies)){

  myspecies <- allspecies[s]

  #get site-level relative abundance predictions (not corrected for detection)
  sitePreds <- list.files("../data/models", full.names = TRUE) %>% # site-level preds - seee HPC_Abundance_models
                str_subset("xggridpreds") %>%
                str_subset(strsplit(myspecies,"/")[[1]][1]) %>%
                readRDS()
  
  #get detection probabilities
  constantSp <- constantESW %>% filter(Species == myspecies)
  bootSp <- bootsESW %>% filter(Species == myspecies)

  #first use preds column x median detection for species
  meanPreds <- sitePreds %>%
                dplyr::select(species:preds) %>%
                mutate(terr_cover = 1- buffer_seawater) %>%
                add_column(ESW = constantSp$ESW) %>%
                mutate(Pop_size = (preds * terr_cover)/(ESW/100)) #direct estimate of pop size

  #add on trait estimate
  meanPreds$predESW <- leave_one_out_results$predicted_ESW[match(meanPreds$species,
                                                                 leave_one_out_results$Species)]
  #and calc trait-based pop size
  meanPreds <- meanPreds %>%
                  mutate(Trait_pop_size = (preds * terr_cover)/(predESW/100))

  #also get bootstrapped samples to calculate the CIs
  bootPreds <- sitePreds %>% dplyr::select(-species, -x, -y, -buffer_seawater, -preds)

  bootPop <- as.numeric()

  for(i in 1:ncol(bootPreds)){

    tempPreds <- sitePreds %>%
      dplyr::select(species:buffer_seawater) %>%
      add_column(preds = bootPreds[,i]) %>%
      mutate(terr_cover = 1- buffer_seawater) %>%
      add_column(ESW = bootSp$vals[i]) %>%
      mutate(Pop_size = (preds * terr_cover)/(ESW/100))

    bootPop[i] <- sum(tempPreds$Pop_size)

  }

  #package up
  popOutput_new <- data.frame(Species = myspecies,
                              Obs_Ab_Abund = sum(meanPreds$Pop_size), #direct estimate of pop size
                              Preds_Ab_Abund = sum(meanPreds$Trait_pop_size), #predicted estimate using traits
                              Obs_Lower_Abund = as.numeric(quantile(bootPop, 0.025)), # lower 95% of direct
                              Obs_Upper_Abund = as.numeric(quantile(bootPop, 0.975))) # upper 95% of direct

  popOutput <- bind_rows(popOutput, popOutput_new)

}
saveRDS(popOutput, file="outputs/popOutput.rds")

popOutput <- readRDS("outputs/popOutput.rds") %>%
            inner_join(.,leave_one_out_results, by="Species") %>%
            mutate(Pop_error = (Obs_Ab_Abund - Preds_Ab_Abund)/Obs_Ab_Abund,
                   Pop_u_error = (Obs_Ab_Abund - Obs_Lower_Abund)/Obs_Ab_Abund, 
                   Pop_l_error = (Obs_Ab_Abund - Obs_Upper_Abund)/Obs_Ab_Abund, 
                   ESW_error = observed_ESW - predicted_ESW,
                   P_error = ESW_error/100) 

#For Table S6:
write.csv(popOutput, file = "outputs/popestimates_lm.csv", row.names=FALSE)

#for how many species is the trait-based estimates within the bootstraped CI?
mean(popOutput$Obs_Upper_Abund > popOutput$Preds_Ab_Abund &
       popOutput$Obs_Lower_Abund < popOutput$Preds_Ab_Abund)
#82%

#which species?
(outSpecies <- popOutput$Species[!(popOutput$Obs_Upper_Abund > popOutput$Preds_Ab_Abund &
                    popOutput$Obs_Lower_Abund < popOutput$Preds_Ab_Abund)])

#plot the error
g_error <- ggplot(popOutput) +
  geom_point(aes(x = ESW_error/100, y = Pop_error*100)) +
  geom_hline(yintercept =  median(popOutput$Pop_u_error)*100, linetype ="dashed", color = "red") + # null error
  geom_hline(yintercept =  median(popOutput$Pop_l_error)*100, linetype ="dashed", color = "red") + # null error
  geom_hline(yintercept = 0, linetype ="dashed") +
  geom_vline(xintercept = 0, linetype ="dashed") +
  xlab("Detection error") +
  ylab("Population error (%)")

g_Error <- ggMarginal(g_error, type="density", fill = "lightblue")

#plot the absolute error
g_abs <- ggplot(popOutput) +
  geom_point(aes(x = abs(ESW_error/100), y = abs(Pop_error*100))) +
  geom_hline(yintercept =  median(abs(c(popOutput$Pop_u_error, popOutput$Pop_l_error))*100), linetype ="dashed", color = "red") + # null error
  xlab("|Detection error|") +
  ylab("|Population error (%)|")

g_Abs <- ggMarginal(g_abs, type="boxplot", fill = "lightblue")

lm_error_abs <- lm(abs(Pop_error) ~ abs(P_error), data=popOutput)
summary(lm_error_abs)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         0.01559    0.01094   1.425    0.159    
# abs(ESW_error/100)  1.71389    0.08752  19.582   <2e-16 ***

summary(abs(popOutput$Pop_error))
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.001078 0.070503 0.125278 0.179357 0.255967 0.737884 

#also look at the correlation between the population size estimates
(gScatter <- ggplot(popOutput, 
       aes(x=Obs_Ab_Abund, y=Preds_Ab_Abund))+
  geom_point()+
  theme_minimal()+
  theme(axis.text=element_text(color="black"),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5))+
  xlab("Estimated population size")+
  ylab("Trait-based population size")+
  #scale_x_log10()+scale_y_log10()+
  geom_abline(intercept = 0, slope = 1, linetype="dashed"))

DescTools::CCC(popOutput$Obs_Ab_Abund,popOutput$Preds_Ab_Abund)$rho.c # 0.96
DescTools::CCC(popOutput$Obs_Ab_Abund,popOutput$Preds_Ab_Abund)$C.b # 0.99
cor(popOutput$Obs_Ab_Abund,popOutput$Preds_Ab_Abund) # 0.97

#combine into a figure 4
combined_fig4 <- patchwork::wrap_elements(gScatter)/((patchwork::wrap_elements(g_Error) +
                                 patchwork::wrap_elements(g_Abs))) + 
  plot_annotation(tag_levels = list(c("a", "b", "c"))) 

ggsave("plots/Figure_4.png",width = 6, height = 6)

# Covariate effects ------------------------------------------------------------

#We can also ask weather we can use traits to predict environmental effects on detection probability
#this output includes the effect of two covariates tested
#forest cover - we expect higher forest cover leads to lower detections
#paths/roads - we expect some species to be attracted to roads and others to be deterred

detDF <- readRDS("outputs/distance_passerines_covs.rds") # see script 01

traitsDF <- detDF %>%
  inner_join(., traits, by="Species") 

# how many significant effects?
mean(abs(traitsDF$Z_pathroad)>1.96)
mean(abs(traitsDF$Z_forest)>1.96)
nrow(subset(traitsDF, abs(Z_pathroad)>1.96))
nrow(subset(traitsDF, abs(Z_forest)>1.96))

subset(traitsDF, Z_pathroad>1.96)$Species
subset(traitsDF, Z_pathroad<(-1.96))$Species 
#"Motacilla alba","Passer domesticus" "Pica pica","Pyrrhula pyrrhula",
#"Sturnus vulgaris"
subset(traitsDF, Z_forest<(-1.96))$Species #"Sturnus vulgaris"
subset(traitsDF, Z_forest>(1.96))$Species # "Erithacus rubecula"

#regression output 
lm1 <- lm(Z_forest ~ Habitat, data=traitsDF)
summary(lm1)
lm1 <- lm(Z_pathroad ~ Habitat, data=traitsDF)
summary(lm1)

g1 <- qplot(Habitat, Z_forest, data=traitsDF, geom="boxplot") + 
  geom_jitter(width = 0.2,colour="black", alpha = 0.5)+
  coord_flip() +
  ylab("Effect of forest cover")

g2 <- qplot(Habitat, Z_pathroad, data=traitsDF, geom="boxplot") + 
  geom_jitter(width = 0.2,colour="black", alpha = 0.5)+
  coord_flip() +
  ylab("Effect of path/road")

g1 + g2


#plot line effects
distance_passerines_covpreds <- readRDS("outputs/distance_passerines_covpreds.rds") #see script 01

temp <- distance_passerines_covpreds  %>%
  filter(Type=="path effects")

g3 <- ggplot(temp) +
  geom_line(data = subset(temp,abs(paths_zval)>2), 
                aes(x = lines_pathroad, y = ESWpreds/100, group = Species))+
    geom_line(data = subset(temp,abs(paths_zval)<2), 
                  aes(x = lines_pathroad, y = ESWpreds/100, group = Species), alpha = 0.2) +
  xlab("log Path/Road cover") + ylab("Detection probability")

#plot forest effects
temp <- distance_passerines_covpreds  %>%
  filter(Type=="forest effects")

g4 <- ggplot(temp) +
  geom_line(data = subset(temp,abs(forest_zval)>2), 
            aes(x = lines_forest, y = ESWpreds/100, group = Species))+
  geom_line(data = subset(temp,abs(forest_zval)<2), 
            aes(x = lines_forest, y = ESWpreds/100, group = Species), alpha = 0.2)+
  xlab("log Forest cover") + ylab("Detection probability")


(g1+g2) / (g3+g4) +   
  plot_annotation(tag_levels = list(c("a", "b", "c", "d"))) +
  plot_layout(guides = 'collect') 

ggsave("plots/Figure_S3.png", width = 6, height = 6)

# end --------------------------------------------------------------------------