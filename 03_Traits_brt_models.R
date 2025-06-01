# README:
# This script includes:
# (1) boosted regression trees to explore the relationship between species traits and mean detection probability
# (2) uses a leave-one-out approach with the model to test predictive ability
# (3) estimates the error from using these predicted estimate over the original estimate

# Libraries --------------------------------------------------------------------

library(tidyverse)
library(broom)
library(lubridate)
library(patchwork)
library(ggExtra)
library(dismo)
library(gbm)
library(rsample)
library(xgboost)
library(fastDummies)
library(ggthemr)
ggthemr("fresh") 

set.seed(1039)


# Get estimated detection probabilities ---------------------------------------

#these are calculated in the script 01
detDF <- readRDS("outputs/distance_passerines_constant.rds")

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

# BRT ------------------------------------------------------------------------

## Hypergrid ---------------------------------------------------------------

#explore the best parameters

#make test and training dataset
cust_split <- traitsDF %>%
        initial_split(prop = 0.7)

train <- training(cust_split)
test <- testing(cust_split)

# Set up the grid for hyperparameters
n.trees <- c(2000, 4000, 6000)
learning.rate <- c(0.0001, 0.0005, 0.001)
tree.complexity <- c(1, 3, 5)
step <- c(10, 50, 100)


# Initialize a dataframe to store results
results <- data.frame(n.trees = integer(),
                      learning.rate = numeric(),
                      tree.complexity = numeric(),
                      step.size = numeric(),
                      RMSE = numeric())

# Cross-validate to find the best parameters
for (n in n.trees) {
    for (lr in learning.rate) {
      for (tc in tree.complexity) {
        for (s in step) {
        
        # Run gbm.step with current parameters
        model <- gbm.step(data = train, 
                         gbm.x = 9:ncol(train), 
                         gbm.y = "ESW", 
                         family = "gaussian",
                         bag.fraction = 0.7,
                         step.size = s,
                         n.trees = n,
                         tree.complexity = tc, 
                         learning.rate = lr)
        
        # Get predictions on test dataset
        predictions <- predict(model, newdata = test, n.trees = model$n.trees)
        
        # Calculate RMSE
        rmse_value <- sqrt(mean((test$ESW - predictions)^2))
        
        # Store the results
        results <- rbind(results, data.frame(n.trees = n,
                                             learning.rate = lr,
                                             step.size = s,
                                             tree.complexity = tc,
                                             RMSE = rmse_value))
      }
    }
    }
}

# Find the best parameters based on RMSE
best_params <- results[which.min(results$RMSE),]
print(best_params)

results %>%
  arrange(RMSE)

## Fit full model ------------------------------------------------------------------------

gbm1 <- gbm.step(data=traitsDF, gbm.x =9:ncol(traitsDF), gbm.y = "ESW", 
                 family = "gaussian",
                 tree.complexity = 1, 
                 learning.rate = 0.0005, 
                 step.size = 100,
                 bag.fraction = 0.7, 
                 n.trees = 4000)

gbm.perf(gbm1)

# we can look at their relative like this:
summary(gbm1)
# var      rel.inf
# Tarsus.Length                 Tarsus.Length 37.160442535
# TailU_MEAN                       TailU_MEAN 13.616194401
# ForStrat.ground             ForStrat.ground 13.018516518
# Wing.Length                     Wing.Length  9.302969143
# Beak.Length_Culmen       Beak.Length_Culmen  8.753716467
# Egg_MASS                           Egg_MASS  4.559070239
# Nest.type                         Nest.type  3.065443495
# Migration                         Migration  2.885368594
# Clutch_MAX                       Clutch_MAX  1.821525639
# Clutch_MEAN                     Clutch_MEAN  1.150871095

png("plots/Figure_S2.png", width=2000, height=2000, res=300)
gbm.plot(gbm1, n.plots=10)
dev.off()

## Leave-one-out predictions -------------------------------------------------------------------

leave_one_out_function <- function(species_name){

  dat_filtered <- traitsDF %>%
    dplyr::filter(Species != species_name)

  gbm1 <- gbm.step(data=dat_filtered, gbm.x =9:ncol(dat_filtered), gbm.y = "ESW", family = "gaussian",
                   tree.complexity = 1, learning.rate = 0.0005, step.size = 100,
                   bag.fraction = 0.7, n.folds = 5, n.trees = 4000)

  new_dat <- traitsDF %>%
    dplyr::filter(Species==species_name)

  predicted_ESW <- predict.gbm(gbm1, new_dat,
                               n.trees=gbm1$gbm.call$best.trees,
                               type="response")

  loo_summary <- data.frame(Species=species_name,
                            observed_ESW=new_dat$ESW,
                            predicted_ESW=predicted_ESW)

  return(loo_summary)

}

# the following code will take over an hour to run
leave_one_out_results <- bind_rows(lapply(unique(traitsDF$Species),
                                         leave_one_out_function))
saveRDS(leave_one_out_results, file="outputs/leave_one_out_results_brt.rds")

## Compare trait predictions and original estimates -----------------------------------------------

leave_one_out_results <- readRDS("outputs/leave_one_out_results_brt.rds")

DescTools::CCC(leave_one_out_results$observed_ESW,
               leave_one_out_results$predicted_ESW)$rho.c
#est    lwr.ci   upr.ci
#1 0.5684925 0.4211185 0.686621

DescTools::CCC(leave_one_out_results$observed_ESW,
               leave_one_out_results$predicted_ESW)$C.b
#0.8933001

summary(abs(leave_one_out_results$observed_ESW-leave_one_out_results$predicted_ESW))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.2758  4.2698  9.2363 10.4287 14.2625 38.1149

cor.test(leave_one_out_results$observed_ESW,
         leave_one_out_results$predicted_ESW)

# Pearson's product-moment correlation
# 
# data:  leave_one_out_results$observed_ESW and leave_one_out_results$predicted_ESW
# t = 6.9027, df = 70, p-value = 1.878e-09
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.4747161 0.7565274
# sample estimates:
#       cor 
# 0.6363959

#for Figure 3
Fig3c <- ggplot(leave_one_out_results, aes(x=observed_ESW/100, y=predicted_ESW/100))+
  geom_point()+
  theme_minimal()+
  theme(axis.text=element_text(color="black"),
        axis.line.x = element_line(color = "black", size = 0.5),
        axis.line.y = element_line(color = "black", size = 0.5))+
  xlab("Detection probability")+
  ylab("Trait-based detection probability")+
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  xlim(0.15,1) + ylim(0.15,1)+
  annotate("text", x=0.40, y=0.95, label= "rho = 0.64, MAD = 10.43 m")
  #geom_smooth(method="lm")
Fig3c

## Fit simplified model --------------------------------------------------------

gbm1_simp <- gbm.simplify(gbm1, n.drops = "auto")
gbm1_simp$deviance.summary

gbm1_reduced <- gbm.step(data=traitsDF, gbm.x =gbm1_simp$pred.list[[54]], gbm.y = "ESW", family = "gaussian",
                         tree.complexity = 1, learning.rate = 0.0005, step.size = 100,
                         bag.fraction = 0.7, n.folds = 5, n.trees = 4000)

summary(gbm1_reduced)
# var   rel.inf
# Tarsus.Length           Tarsus.Length 34.115895
# TailU_MEAN                 TailU_MEAN 13.893536
# ForStrat.ground       ForStrat.ground 12.957055
# Wing.Length               Wing.Length  9.619072
# Beak.Length_Culmen Beak.Length_Culmen  9.419607
# Nest.type                   Nest.type  5.323355
# Migration                   Migration  5.235285
# Egg_MASS                     Egg_MASS  5.015152
# Clutch_MEAN               Clutch_MEAN  3.070287
# LengthU_MEAN             LengthU_MEAN  1.350757

## Leave-one-out predictions -------------------------------------------------------------------

leave_one_out_function <- function(species_name){

  dat_filtered <- traitsDF %>%
    dplyr::filter(Species != species_name)

  gbm1 <- gbm.step(data=dat_filtered, gbm.x =gbm1_simp$pred.list[[54]], gbm.y = "ESW", family = "gaussian",
                   tree.complexity = 1, learning.rate = 0.0005, step.size = 100,
                   bag.fraction = 0.7, n.folds = 5, n.trees = 4000)

  new_dat <- traitsDF %>%
    dplyr::filter(Species==species_name)

  predicted_ESW <- predict.gbm(gbm1, new_dat,
                               n.trees=gbm1$gbm.call$best.trees,
                               type="response")

  loo_summary <- data.frame(Species=species_name,
                            observed_ESW=new_dat$ESW,
                            predicted_ESW=predicted_ESW)

  return(loo_summary)

}

# the following code will take over an hour to run
leave_one_out_results_reduced <- bind_rows(lapply(unique(traitsDF$Species),
                                           leave_one_out_function))
saveRDS(leave_one_out_results_reduced, file="outputs/leave_one_out_results_brt_reduced.rds")

## Compare trait predictions and original estimates -----------------------------------------------

leave_one_out_results_reduced <- readRDS("outputs/leave_one_out_results_brt_reduced.rds")

DescTools::CCC(leave_one_out_results_reduced$observed_ESW,
               leave_one_out_results_reduced$predicted_ESW)$rho.c
#est    lwr.ci   upr.ci
#1 0.6204263 0.4765518 0.731935

DescTools::CCC(leave_one_out_results_reduced$observed_ESW,
               leave_one_out_results_reduced$predicted_ESW)$C.b
# 0.9305272

summary(abs(leave_one_out_results_reduced$observed_ESW-leave_one_out_results_reduced$predicted_ESW))
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.01174  3.96547  9.10580 10.06412 13.14698 39.01209

cor.test(leave_one_out_results_reduced$observed_ESW,
         leave_one_out_results_reduced$predicted_ESW)

# Pearson's product-moment correlation
# 
# data:  leave_one_out_results_reduced$observed_ESW and leave_one_out_results_reduced$predicted_ESW
# t = 7.4849, df = 70, p-value = 1.619e-10
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.5145596 0.7782100
# sample estimates:
#       cor 
# 0.6667471

#Figure 3
Fig3d <- ggplot(leave_one_out_results_reduced, aes(x=observed_ESW/100, y=predicted_ESW/100))+
  geom_point()+
  theme_minimal()+
  theme(axis.text=element_text(color="black"),
        axis.line.x = element_line(color = "black", size = 0.5),
        axis.line.y = element_line(color = "black", size = 0.5))+
  xlab("Detection probability")+
  ylab("Trait-based detection probability")+
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  xlim(0.15,1) + ylim(0.15,1)+
  annotate("text", x=0.40, y=0.95, label= "rho = 0.67, MAD = 10.06 m")
  #geom_smooth(method="lm")
Fig3d

## combine Fig 3 -------------------------------------------------------

Fig3a <- readRDS("plots/Figure_3a.rds") # script 02
Fig3b <- readRDS("plots/Figure_3b") # script 02

g1 <- Fig3a + Fig3b + plot_layout(axes = 'collect') 
g3 <- Fig3c + Fig3d + plot_layout(axes = 'collect')  

g1 / g3 + plot_layout(axes = 'collect') +
  plot_annotation(tag_levels = list(c("a", "b", "c", "d")))
  
ggsave("plots/Figure_3.png", width=8, height=7, units="in")

# Population size estimate --------------------------------------------

#We now want to understand the implications of using the predicted detection probability 
#instead of direct estimates based on observed data 

constantESW <- readRDS("outputs/distance_passerines_constant.rds") # see script 01
bootsESW <- readRDS("outputs/distance_passerines_constant_boot.rds") #see HPC_Detection_bootstraps (used to calc CIs)
leave_one_out_results <- readRDS("outputs/leave_one_out_results_brt.rds") # made above

# correct all the site level predictions by detection and sum abundances across all sites (= total pop size)
popOutput <- data.frame()
allspecies <- sort(unique(constantESW$Species))
 
for(s in 1:length(allspecies)){

  myspecies <- allspecies[s]

  #get site-level relative abundance predictions (not corrected for detection)
  sitePreds <- list.files("../data/models", full.names = TRUE) %>% # site-level preds - seee HPC_Abundance_models
    str_subset("xggridpreds") %>%
    str_subset(strsplit(myspecies,"/")[[1]][1]) %>%
    readRDS()

  #get esws
  constantSp <- constantESW %>% filter(Species == myspecies)
  bootSp <- bootsESW %>% filter(Species == myspecies)

  #first use preds column x median detection for species
  meanPreds <- sitePreds %>%
    dplyr::select(species:preds) %>%
    mutate(terr_cover = 1- buffer_seawater) %>%
    add_column(ESW = constantSp$ESW) %>%
    mutate(Pop_size = (preds * terr_cover)/(ESW/100))

  #add on trait estimate
  meanPreds$predESW <- leave_one_out_results$predicted_ESW[match(meanPreds$species,
                                                                 leave_one_out_results$Species)]
  meanPreds <- meanPreds %>%
    mutate(Trait_pop_size = (preds * terr_cover)/(predESW/100))

  #also get range of predictions using bootstrapped sampled - sum each one
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
                              Obs_Ab_Abund = sum(meanPreds$Pop_size),
                              Preds_Ab_Abund = sum(meanPreds$Trait_pop_size),
                              Obs_Lower_Abund = as.numeric(quantile(bootPop, 0.025)),
                              Obs_Upper_Abund = as.numeric(quantile(bootPop, 0.975)))

  popOutput <- bind_rows(popOutput, popOutput_new)

}

saveRDS(popOutput, file="outputs/popOutput_brt.rds")

# Compare original and trait-based, like in script 02

popOutput <- readRDS("outputs/popOutput_brt.rds") %>%
  inner_join(.,leave_one_out_results, by="Species") %>%
  mutate(Pop_error = (Obs_Ab_Abund - Preds_Ab_Abund)/Obs_Ab_Abund,
         Pop_u_error = (Obs_Ab_Abund - Obs_Lower_Abund)/Obs_Ab_Abund, 
         Pop_l_error = (Obs_Ab_Abund - Obs_Upper_Abund)/Obs_Ab_Abund, 
         ESW_error = observed_ESW - predicted_ESW,
         P_error = ESW_error/100) 

#for Table S6:
write.csv(popOutput, file = "outputs/popestimates_brt.csv", row.names=FALSE)

#for how many species is the trait estimates within the bootstraped CI
mean(popOutput$Obs_Upper_Abund > popOutput$Preds_Ab_Abund &
       popOutput$Obs_Lower_Abund < popOutput$Preds_Ab_Abund)
#75%

(outSpecies <- popOutput$Species[!(popOutput$Obs_Upper_Abund > popOutput$Preds_Ab_Abund &
                                     popOutput$Obs_Lower_Abund < popOutput$Preds_Ab_Abund)])

g_error <- ggplot(popOutput) +
  geom_point(aes(x = ESW_error/100, y = Pop_error*100)) +
  geom_hline(yintercept =  median(popOutput$Pop_u_error)*100, linetype ="dashed", color = "red") + # null error
  geom_hline(yintercept =  median(popOutput$Pop_l_error)*100, linetype ="dashed", color = "red") + # null error
  geom_hline(yintercept = 0, linetype ="dashed") +
  geom_vline(xintercept = 0, linetype ="dashed") +
  xlab("Detection error") +
  ylab("Population error (%)")

g_Error <- ggMarginal(g_error, type="density", fill = "lightblue")

g_abs <- ggplot(popOutput) +
  geom_point(aes(x = abs(ESW_error/100), y = abs(Pop_error*100))) +
  geom_hline(yintercept =  median(abs(c(popOutput$Pop_u_error, popOutput$Pop_l_error))*100), linetype ="dashed", color = "red") + # null error
  xlab("|Detection error|") +
  ylab("|Population error (%)|")

g_Abs <- ggMarginal(g_abs, type="boxplot", fill = "lightblue")

lm_error_abs <- lm(abs(Pop_error) ~ abs(P_error), data=popOutput)
summary(lm_error_abs)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.026936   0.007988   3.372  0.00122 ** 
#  abs(P_error) 1.535965   0.060630  25.334  < 2e-16 ***
  
confint(lm_error_abs)

summary(abs(popOutput$Pop_error))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00540 0.08189 0.16587 0.18712 0.26916 0.54932 


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

DescTools::CCC(popOutput$Obs_Ab_Abund,popOutput$Preds_Ab_Abund)$rho.c # 0.97
DescTools::CCC(popOutput$Obs_Ab_Abund,popOutput$Preds_Ab_Abund)$C.b # 0.99
cor(popOutput$Obs_Ab_Abund,popOutput$Preds_Ab_Abund) # 0.97

# saving alone
patchwork::wrap_elements(gScatter)/((patchwork::wrap_elements(g_Error) +
                                                          patchwork::wrap_elements(g_Abs))) + 
  plot_annotation(tag_levels = list(c("a", "b", "c"))) 

ggsave("plots/Figure_S4.png",width = 6, height = 6)

# end --------------------------------------------------------------------------