# README:
# Script for HPC to run models to estimate site-level relative abundance estimates

library(tidyverse)
library(gbm)
library(MASS)
library(boot)
library(xgboost)
library(rsample)

# read in data files ---------------------------------------------

myfolder <- "/data/idiv_ess/DOF/"

data <- readRDS(paste0(myfolder, "data.rds")) # bird data
info <- readRDS(paste0(myfolder, "info.rds")) # survey data
gridData <- readRDS(paste0(myfolder, "gridData.rds")) # environ info for all grids
coordsDF <- readRDS(paste0(myfolder, "coordsDF.rds")) # corodinates of survey locations
environData <- readRDS(paste0(myfolder, "environData.rds")) # environ info for sampling locations
# species list ---------------------------------------------------

species <- c("Acanthis flammea/A. cabaret","Accipiter nisus","Acrocephalus palustris",
             "Acrocephalus schoenobaenus","Acrocephalus scirpaceus","Aegithalos caudatus",
             "Alauda arvensis","Anser anser","Anthus pratensis","Anthus trivialis","Buteo buteo",
             "Carduelis carduelis","Certhia familiaris","Chloris chloris",
             "Coccothraustes coccothraustes","Coloeus monedula","Columba livia",
             "Columba oenas","Columba palumbus","Corvus corax","Corvus cornix",
             "Corvus corone","Corvus frugilegus","Curruca communis","Curruca curruca",
             "Cyanistes caeruleus","Dendrocopos major","Emberiza calandra","Emberiza citrinella",
             "Emberiza schoeniclus","Erithacus rubecula","Fringilla coelebs",
             "Gallinula chloropus","Garrulus glandarius","Hippolais icterina",
             "Lanius collurio","Linaria cannabina","Lophophanes cristatus",
             "Loxia curvirostra","Luscinia luscinia","Motacilla alba","Motacilla flava",
             "Muscicapa striata","Oenanthe oenanthe","Parus major","Passer domesticus",
             "Passer montanus","Perdix perdix","Periparus ater","Phoenicurus ochruros",
             "Phoenicurus phoenicurus","Phylloscopus collybita","Phylloscopus sibilatrix",
             "Phylloscopus trochilus","Pica pica","Poecile palustris","Prunella modularis", 
             "Pyrrhula pyrrhula","Regulus regulus","Saxicola rubetra","Sitta europaea",
             "Spinus spinus","Streptopelia decaocto","Sturnus vulgaris","Sylvia atricapilla",
             "Sylvia borin","Troglodytes troglodytes","Turdus merula","Turdus philomelos",
             "Turdus pilaris","Turdus viscivorus","Vanellus vanellus")

# functions ------------------------------------------------------

addCoords <- function(df){
  
  coordsDF <- coordsDF %>%
    as_tibble() %>%
    rename(X = x, Y = y)
  
  df <- inner_join(df, coordsDF, by="kvadratnr")
  
  return(df)
  
}

getSpeciesData <- function(myspecies, data){
  
  #and just one species
  dataS <- data %>% 
    filter(Species==myspecies) 
  
  infoS <- info %>%
    filter(type !="vinter")
  
  #get all transects for the transe
  allData <- left_join(infoS, dataS) %>%
    mutate(total = ifelse(is.na(total), 0, total),
           X.0 = ifelse(is.na(X.0), 0, X.0),
           X.1 = ifelse(is.na(X.1), 0, X.1),
           X.2 = ifelse(is.na(X.2), 0, X.2)) %>%
    mutate(Species = ifelse(is.na(Species), myspecies, Species))
  
  #cap each species count at 99%
  allData <- allData %>%
    group_by(Species) %>%
    mutate(total = ifelse(total > quantile(total,0.99),
                          round(as.numeric(quantile(total, 0.99))),
                          total)) %>%
    ungroup()
  
  #add on environmental data
  allData <- inner_join(allData, environData)
  
  #add coordinates
  allData <- addCoords(allData)
  
  return(allData)
  
}

# species of this run --------------------------------------------

task.id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

print(task.id)

myspecies <- species[task.id]
  
print(myspecies)

# format data -----------------------------------------------

allData <- getSpeciesData(myspecies, data) 
  
allDataS <- allData %>%
    dplyr::select(total, yday, sinceSunrise, 
                  X, Y, starts_with("buffer")) %>%
    dplyr::select(-buffer_mapped, -buffer500_mapped, 
                  -buffer500_total,-buffer_total) #remove irrelevant land cover categories

# develop hypergrid ------------------------------------------------------

#identify best parameters

#make test and training dataset
cust_split <- allDataS  %>%
  initial_split(prop = 0.7)

train <- training(cust_split)
test <- testing(cust_split)

# Set up the grid for hyperparameters
nrounds = c(50, 100, 150)
eta = c(0.01, 0.1, 0.3)
max_depth = c(3, 5, 7)
colsample_bytree = c(0.5, 0.7, 1)
min_child_weight = c(1, 3, 5)

# Initialize a dataframe to store results
results <- data.frame(nrounds = integer(),
                      eta = numeric(),
                      max_depth = numeric(),
                      colsample_bytree = numeric(),
                      min_child_weight = integer(),
                      RMSE = numeric())

# Cross-validate to find the best parameters
for (n in nrounds){
        for (e in eta) {
            for (md in max_depth){
              for (ct in colsample_bytree){
                  for(mw in min_child_weight){
        
        # Run gbm.step with current parameters
        xgbm1 <- xgboost(data = as.matrix(train[,-1]), label = train$total, 
                         subsample = 0.7,
                         max.depth = md, # Lower values avoid over-fitting.
                         min_child_weight = mw, #Larger values avoid over-fitting.
                         eta = e, # Lower values avoid over-fitting.
                         nthread = 3, 
                         nrounds = n, 
                         colsample_bytree = ct, # Lower ratios avoid over-fittin
                         objective = "count:poisson", 
                         booster = 'gbtree')
        
        # Get predictions on test dataset
        predictions <- predict(xgbm1, newdata = as.matrix(test[,-1]))
        
        # Calculate RMSE
        rmse_value <- sqrt(mean((test$total - predictions)^2))
        
        # Store the results
        results <- rbind(results, data.frame(nrounds = n,
                                             eta = e,
                                             max_depth = md,
                                             colsample_bytree = ct,
                                             min_child_weight = mw,
                                             RMSE = rmse_value))
      }
    }
   }
  }
}

# Find the best parameters based on RMSE
best_params <- results[which.min(results$RMSE),]
print(best_params)

# fit final model --------------------------------------------------------

xgbm1 <- xgboost(data = as.matrix(allDataS[,-1]), label = allDataS$total, 
                 subsample = 0.7, # Lower ratios avoid over-fitting.
                 max.depth = best_params$max_depth, # Lower values avoid over-fitting.
                 min_child_weight = best_params$min_child_weight, #Larger values avoid over-fitting.
                 eta = best_params$eta, # Lower values avoid over-fitting.
                 nthread = 3, 
                 nrounds = best_params$nrounds, 
                 colsample_bytree = best_params$colsample_bytree, # Lower ratios avoid over-fittin
                 objective = "count:poisson", 
                 booster = 'gbtree')
  
# prediction data frame --------------------------------------------

#get prediction dataframe for whole country

#best season for the species to predict
bestSeason <- unique(allData$bestSeason[!is.na(allData$bestSeason)])

gridDataS <- gridData %>%
  mutate(island = factor(island)) %>%
  dplyr::select(X, Y, island, starts_with("buffer")) %>%
  dplyr::select(-buffer_mapped, -buffer500_mapped, 
                -buffer500_total,-buffer_total) %>%
  add_column(yday = median(info$yday[info$type==bestSeason]), 
             sinceSunrise = median(info$sinceSunrise[info$type==bestSeason]))

#exclude islands for this species
if(myspecies %in% excludeIslands$Species){
  speciesExclusions <- excludeIslands %>%
    filter(Species == myspecies) %>%
    pull("island")
  
  gridDataS <- gridDataS %>%
    filter(!island %in% speciesExclusions)
}

#remove island as a predictor
gridDataS <- gridDataS %>%
  dplyr::select(-island)

#remove first column
gridDataS <- gridDataS[,names(allDataS)[-1]]

preds<- predict(xgbm1, newdata = as.matrix(gridDataS), type="response")

# bootstrap site predictions -------------------------------------------------

n.boot <- 1000

allPreds <- matrix(data = NA, ncol = n.boot, nrow = nrow(gridDataS))

for(i in 1:n.boot) {
  
  allData <- getSpeciesData(myspecies, data) 
  
  allDataS <- allData %>%
    dplyr::select(total, yday, sinceSunrise, 
                  X, Y, starts_with("buffer")) %>%
    dplyr::select(-buffer_mapped, -buffer500_mapped, 
                  -buffer500_total,-buffer_total)
  
  #sample detections    
  d1 <- allDataS %>% filter(total>0) 
  d1 <-  d1[sample(nrow(d1), replace = TRUE),] 
  
  #and non-detections
  d2 <- allDataS %>% filter(total==0) 
  d2 <-  d2[sample(nrow(d2), replace = TRUE),] 
  
  d <- bind_rows(d1, d2)
  
  #model
  xgbm1 <- xgboost(data = as.matrix(allDataS[,-1]), label = allDataS$total, 
                   subsample = 0.7,
                   max.depth = best_params$max_depth, # Lower values avoid over-fitting.
                   min_child_weight = best_params$min_child_weight, #Larger values avoid over-fitting.
                   eta = best_params$eta, # Lower values avoid over-fitting.
                   nthread = 3, 
                   nrounds = best_params$nrounds, 
                   colsample_bytree = best_params$colsample_bytree, # Lower ratios avoid over-fitting.
                   objective = "count:poisson", 
                   booster = 'gbtree')
  
  allPreds[,i] <- predict(xgbm1, newdata = as.matrix(gridDataS), type="response")
  
}
  
#add predictions
gridDataS <- gridDataS %>%
  add_column(Species = myspecies) %>%
  add_column(preds = preds) %>%
  dplyr::select(Species, X, Y, buffer_seawater, preds) %>%
  bind_cols(., allPreds) %>%
  janitor::clean_names()

gridDataS   %>%
    saveRDS(., file = paste0("xggridpreds_samples_",filename,".rds"))

print('DONE')