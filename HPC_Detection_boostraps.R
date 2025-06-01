# README:
# Script for HPC to run models to boostrap the detection probability estimates

library(tidyverse)
library(gbm)
library(MASS)
library(boot)
library(Distance)

# read in data files --------------------------------

myfolder <- "/data/idiv_ess/DOF/"

data <- readRDS(paste0(myfolder, "data.rds")) # bird data
info <- readRDS(paste0(myfolder, "info.rds")) # survey data

# species list -------------------------------------

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
             "Lanius collurio","Larus canus","Linaria cannabina","Lophophanes cristatus",
             "Loxia curvirostra","Luscinia luscinia","Motacilla alba","Motacilla flava",
             "Muscicapa striata","Oenanthe oenanthe","Parus major","Passer domesticus",
             "Passer montanus","Perdix perdix","Periparus ater","Phoenicurus ochruros",
             "Phoenicurus phoenicurus","Phylloscopus collybita","Phylloscopus sibilatrix",
             "Phylloscopus trochilus","Pica pica","Poecile palustris","Prunella modularis", 
             "Pyrrhula pyrrhula","Regulus regulus","Saxicola rubetra","Sitta europaea",
             "Spinus spinus","Streptopelia decaocto","Sturnus vulgaris","Sylvia atricapilla",
             "Sylvia borin","Troglodytes troglodytes","Turdus merula","Turdus philomelos",
             "Turdus pilaris","Turdus viscivorus","Vanellus vanellus")

# species of this run -----------------------------

task.id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

print(task.id)

myspecies <- species[task.id]
  
print(myspecies)

# get data ---------------------------------------
  
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
  
# Resamples ----------------------------------  
  
  vals <- as.numeric()
  
  for(i in (length(vals)+1):1000){
    
    d <-   detectionData[sample(nrow(detectionData), replace=TRUE),]
    
    hr.model <- ds(d, truncation=100,
                   transect="line",
                   formula = ~1,
                   cutpoints = c(0,25,50,100),
                   key = "hn")
    
    temp <- summary(hr.model)
    vals[i] <- temp$ds$average.p * 100
    
  }
  
  temp <- data.frame(Species = myspecies,
                     vals = vals)
  
saveRDS(temp, file=paste0("ESW_constant_bootstrapsamples_",myspecies,".rds")) 

# # Post HPC processing to combine outputs across species ----------------------
# 
# distanceSp <- list.files("../DOF_esw_bootstraps", full.names = TRUE) %>%
#   map_dfr(readRDS)
# 
# saveRDS(distanceSp, file="outputs/distance_passerines_constant_boot.rds")
