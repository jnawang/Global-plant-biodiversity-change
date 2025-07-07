#####predict latitudinal and longitudinal migration with lmer
library(tidyverse)
library(caret)
library(ggplot2)
library(xgboost)
library(Matrix)
library(data.table)
library(vcd)
library(jtools)
library('lme4')
library(MuMIn)
library(jtools)
rm(list=ls())

# data.table::setDTthreads(2)
setwd('/Users/junnawang/UCDLab/Biodiversity')
# mean, sd for different growth forms, are they significantly different?
data <- read.csv('Bioshift/Lenoir_et_al/Analysis/Table_S1.csv', sep=";")
data_plants <- data %>% filter(Kingdom == 'Plantae' & Ecosystem == 'Terrestrial')   # we have 80 studies. 

# get growth form of these species. Still, so many species we do not know their growth form
plants_gf   <- read.csv('pre_post_processing/traits_phylo/plants_gf.csv')
plants_disp <- read.csv('pre_post_processing/traits_phylo/plants_disp.csv')
data_plants <- data_plants %>% left_join(plants_gf, by = c('Species' = 'species')) %>% left_join(plants_disp, by = c('Species' = 'species')) %>% filter(!is.na(disp))     # I am not sure if I should delete the 457 species [either cannot be binded to gene tree, or name has error]?

# add plant height and seed mass
traits      <- read.csv('pre_post_processing/traits_phylo/imputed_traits.csv')
data_plants <- data_plants %>% left_join(traits, by = c('Species' = 'species'))

# add barrier information
barrier <- read.csv('pre_post_processing/02_predict_shift_rate/02_02_select_best_shift_rate_prediction_model/bioshift_HFI_mountain.csv')        # Need to change A66_P1; 
barrier$high_mountainT <- barrier$high_mountain + barrier$high_mountain_scat
barrier$mountainT <- barrier$high_mountainT + barrier$low_mountain + barrier$low_mountain_scat
barrier <- barrier[,-1]
rSdata  <- data_plants %>% left_join(barrier, by = c('Source' = 'name'))

# factorize categorical variables: areaF, sampling, startF
#transforming continuous method variables in qualitative variables    # I see, learned, this is one way to add random effect; Junna
rSdata$NtaxaF = as.numeric(cut(rSdata$Ntaxa,quantile(rSdata$Ntaxa,c(0,0.25,0.5,0.65,1)),include.lowest = T))  #75th quantile is also the max
# threshold of cut-up: 1, 183, 631, 1334, 4426
rSdata$StartF = as.numeric(cut(rSdata$Start,quantile(rSdata$Start,c(0,0.25,0.5,0.75,1)),include.lowest = T))  #75th quantile is also the max
# threshold of cut-up: 1802,1895,1930,1978,2001
rSdata$AreaF = as.numeric(cut(rSdata$areaT,quantile(rSdata$areaT,c(0,0.25,0.5,0.75,1)),include.lowest = T))  #75th quantile is also the max
# threshold of cut-up: 0.56524, 21673.87, 151136.0, 409778.0, 3090540.0

# relevant variable to limit singularity issue
rSdata$Sampling = ifelse(rSdata$Sampling == "TWO","TWO","MULT")
rSdata$Quality = ifelse(rSdata$Quality == "BALANCED","RESURVEYED",rSdata$Quality)

####################################################################Fit lat models###################################################
#####################################################################################################################################
# obtain latitudinal shift model first
data_lat <- rSdata %>% filter(!is.na(LatVeloT)) # %>% filter(!Species %in% species_both)    # One LatVeloT is missing!
# find modeled species in the lat bioshift dataset.
plants_modeled   <- read.csv('pre_post_processing/traits_phylo/plant_taxonomy_nloc.csv')
sps.bioshift.lat <- intersect(plants_modeled$sps, data_lat$Species)      # 1107 species

# I need to standarize these variables first; meaning of gscale = (x-mean)/sd/2
data_lat$g_ShiftR    <- gscale(data_lat$ShiftR)              #mean: -0.07615671; std: 1.306827;
data_lat$g_LatveloT  <- gscale(data_lat$LatVeloT)            #mean: 0.7524375   ; std: 0.8338242;
data_lat$g_BaseT     <- gscale(data_lat$temp.mean.lat)       #mean: 9.618638    ; std: 4.184582;
data_lat$g_BaseT2    <- gscale(data_lat$temp.mean.lat)^2     #mean: 9.618638    ; std: 4.184582;  The same as above;
data_lat$g_varT      <- gscale(data_lat$temp.var.lat)        #mean: 552.77     ; std: 257.17;
data_lat$g_HFI       <- gscale(data_lat$HFI.lat)             #mean: 19.14025   ; std: 5.871376;
data_lat$g_Height    <- gscale(data_lat$Height)              #mean: 2.87087    ; std: 6.263124;
data_lat$g_Seedmass  <- gscale(data_lat$SeedMass_log2)       #mean: 0.689      ; std: 3.684308;
data_lat$g_Mount     <- gscale(data_lat$mountainT)           #mean: 0.3383832  ; std: 0.2572103;

# make categorical variables to be factors
data_lat$StartF  <- as.factor(data_lat$StartF)
data_lat$NtaxaF  <- as.factor(data_lat$NtaxaF)
data_lat$AreaF   <- as.factor(data_lat$AreaF)
data_lat$Signif  <- as.factor(data_lat$Signif)
data_lat$Species <- as.factor(data_lat$Species)
data_lat$disp    <- as.factor(data_lat$disp)

# train the model first: train a species_model [use species identity] and a trait_model [use dispersal syndrome].
# species_level model first;
mod_lat_species = lmer(g_ShiftR ~ g_Height + g_BaseT2 + g_Mount + g_varT + g_LatveloT * g_BaseT + g_LatveloT * g_HFI + (1|AreaF) + (1|Species), data_lat, REML=TRUE, na.action="na.fail", control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
summary(mod_lat_species)
# How to know significance of each variable? check functions in the package
confint(mod_lat_species, method = "Wald", level = 0.95)   # c("profile", "Wald", "boot")
# Check effect size of each random variable, including random effects.
r.squaredGLMM(mod_lat_species)    # R2m: 0.07437664; R2c: 0.2023812

#####alternative model start?
# mod_lat_species = lmer(g_ShiftR ~ g_varT + g_LatveloT * g_BaseT + g_LatveloT * g_HFI  + (1|AreaF) + (1|Species) , data_lat, REML=TRUE, na.action="na.fail", control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
# summary(mod_lat_species)
# # How to know significance of each variable? check functions in the package
# confint(mod_lat_species, method = "Wald", level = 0.95)   # c("profile", "Wald", "boot")
# # Check effect size of each random variable, including random effects.
# r.squaredGLMM(mod_lat_species)    # R2m: 0.1203691 ; R2c: 0.3111541
#######
#######
# trait_level model next;
mod_lat_trait = lmer(g_ShiftR ~ g_BaseT2 + g_Mount + g_varT + g_LatveloT * g_BaseT + g_LatveloT * g_HFI + (1|AreaF) + (1|disp), data_lat, REML=TRUE, na.action="na.fail", control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
summary(mod_lat_trait)
# How to know significance of each variable? check functions in the package
confint(mod_lat_trait, method = "Wald", level = 0.95)   # c("profile", "Wald", "boot")
# Check effect size of each random variable, including random effects.
r.squaredGLMM(mod_lat_trait)    # R2m: 0.0732952; R2c: 0.1508451


#####alternative model start?
# mod_lat_trait = lmer(g_ShiftR ~ g_varT + g_LatveloT * g_BaseT + g_LatveloT * g_HFI + (1|AreaF) + (1|disp), data_lat, REML=TRUE, na.action="na.fail", control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
# summary(mod_lat_trait)
# # How to know significance of each variable? check functions in the package
# confint(mod_lat_trait, method = "Wald", level = 0.95)   # c("profile", "Wald", "boot")
# # Check effect size of each random variable, including random effects.
# r.squaredGLMM(mod_lat_trait)    # R2m: 0.1162595; R2c: 0.2700121

ranef(mod_lat_trait)
##effect size of random effects (ranef(mod_lat_trait));
##disp: Anemochor   0.083379554; Autochor    0.013394151; Chamaechor -0.005204868; Hemerochor -0.013437598;
##disp: Nautochor   0.020228425; Ombrochor  -0.144021936; Zoochor     0.045662272;
##AreaF: 1 -0.04299852; 2  0.17025040; 3 -0.02013603; 4 -0.10711585; 
###################################################################Fit elev models#########################################################
###########################################################################################################################################
# obtain the elevation shift model next
data_ele <- rSdata %>% filter(!is.na(EleVeloT) & Gradient == "Elevation")
# find modeled species in the ele bioshift dataset.
sps.bioshift.ele <- intersect(plants_modeled$sps, data_ele$Species) # 4546 species
# deal with na HFI
data_ele$HFI.ele[is.na(data_ele$HFI.ele)] <- 0.34            # give a very small value. 

# I need to standarize these variables first; meaning of gscale = (x-mean)/sd/2
data_ele$g_ShiftR    <- gscale(data_ele$ShiftR)              #mean: 1.150135; std: 5.570041;
data_ele$g_EleveloT  <- gscale(data_ele$EleVeloT)            #mean:    ; std: ;
data_ele$g_BaseT     <- gscale(data_ele$temp.mean.ele)       #mean:    ; std: ;
data_ele$g_BaseT2    <- gscale(data_ele$temp.mean.ele)^2     #mean:    ; std: ;  The same as above;
data_ele$g_varT      <- gscale(data_ele$temp.var.ele)        #mean:    ; std: ;
data_ele$g_HFI       <- gscale(data_ele$HFI.ele)             #mean:    ; std: ;
data_ele$g_mount     <- gscale(data_ele$mountainT)           #this is already mountain fraction!
data_ele$g_Height    <- gscale(data_ele$Height)              #mean:     ; std: ;
data_ele$g_Seedmass  <- gscale(data_ele$SeedMass_log2)       #mean:     ; std: ;

# make categorical variables to be factors
data_ele$StartF   <- as.factor(data_ele$StartF)
data_ele$Sampling <- as.factor(data_ele$Sampling)
data_ele$AreaF    <- as.factor(data_ele$AreaF)
data_ele$Quality  <- as.factor(data_ele$Quality)
data_ele$Species  <- as.factor(data_ele$Species)
data_ele$disp     <- as.factor(data_ele$disp)

# train the model first: train a species_model [use species identity] and a trait_model [use dispersal syndrome].
# species_level model first; g_BaseT2 + 
mod_ele_species = lmer(g_ShiftR ~ g_Seedmass + g_varT + g_mount + g_EleveloT * g_BaseT + (1|StartF) + (1|Sampling) + (1|Quality) + (1|Species), data_ele, REML=TRUE, na.action="na.fail", control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
summary(mod_ele_species)
# How to know significance of each variable? check functions in the package
confint(mod_ele_species, method = "Wald", level = 0.95)   # c("profile", "Wald", "boot")
# Check effect size of each random variable, including random effects.
r.squaredGLMM(mod_ele_species)    # R2m: 0.020; R2c: 0.286
##effect size of random effects (ranef(mod_ele_species))
##StartF: 1 -0.07427801; 2 -0.08834376; 3  0.05880334; 4  0.10381843;
##Quality: LOW -0.06129444; MODELED 0.02810256; RESURVEYED  0.03319187;
##Sampling: MULT: 0.1419091; TWO:   -0.1419091;
#######
#######
# trait_level model next;
mod_ele_trait = lmer(g_ShiftR ~ g_varT + g_mount + g_EleveloT * g_BaseT + (1|StartF) + (1|Sampling) + (1|Quality), data_ele, REML=TRUE, na.action="na.fail", control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
summary(mod_ele_trait)
# How to know significance of each variable? check functions in the package
confint(mod_ele_trait, method = "Wald", level = 0.95)   # c("profile", "Wald", "boot")
# Check effect size of each random variable, including random effects.
r.squaredGLMM(mod_ele_trait)    # R2m: 0.021; R2c: 0.20
##effect size of random effects (ranef(mod_ele_trait))
##StartF: 1 -0.07170338; 2 -0.09096267; 3  0.05366448; 4  0.10900158;
##Quality: LOW -0.05519604; MODELED 0.02048866; RESURVEYED  0.03470738;
##Sampling: MULT: 0.1374191; TWO:   -0.1374191;
#######
#######
########################################################################Predict lat shift###################################################
############################################################################################################################################
# read through all the data
data_env <- read.csv('pre_post_processing/02_predict_shift_rate/02_03_perform_shift_rate_prediction/env_data_predict_bioshift2025.csv')
# remove the first column and transpose; no longer used
# data_env <- data_env %>% dplyr::select(-1) %>% t %>% as.data.frame()

# the col names are: area.in.mount, area.on.flat, frac.mount, HFI.ele, HFI.lat, temp.mean.ele, temp.mean.lat, temp.var.ele, temp.var.lat, elev.het.ele, elev.het.lat, EleVeloT (16 col; 12-27), LatVeloT (16 col; 28-43)
colnames(data_env)[1:12] = c('Species', 'area.mount', 'area.flat', 'frac.mount', 'HFI.ele', 'HFI.lat', 'temp.mean.ele', 'temp.mean.lat', 'temp.var.ele', 'temp.var.lat', 'elev.ele', 'elev.lat')
data_env <- data_env %>% left_join(plants_disp, by = c('Species' = 'species')) %>% filter(!is.na(disp)) %>% left_join(traits, by = c('Species' = 'species'))

# separate lat and ele
data_env_lat          <- data_env[, c('Species', 'disp')]
#
# we cannot use this method for AreaF because most of our species have large areas which corresponding the lowest shift rates, use specific values
# data_env_lat$AreaF    <- as.factor(ifelse(data_env$area.mount + data_env$area.flat < 151136.0,
#                                           ifelse(data_env$area.mount + data_env$area.flat < 21673.87, 1, 2),
#                                           ifelse(data_env$area.mount + data_env$area.flat < 409778.0, 3, 4)))
data_env_lat$AreaF    <- as.factor(2)       # M: 1; F: 2; L: 4;
#
data_env_lat$g_Height <- (data_env$Height - mean(data_lat$Height)) / 2 / sd(data_lat$Height)
data_env_lat$g_BaseT  <- (data_env$temp.mean.lat - mean(data_lat$temp.mean.lat)) / 2 / sd(data_lat$temp.mean.lat)
data_env_lat$g_BaseT2 <- data_env_lat$g_BaseT^2
data_env_lat$g_Mount  <- (data_env$frac.mount - mean(data_lat$mountainT)) / 2 / sd(data_lat$mountainT)
data_env_lat$g_varT   <- (data_env$temp.var.lat - mean(data_lat$temp.var.lat)) / 2 / sd(data_lat$temp.var.lat)
data_env_lat$g_HFI    <- (data_env$HFI.lat - mean(data_lat$HFI.lat)) / 2 / sd(data_lat$HFI.lat)

# repeat this 16 times [scenarios]
data_env_lat <- data_env_lat[rep(seq_len(nrow(data_env_lat)), each = 16), ]
data_env_lat$g_LatveloT <- (c(t(as.matrix(data_env[, 29:44]))) - mean(data_lat$LatVeloT)) / 2 / sd(data_lat$LatVeloT)
data_env_lat$scenario <- rep(1:16, nrow(data_env))

# filter out NA rows
data_env_lat_species <- data_env_lat  %>% filter(!is.na(g_HFI) & Species %in% sps.bioshift.lat)
data_env_lat_trait   <- data_env_lat  %>% filter(!is.na(g_HFI) & !is.na(disp) & !(Species %in% sps.bioshift.lat))

# do the prediction for species level and traits level separately, and then merge.
# Still need to do: mean values and different credible values!
data_env_lat_species$LatShiftR <- predict(mod_lat_species, newdata = data_env_lat_species)  
data_env_lat_trait$LatShiftR   <- predict(mod_lat_trait,   newdata = data_env_lat_trait)    

# another prediction method with credible interval
mySumm_lat_species <- function(.) {
  predict(., newdata=data_env_lat_species, re.form=NULL)
}
boot_lat_species <- lme4::bootMer(mod_lat_species, mySumm_lat_species, nsim=200, use.u=TRUE, type="parametric", verbose = TRUE) # use.u whether random effects are simulated/bootstrapped? should I use.u = TRUE? type = "parametric", "semiparametric"
mySumm_lat_trait <- function(.) {
  predict(., newdata=data_env_lat_trait, re.form=NULL)
}
boot_lat_trait <- lme4::bootMer(mod_lat_trait, mySumm_lat_trait, nsim=200, use.u=TRUE, type="parametric", verbose = TRUE) # use.u whether random 

####Collapse bootstrap into median, 95% PI, 50% PI
sumBoot <- function(merBoot) {
  return(
    data.frame(fit = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.5, na.rm=TRUE))),
               lwr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.025, na.rm=TRUE))),
               upr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.975, na.rm=TRUE))),
               lwr0 = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.25, na.rm=TRUE))),
               upr0 = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.75, na.rm=TRUE)))
    )
  )
}
#
data_env_lat_species[,c('Lat50', 'Lat025', 'Lat975', 'Lat25', 'Lat75')] <- sumBoot(boot_lat_species)
data_env_lat_trait[,c('Lat50', 'Lat025', 'Lat975', 'Lat25', 'Lat75')]   <- sumBoot(boot_lat_trait)

# merge
data_env_lat <- rbind(data_env_lat_species, data_env_lat_trait)
# convert to normal values using mean(data_lat$ShiftR) and std()
data_env_lat[, c('LatShiftR', 'Lat50', 'Lat025', 'Lat975', 'Lat25', 'Lat75')] <- data_env_lat[, c('LatShiftR', 'Lat50', 'Lat025', 'Lat975', 'Lat25', 'Lat75')] * 2 * sd(data_lat$ShiftR) + mean(data_lat$ShiftR)

########################################################################Predict elevational shift############################
#############################################################################################################################
#Need to use the area I calculated. 
data_env_ele          <- data_env[, c('Species', 'disp')]
data_env_ele$Sampling <- as.factor("MULT")         # M:MULT; L:TWO; F:MULT
data_env_ele$StartF   <- as.factor(4)              # M: 3; L: 2; F: 4;
data_env_ele$Quality  <- as.factor("RESURVEYED")      # M: MODELED; L: LOW; F: RESURVEYED
# how to choose their values for fast and slow prediction
#
#
data_env_ele$g_Seedmass  <- (data_env$SeedMass_log2 - mean(data_ele$SeedMass_log2)) / 2 / sd(data_ele$SeedMass_log2)
data_env_ele$g_BaseT  <- (data_env$temp.mean.ele - mean(data_ele$temp.mean.ele)) / 2 / sd(data_ele$temp.mean.ele)
data_env_ele$g_BaseT2 <- data_env_ele$g_BaseT^2
data_env_ele$g_varT   <- (data_env$temp.var.ele - mean(data_ele$temp.var.ele)) / 2 / sd(data_ele$temp.var.ele)
data_env_ele$g_HFI    <- (data_env$HFI.ele      - mean(data_ele$HFI.ele)) / 2 / sd(data_ele$HFI.ele)
data_env_ele$g_mount  <- (data_env$frac.mount   - mean(data_ele$mountainT)) / 2 / sd(data_ele$mountainT)
# repeat this 16 times [scenarios]
data_env_ele            <- data_env_ele[rep(seq_len(nrow(data_env_ele)), each = 16), ]
data_env_ele$g_EleveloT <- (c(t(as.matrix(data_env[, 13:28]))) - mean(data_ele$EleVeloT)) / 2 / sd(data_ele$EleVeloT)
data_env_ele$scenario <- rep(1:16, nrow(data_env))

# filter out NA rows
data_env_ele_species <- data_env_ele  %>% filter(!is.na(g_HFI) & Species %in% sps.bioshift.ele)
data_env_ele_trait   <- data_env_ele  %>% filter(!is.na(g_HFI) & !(Species %in% sps.bioshift.ele))

# do the prediction for species level and traits level separately, and then merge.
# Still need to do: mean values and different credible values!
data_env_ele_species$EleShiftR <- predict(mod_ele_species, newdata = data_env_ele_species)
data_env_ele_trait$EleShiftR   <- predict(mod_ele_trait,   newdata = data_env_ele_trait)

# use a new way of doing prediction
mySumm_ele_species <- function(.) {
  predict(., newdata=data_env_ele_species, re.form=NULL)
}
boot_ele_species <- lme4::bootMer(mod_ele_species, mySumm_ele_species, nsim=200, use.u=TRUE, type="parametric", verbose = TRUE) # use.u whether random effects are simulated/bootstrapped? should I use.u = TRUE? type = "parametric", "semiparametric"
mySumm_ele_trait <- function(.) {
  predict(., newdata=data_env_ele_trait, re.form=NULL)
}
boot_ele_trait <- lme4::bootMer(mod_ele_trait, mySumm_ele_trait, nsim=200, use.u=TRUE, type="parametric", verbose = TRUE) # use.u whether random 

#
data_env_ele_species[,c('Ele50', 'Ele025', 'Ele975', 'Ele25', 'Ele75')] <- sumBoot(boot_ele_species)
data_env_ele_trait[,c('Ele50', 'Ele025', 'Ele975', 'Ele25', 'Ele75')]   <- sumBoot(boot_ele_trait)


# merge
data_env_ele <- rbind(data_env_ele_species, data_env_ele_trait)
# convert to normal values using mean(data_lat$ShiftR) and std()
data_env_ele[,c('EleShiftR', 'Ele50', 'Ele025', 'Ele975', 'Ele25', 'Ele75')]  <- data_env_ele[,c('EleShiftR', 'Ele50', 'Ele025', 'Ele975', 'Ele25', 'Ele75')] * 2 * sd(data_ele$ShiftR) + mean(data_ele$ShiftR)

###########################################merge the two data frames and write down the results########################################
data_env_lat <- data_env_lat[, c('Species', 'disp', 'scenario', 'LatShiftR', 'Lat50', 'Lat025', 'Lat975', 'Lat25', 'Lat75')]
data_env_ele <- data_env_ele[, c('Species', 'disp', 'scenario', 'EleShiftR', 'Ele50', 'Ele025', 'Ele975', 'Ele25', 'Ele75')]
ShiftR_Final <- merge(data_env_lat, data_env_ele, by = c('Species', 'disp', 'scenario'), all = TRUE)
write.csv(ShiftR_Final, 'pre_post_processing/02_predict_shift_rate/02_03_perform_shift_rate_prediction/Final_prediction_shift_rate_fast_2025.csv')
