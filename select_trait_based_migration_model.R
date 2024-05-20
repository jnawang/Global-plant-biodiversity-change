library(MuMIn) #version 1.43.6
library(lme4) #version 1.1.21
library(optimx) #version 2018-7.10
library(jtools) #version 2.0.1
library(tidyverse, caret)
#
rm(list=ls()) 

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#############################
#### DATASET preparation ####
#############################
rSdata = read.table("../Table_S1.csv",sep=";",h=T,dec=".",stringsAsFactors = FALSE) #exported from Table_S1.xlsx with separator=";" and decimal="."
rSdata = rSdata %>% filter(Kingdom == 'Plantae' & Ecosystem == 'Terrestrial')
rSdata = rSdata[order(rSdata$n_sp),]     # order the data based on n_sp, why?  Junna

# Attach barrier information
# add barrier information
barrier <- read.csv('../bioshift_HFI_mountain.csv')
barrier$high_mountainT <- barrier$high_mountain + barrier$high_mountain_scat
barrier$mountainT <- barrier$high_mountainT + barrier$low_mountain + barrier$low_mountain_scat
barrier <- barrier[,-1]
rSdata  <- rSdata %>% left_join(barrier, by = c('Source' = 'name'))

# add growth form, dispersal syndrome, and other traits; who knows if these traits are useful?
plants_gf   <- read.csv('../../traits_phylo/plants_gf.csv')
plants_disp <- read.csv('../../traits_phylo/plants_disp.csv')
rSdata      <- rSdata %>% left_join(plants_gf, by = c('Species' = 'species')) %>% left_join(plants_disp, by = c('Species' = 'species')) %>% filter(!is.na(disp))     # I am not sure if I should delete the 457 species [either cannot be binded to gene tree, or name has error]?

# add plant height and seed mass
traits      <- read.csv('../../traits_phylo/imputed_traits.csv')
rSdata      <- rSdata %>% left_join(traits, by = c('Species' = 'species'))


# species have both latitudinal and elevational shifts
species_both <- intersect(rSdata$Species[is.na(rSdata$LatVeloT)], rSdata$Species[is.na(rSdata$EleVeloT)])   # 748 species
rSdata$inMount[!is.na(rSdata$EleVeloT)] = 'YES'
rSdata$inMount[!is.na(rSdata$LatVeloT)] = 'NO'
rSdata$inMount[rSdata$Species %in% species_both] = 'BOTH'

#transforming continuous method variables in qualitative variables    # I see, learned, this is one way to add random effect; Junna
rSdata$NtaxaF = as.numeric(cut(rSdata$Ntaxa,quantile(rSdata$Ntaxa,c(0,0.25,0.5,0.65,1)),include.lowest = T))  #75th quantile is also the max
rSdata$StartF = as.numeric(cut(rSdata$Start,quantile(rSdata$Start,c(0,0.25,0.5,0.75,1)),include.lowest = T))  #75th quantile is also the max
rSdata$AreaF = as.numeric(cut(rSdata$Area,quantile(rSdata$Area,c(0,0.25,0.5,0.75,1)),include.lowest = T))  #75th quantile is also the max

#relevel variable to limit singularity issue; got the meaning!
rSdata$Sampling = ifelse(rSdata$Sampling == "TWO","TWO","MULT")
rSdata$Quality = ifelse(rSdata$Quality == "BALANCED","RESURVEYED",rSdata$Quality)

# variable pre-selection     # list all the variables to be selected, Junna
chosen_varlatT=c("ShiftR", "LatVeloT", "temp.mean.lat", "HFI.lat", "temp.var.lat", "prec.mean.lat", "prec.var.lat", "elev.het.lat", "high_mountainT", "low_mountain", "mountainT", "SeedMass_log2", "Height", "inMount", "NtaxaF", "StartF", "AreaF", "PrAb", "Sampling", "Grain", "Quality", "Signif", "Source", "Class", "Family", "Genus", "Species", "gf", "disp", "Gradient", "Ecosystem", "LifeForm")
chosen_varele =c("ShiftR", "EleVeloT", "temp.mean.ele", "HFI.ele", "temp.var.ele", "prec.mean.ele", "prec.var.ele", "elev.het.ele", "high_mountainT", "low_mountain", "mountainT", "SeedMass_log2", "Height", "inMount", "NtaxaF", "StartF", "AreaF", "PrAb", "Sampling", "Grain", "Quality", "Signif", "Source", "Class", "Family", "Genus", "Species", "gf", "disp", "Gradient", "Ecosystem", "LifeForm")

#####################################################
#### MODELING: random effect structure selection ####
#####################################################
#Elevational range shift; I have to increase my working efficiency! Use as many packages as possible; Junna
Class <- as.data.frame(table(rSdata$Class[rSdata$Gradient == "Elevation"]))
names(Class) <- c("Class", "Freq")
Class


# Selecting observation and variables to analyse range shifts
data <- rSdata[rSdata$Gradient == "Elevation", chosen_varele]
#data <- data[which(data$Class %in% unique(Class$Class)), ]
data$Class  <- as.factor(as.character(data$Class))
data$Family <- as.factor(as.character(data$Family))
data$Genus  <- as.factor(as.character(data$Genus))
data$LifeForm <- as.factor(as.character(data$LifeForm))
data$StartF=as.factor(data$StartF)
data$NtaxaF=as.factor(data$NtaxaF)
data$AreaF=as.factor(data$AreaF)
data$inMount = as.factor(data$inMount)
data$gf      = as.factor(data$gf)
data$disp    = as.factor(data$disp)
data <- na.omit(data)
data=droplevels(data)
dim(data) # 8880 obs

#scaling explaining continuous variables
data$g_velT=gscale(data$EleVeloT)
data$BaselineT2=gscale(data$temp.mean.ele)^2
data$g_basT=gscale(data$temp.mean.ele)
data$g_varT=gscale(data$temp.var.ele)
data$g_basP=gscale(data$prec.mean.ele)
data$g_varP=gscale(data$prec.var.ele)
data$g_hetE=gscale(data$elev.het.ele)
data$g_HFI=gscale(data$HFI.ele)
data$g_height=gscale(data$Height)
data$g_seedmass=gscale(data$SeedMass_log2)
#
data$g_highMount=gscale(data$high_mountainT)
data$g_lowMount =gscale(data$low_mountain)
data$g_totMount =gscale(data$mountainT)

# indicate which mountain data to use
data$g_Mount = data$g_totMount

# not very sure why it is necessary to do this! Junna
t1=data.frame(table(data$LifeForm))
t1=t1[order(t1$Freq,decreasing=T),]
ref1=as.character(t1[1,1])
data$LifeForm=relevel(data$LifeForm,ref=ref1) #Intercept value of the model will be the range shift estimation of the life form having the greatest number of observation

# this is a key function; 
testSingAl=function(x,f,data,n=1:length(x)){
  #function testing for random effect structure
  #x=set of variables to test as random effect
  #f=formula of the fixed structure of the model
  #data= data base with all variables
  #n= levels of variable combination to test (start from 1 when signle variable are tested to the size of the set of variables to test in combination)
  b=1
  for(i in n){
    v1=combn(x,i)
    if(i>1){
      v2=apply(t(v1),1,paste,sep="",collapse="+")
    }else{
      v2=t(v1)[,1]
    }
    for(j in 1:length(v2)){
      print(paste(i," : ",j,sep=""))
      f1=paste(f,"+",v2[j],sep="")
      lme1=lmer(as.formula(f1),data,REML=F,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      a=summary(lme1)$optinfo$conv$lme4$messages
      if(is.null(a)==F){
        a=paste(a,sep="",collapse="_")
      }else{
        a=NA
      }
      res=data.frame(test=v2[j],R2m=r.squaredGLMM(lme1)[[1]],R2c=r.squaredGLMM(lme1)[[2]],aic=AIC(lme1)[[1]],aicc=AICc(lme1)[[1]],singular=isSingular(lme1)[[1]],nV=i,source=grep("Source",v2[j])[1],warning=a[[1]])
      if(b==1){
        resOK=res
      }else{
        resOK=rbind(resOK,res)
      }
      b=b+1
    }
  }
  return(resOK)
  #output:
  #test= random structure
  #R2m= marginal R2 of the model
  #R2c= conditional R2 of the model
  #aic= Akaike information criterion of the model
  #aicc= Akaike information criterion with small-sample correction
  #singular= output of the model singularity test (FALSE= no signularity; TRUE= singularity issue)
  #nV= number of variables tested as random effect
  #warning= inform for warning during the model fit (such as convergence issue)
}

#selection of the set of qualitative variables to test as random effect (based on number of observations, and correlation among variables); unclear Junna
table(data$AreaF)
table(data$StartF)
table(data$NtaxaF)
table(data$PrAb)
table(data$Sampling)  
table(data$Quality)
table(data$Signif)
table(data$Grain)
table(data$gf)
table(data$disp)

#selection of the best random effect structure
x=c("(1|AreaF)","(1|StartF)","(1|NtaxaF)","(1|PrAb)","(1|Sampling)","(1|Grain)","(1|Quality)","(1|inMount)", "(1|gf)", "(1|disp)")  #remove: "(1|Signif)", 12/6/2022
f="ShiftR ~ BaselineT2 + g_Mount + g_varT + g_hetE + g_height + g_seedmass + g_velT * g_basT + g_velT * g_HFI"
tx=testSingAl(x,f,data,n=1:length(x)) #take some time
write.table(tx,"sing_ELE.csv",sep=",",dec=".",row=F)
tx$Issue=grepl("failed to converge",as.character(tx$warning))
tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue; the meaning of is.na(source)? Junna
tx1=tx1[order(tx1$aicc,decreasing=F),] 

if(nrow(tx1)==0){ #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|Family/Genus) + (1|Source)
  f1=paste(f,"+(1|Source)",sep="")
  mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
}else{
  z=1
  iS=T
  while(iS==T & z<=nrow(tx1)){
    f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))   # why do this calculation again, because these formula are already not irregular? Junna. This is commonly used. 
    iS=isSingular(mod1)[[1]]
    z=z+1
  }
}


### sources and genus may cause singularity!
if(isSingular(mod1)==T){ #in case of the selected model (with REML=T) is singular, a new random model structure selection is conducted without any phylogenetic effect
  f="ShiftR ~ ShiftR ~ BaselineT2 + g_Mount + g_varT + g_hetE + g_height + g_seedmass + g_velT * g_basT + g_velT * g_HFI"
  tx=testSingAl(x,f,data,n=1:length(x)) #take some time
  write.table(tx,"sing_ELE.csv",sep=",",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),]
  if(nrow(tx1)>0){
    z=1
    iS=T
    while(iS==T & z<=nrow(tx1)){
      f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      iS=isSingular(mod1)[[1]]
      z=z+1
    }
  }
}

if(isSingular(mod1)==F){
  mod2=lmer(formula(mod1),data,REML=F,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
  ms1=dredge(mod2) #take some time: ms1 is the table reporting results of the model selection through AICc
  write.table(ms1,"elevation_selAIC.csv",sep=",",dec=".",row=F)
  ms1E=ms1
  
  #best model selection with no singularity
  a=1
  iS=T
  while(iS==T){
    print(a)
    m1=get.models(ms1E,subset=a)[[1]]
    mbe=lmer(formula(m1), data=data,REML=T,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    iS=isSingular(mbe)[[1]]
    a=a+1
  }
  summary(mbe)
  r.squaredGLMM(mbe)
}

#Latitudinal range shifts
Class <- as.data.frame(table(rSdata$Class[rSdata$Gradient == "Latitudinal"]))
names(Class) <- c("Class", "Freq")
Class=Class[which(Class$Freq>30), ]       # Criteria: Class > 30 obs
Class


# Selecting observation and variables to analyse range shifts
data <- rSdata[rSdata$Gradient == "Latitudinal", chosen_varlatT]
data <- data[which(data$Class %in% unique(Class$Class)), ]
data$Class <- as.factor(as.character(data$Class))
data$Family <- as.factor(as.character(data$Family))
data$Genus <- as.factor(as.character(data$Genus))
data$LifeForm <- as.factor(as.character(data$LifeForm))
data$StartF=as.factor(data$StartF)
data$NtaxaF=as.factor(data$NtaxaF)
data$AreaF=as.factor(data$AreaF)
data$inMount = as.factor(data$inMount)
data$gf      = as.factor(data$gf)
data$disp    = as.factor(data$disp)
data <- na.omit(data)
data=droplevels(data)
dim(data) # 15118 obs

#scaling explaining continuous variables
data$g_velT=gscale(data$LatVeloT)
data$BaselineT2=gscale(data$temp.mean.lat)^2
data$g_basT=gscale(data$temp.mean.lat)
data$g_varT=gscale(data$temp.var.lat)
data$g_basP=gscale(data$prec.mean.lat)
data$g_varP=gscale(data$prec.var.lat)
data$g_hetE=gscale(data$elev.het.lat)
data$g_HFI=gscale(data$HFI.lat)
data$g_height=gscale(data$Height)
data$g_seedmass=gscale(data$SeedMass_log2)
#
data$g_highMount=gscale(data$high_mountainT)
data$g_lowMount =gscale(data$low_mountain)
data$g_totMount =gscale(data$mountainT)

# indicate which mountain data to use
data$g_Mount = data$g_totMount

t1=data.frame(table(data$LifeForm))
t1=t1[order(t1$Freq,decreasing=T),]
ref1=as.character(t1[1,1])
data$LifeForm=relevel(data$LifeForm,ref=ref1) #Intercept value of the model will be the range shift estimation of the life form having the greatest number of observation

#selection of the set of qualitative variables to test as random effect (based on number of observations, and correlation among variables)
table(data$AreaF)
table(data$StartF)
table(data$NtaxaF)
table(data$PrAb)
table(data$Sampling)  
table(data$Quality)
table(data$Signif)
table(data$Grain)

#selection of the best random effect structure
x=c("(1|AreaF)","(1|StartF)","(1|NtaxaF)","(1|PrAb)","(1|Sampling)","(1|Grain)","(1|Quality)","(1|inMount)", "(1|gf)", "(1|disp)") # set of random effects; remove: "(1|Signif)", 12/6/2022
f="ShiftR ~ BaselineT2 + g_Mount + g_varT + g_height + g_seedmass + g_velT * g_basT + g_velT * g_HFI"
tx=testSingAl(x,f,data,n=1:length(x)) #take some time
write.table(tx,"sing_LATM.csv",sep=",",dec=".",row=F)
tx$Issue=grepl("failed to converge",as.character(tx$warning))
tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
tx1=tx1[order(tx1$aicc,decreasing=F),] 

if(nrow(tx1)==0){ #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|Family/Genus) + (1|Source)
  f1=paste(f,"+(1|Source)",sep="")
  mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
}else{
  z=1
  iS=T
  while(iS==T & z<=nrow(tx1)){
    f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    iS=isSingular(mod1)[[1]]
    z=z+1
  }
}

### sources and genus may cause singularity!
if(isSingular(mod1)==T){ #in case of the selected model (with REML=T) is singular, a new random model structure selection is conducted without any phylogenetic effect
  f="ShiftR ~ ShiftR ~ BaselineT2 + g_Mount + g_varT + g_height + g_seedmass + g_velT * g_basT + g_velT * g_HFI"
  tx=testSingAl(x,f,data,n=1:length(x)) #take some time
  write.table(tx,"sing_LATM.csv",sep=",",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),]
  if(nrow(tx1)>0){
    z=1
    iS=T
    while(iS==T & z<=nrow(tx1)){
      f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      iS=isSingular(mod1)[[1]]
      z=z+1
    }
  }
}


if(isSingular(mod1)==F){
  mod2=lmer(formula(mod1),data,REML=F,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
  ms1=dredge(mod2) #take some time: ms1 is the table reporting results of the model selection through AICc
  write.table(ms1,"latT_selAIC.csv",sep=",",dec=".",row=F)
  ms1Lm=ms1
  
  #best model selection with no singularity
  a=1
  iS=T
  while(iS==T){
    print(a)
    m1=get.models(ms1Lm,subset=a)[[1]]
    mblm=lmer(formula(m1), data=data,REML=T,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    iS=isSingular(mblm)[[1]]
    a=a+1
  }
  summary(mblm)
  r.squaredGLMM(mblm)
}
