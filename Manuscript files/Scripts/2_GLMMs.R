##### Using data outputs from DADA2 pipeline to run mixed models on read abundances 

rm(list=ls())

# packages required  ------------------------------------------------------
# manipulation
library(tidyverse)
library(metagMisc)

#models 
library(MCMCglmm) #option4
library(ggregplot) #variance props model output 

# transformations 
library(ALDEx2)


# read in data  -----------------------------------------------------------

dada16<- read.csv("Manuscript files/Data/modFiltered100DF16.csv") %>% 
  select(-c("X", "x")) 

dada13<- read_csv("Manuscript files/Data/modFiltered100DF13.csv") %>% 
  select(-c(`...1`, "x1"))


#### GLMM --------- 
resp <- "abundance"

fixed16<- c("season")
fixed13<- c("age_class")

mcmc16<- as.formula(paste(resp, "~", paste(fixed16)))
mcmc13<- as.formula(paste(resp, "~", paste(fixed13)))

MicMix16<-list()
MicMix13<-list()

mf <- 20 # multiplying factor for iterations 


Prior2 <- list(R = list(V = diag(1), nu = 0.002),
               G = list(G1 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100),
                        G2 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100)))


Prior3 <- list(R = list(V = diag(1), nu = 0.002),
               G = list(G1 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100),
                        G2 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100),
                        G3 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100)))


Prior4 <- list(R = list(V = diag(1), nu = 0.002),
               G = list(G1 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100),
                        G2 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100),
                        G3 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100),
                        G4 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100)))


# mcmc poisson  -------------------------------------------------------------------

MicMix16[[1]] <- MCMCglmm(fixed= mcmc16, 
                          random= ~  sample + asv + season:asv + id:asv, 
                          prior=Prior4,
                          data=dada16,
                          verbose=T, pr=T, pl=T, # will save posteriors for random effect levels & model predictions for each iteration 
                          family = "poisson",
                          nitt = 13000*mf,
                          thin = 10*mf,burnin=3000*mf)

MicMix13[[1]] <- MCMCglmm(fixed= mcmc13, 
                          random= ~  sample + asv + age_class:asv, 
                          prior=Prior3,
                          data=dada13,
                          verbose=T, pr=T, pl=T, # will save posteriors for random effect levels & model predictions for each iteration 
                          family = "poisson",
                          nitt = 13000*mf,
                          thin = 10*mf,burnin=3000*mf)

MCMCRep(MicMix16[[1]], scale="link") # proportion variance 
MCMCRep(MicMix13[[1]], scale="link") # propoetion variance 



# CLR alternative ---------------------------------------------------------

## clr transformation 

#2016 
df_t16 <- dada16  %>% 
  dplyr::select(asv:abundance) %>% 
  spread(sample, abundance) # spread data to wide format 

df_m16 <-df_t16 %>% 
  column_to_rownames("asv") %>% 
  as.matrix() # turn into a matrix 

#2013 
df_t13 <- dada13  %>% 
  dplyr::select(asv:abundance) %>% 
  spread(sample, abundance) # spread data to wide format 

df_m13 <-df_t13 %>% 
  column_to_rownames("asv") %>% 
  as.matrix() # turn into a matrix 

#2016 
n_sample16<-length(unique(dada16$sample)) #72 
x16<-aldex.clr(df_m16,as.character(seq(1,n_sample16)), mc.samples=128, denom="all", verbose=F) # run ALDEx2

y16<-reshape::melt(getMonteCarloInstances (x16)) # extract the data
colnames(y16)<-c("asv","mcmc_instance","clr","sample") # tidy for df 

z16<-aggregate(y16$clr,by=list(y16$asv,y16$sample),FUN=mean) # take mean across monte-carlo samples
meas_var16<-aggregate(y16$clr,by=list(y16$asv,y16$sample),FUN=var) # calculate variance across monte-carlo samples
z16$sd<-meas_var16$x # add variances to data frame
colnames(z16)<-c("asv","sample","clr","meas_var") # correct dimensions 
mean(z16$clr) #  will be very close to zero

meta16<- dada16 %>% 
  dplyr::select(sample, asv, phylum, abundance, id, season, sex) 

CLR16 <- merge(z16, meta16, by=c("asv", "sample"), all.x=TRUE) #looks good 

#2013 
n_sample13<-length(unique(dada13$sample)) #58 
x13<-aldex.clr(df_m13,as.character(seq(1,n_sample13)), mc.samples=128, denom="all", verbose=F) # run ALDEx2

y13<-reshape::melt(getMonteCarloInstances (x13)) # extract the data
colnames(y13)<-c("asv","mcmc_instance","clr","sample") # tidy for df 

z13<-aggregate(y13$clr,by=list(y13$asv,y13$sample),FUN=mean) # take mean across monte-carlo samples
meas_var13<-aggregate(y13$clr,by=list(y13$asv,y13$sample),FUN=var) # calculate variance across monte-carlo samples
z13$sd<-meas_var13$x # add variances to data frame
colnames(z13)<-c("asv","sample","clr","meas_var") # correct dimensions 
mean(z13$clr) #  will be very close to zero

meta13<- dada13 %>% 
  dplyr::select(sample, asv, phylum, abundance, id, age_class, sex) 

CLR13 <- merge(z13, meta13, by=c("asv", "sample"), all.x=TRUE) #looks good 

## clr mcmc 

MicMix16[[2]] <-MCMCglmm(clr~ 1,
                     random=~ season:asv + id:asv + asv, # note sample does not need to be included here 
                     mev=CLR$meas_var, # measurement error variance 
                     data=CLR16,
                     prior=Prior3,
                     verbose=T, 
                     pr=F, pl=F, 
                     nitt= 13000*mf, 
                     thin=10*mf, burnin=3000*mf) 

MicMix13[[2]] <-MCMCglmm(clr~ 1,
                         random=~ age_class:asv + asv, # note sample does not need to be included here 
                         mev=CLR$meas_var, # measurement error variance 
                         data=CLR13,
                         prior=Prior2,
                         verbose=T, 
                         pr=F, pl=F, 
                         nitt= 13000*mf, 
                         thin=10*mf, burnin=3000*mf) 

MCMCRep(MicMix16[[2]]) # proportion variance 
MCMCRep(MicMix13[[2]]) # propoetion variance 


# hierarchical taxonomic effect variations  ----------------------------------


Prior9 <- list(R = list(V = diag(1), nu = 0.002),
               G = list(G1 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100),
                        G2 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100),
                        G3 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100),
                        G4 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100),
                        G5 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100),
                        G6 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100),
                        G7 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100), 
                        G8 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100),
                        G9 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100)))

Prior10 <- list(R = list(V = diag(1), nu = 0.002),
               G = list(G1 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100),
                        G2 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100),
                        G3 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100),
                        G4 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100),
                        G5 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100),
                        G6 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100),
                        G7 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100), 
                        G8 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100),
                        G9 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100),
                        G10 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100))) 


MicMix16[[3]]<- MCMCglmm(fixed= mcmc13, 
                         random= ~ sample + asv + id:asv + season:asv +  
                           phylum + phylum:sample + season:phylum + 
                           fam + fam:sample + season:fam, 
                         prior=Prior10,
                         data=dada16,
                         verbose=T, pr=F, pl=F, # will save posteriors for random effect levels & model predictions for each iteration 
                         family = "poisson",
                         nitt = 13000*mf,
                         thin = 10*mf,burnin=3000*mf) 



MicMix13[[3]]<- MCMCglmm(fixed= mcmc13, 
                     random= ~ sample + asv + age_class:asv +  
                       phylum + phylum:sample + age_class:phylum + 
                       family2 + family2:sample + age_class:family2, 
                     prior=Prior9,
                     data=dada13,
                     verbose=T, pr=F, pl=F, # will save posteriors for random effect levels & model predictions for each iteration 
                     family = "poisson",
                     nitt = 13000*mf,
                     thin = 10*mf,burnin=3000*mf) 
