library(MCMCglmm);library(tidyverse); library(magrittr)

# model outputs -----------------------------------------------------------

model13 <- MicMix13[[1]]
model16 <- MicMix16[[1]]

data13 <- read.csv("PilotAnalysisDataObjects/modFiltered100DF13.csv") %>% 
  mutate(age_class=as.factor(age_class) %>% fct_relevel(., "L", "A")) %>% 
  select(-family2)
data16 <- read.csv("PilotAnalysisDataObjects/modFiltered100DF16.csv")

# extract solutions -------------------------------------------------------

Solutions13<-model13$Sol[, 1:dim(model13$Z)[2]]
Solutions16<-model16$Sol[, 1:dim(model16$Z)[2]]

Solutions13DF<- Solutions13 %>% as.data.frame()
Solutions16DF<- Solutions16 %>% as.data.frame()

# new models saving latent variables --------------------------------------

# set-up ------------------------------------------------------------------

resp1<-"abundance"
fixed16<- c("season")
fixed13<- c("age_class")

mcmc16<- as.formula(paste(resp1, "~", paste(fixed16)))
mcmc13<- as.formula(paste(resp1, "~", paste(fixed13)))

mf <- 20 

Prior3 <- list(R = list(V = diag(1), nu = 0.002),
               G = list(G1 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100),
                        G2 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100),
                        G3 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100)))


Prior4 <- list(R = list(V = diag(1), nu = 0.002),
               G = list(G1 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100),
                        G2 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100),
                        G3 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100),
                        G4 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100)))


# 2016 season  ------------------------------------------------------------

ptm <- proc.time()
Latent16 <- MCMCglmm(fixed= mcmc16, 
                     random= ~  sample_id + asv + id:asv + season:asv, 
                     prior=Prior4,
                     data=data16,
                     family = "poisson",
                     pr=TRUE,
                     pl=TRUE, 
                     nitt = 13000*mf,#REMEMBER YOU'VE DONE THIS
                     thin = 10*mf,burnin=3000*mf)
proc.time() -ptm  # it took 100 hours lol oop 

saveRDS(Latent16, "/Volumes/EcoWithin1/ModelOutputs/PilotAnalysis/LatentModel16.rds")


LatentM16 <- Latent16$Liab %>% as.matrix 
ActualZero16 <- sum(data16$abundance==0) 

saveRDS(LatentM16, "/Volumes/EcoWithin1/ModelOutputs/PilotAnalysis/LatentMatrix16.rds")
LatentM16 <- readRDS("/Volumes/EcoWithin1/ModelOutputs/PilotAnalysis/LatentMatrix16.rds")

lapply(1:dim(LatentM16)[1], 
       function(i){
         print(paste0("Iteration", i))
         LatentM16[i,] %>% map(function(x) exp(-(exp(x)))) %>% 
           unlist %>% sum 
       }) -> PredictedZeroes16  

PredictedZeroes16 %>% unlist %>% 
  map(function(x) x-ActualZero16) %>% unlist -> ExcessZeroes16

PredictedZeroes16 %>% unlist %>% 
  map(function(x) ActualZero16/x) %>% unlist -> PropExcess16 


ggplot(ExcessZeroes16 %>% as.data.frame, aes(ExcessZeroes16)) + 
  geom_histogram() + 
  theme_bw(base_size = 16) + 
  ggtitle('Predicted-Actual Zeroes 2016\n(Actual Zeroes=118,622)') -> Excess16
# there are p reliably 8000 more actual zeroes and excess distribution doesn't span zero which isn't great
# but the excess are less than 10% (6.7%) of the actual 

ggplot(PropExcess16 %>% as.data.frame, aes(PropExcess16)) + 
  geom_histogram() + 
  theme_bw(base_size = 16) + 
  ggtitle('Actual:Predicted Zeroes 2016') +
  labs(x='Proportion Excess') -> ExcessProp16

lapply(1:dim(LatentM16)[1], 
       function(i){
         print(paste0("Iteration", i))
         LatentM16[i,] %>% map(function(x) exp(x)) %>% unlist %>% mean 
       }) -> PredictedMeans16 

PredictedMeans16 %>% unlist %>% mean # 22.33043 
data16$abundance %>% mean # 22.33029 
#these are stupidly close so model isn't being too weird 


# 2013 age  ---------------------------------------------------------------

ptm <- proc.time()
Latent13 <- MCMCglmm(fixed= mcmc13, 
                     random= ~  sample_id + asv + age_class:asv, 
                     prior=Prior3,
                     data=data13, 
                     family = "poisson",
                     pr=TRUE,
                     pl=TRUE, 
                     nitt = 13000*mf,#REMEMBER YOU'VE DONE THIS
                     thin = 10*mf,burnin=3000*mf)
proc.time() -ptm  # 16 hours 

saveRDS(Latent13, "/Volumes/EcoWithin1/ModelOutputs/PilotAnalysis/LatentModel13.rds")


LatentM13 <- Latent13$Liab %>% as.matrix 
ActualZero13 <- sum(data13$abundance==0) 

saveRDS(LatentM13, "/Volumes/EcoWithin1/ModelOutputs/PilotAnalysis/LatentMatrix13.rds") 
LatentM13 <- readRDS("/Volumes/EcoWithin1/ModelOutputs/PilotAnalysis/LatentMatrix13.rds")


lapply(1:dim(LatentM13)[1], 
       function(i){
         print(paste0("Iteration", i))
         LatentM13[i,] %>% map(function(x) exp(-(exp(x)))) %>% 
           unlist %>% sum 
       }) -> PredictedZeroes13  

PredictedZeroes13 %>% unlist %>% 
map(function(x) x-ActualZero13) %>% unlist -> ExcessZeroes13

PredictedZeroes13 %>% unlist %>% 
  map(function(x) ActualZero13/x) %>% unlist -> PropExcess13 


ggplot(ExcessZeroes13 %>% as.data.frame, aes(ExcessZeroes13)) + 
  geom_histogram() + 
  theme_bw(base_size = 16) + 
  ggtitle('Predicted-Actual Zeroes 2013\n(Actual Zeroes=90,972)') -> Excess13
# there are p reliably 4000 more actual zeroes and excess distribution doesn't span zero which isn't great
# but the excess are less than 10% (4.7%) of the actual 

ggplot(PropExcess13 %>% as.data.frame, aes(PropExcess13)) + 
  geom_histogram() + 
  theme_bw(base_size = 16) + 
  ggtitle('Actual:Predicted Zeroes 2013') + 
  labs(x='Proportion Excess') -> ExcessProp13


lapply(1:dim(LatentM13)[1], 
       function(i){
         print(paste0("Iteration", i))
         LatentM13[i,] %>% map(function(x) exp(x)) %>% unlist %>% mean 
       }) -> PredictedMeans13 

PredictedMeans13 %>% unlist %>% mean # 20.07958 
data13$abundance %>% mean # 20.07891  
#these are stupidly close so model isn't being too weird 

require(patchwork)

Excess13 + Excess16

ExcessProp13 + ExcessProp16 
