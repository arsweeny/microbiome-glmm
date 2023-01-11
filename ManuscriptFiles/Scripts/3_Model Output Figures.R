##### figure creation from model outputs 
rm(list=ls())

library(tidyverse); library(ggregplot); library(MCMCglmm); library(janitor)
library(pals); library(hrbrthemes); library(ggthemes)

source("ManuscriptFiles/Functions/CleanMCMC.R")
source("ManuscriptFiles/Functions/MicrobiomeFunctions.R")


MicroColours <- tol(n=12)[c(1:5, 8:12)]
MicroColoursRamp <- colorRampPalette(MicroColours)
MicroColoursExpand <- MicroColoursRamp(n=30)
MicroColoursExpand2 <- MicroColoursRamp(n=80)
MicroColoursCat <- list(
  Age = tol(12)[c(1,4)],
  Season= tol(12)[c(8,10)]
) %>% unlist()

dada16<- read.csv("Manuscript files/Data/modFiltered100DF16.csv") %>% 
  select(-c("X", "x")) 

dada13<- read_csv("Manuscript files/Data/modFiltered100DF13.csv") %>% 
  select(-c(`...1`, "x1"))


## model objects 
model13 <- MicMix13[[1]]
model16 <- MicMix16[[1]]
   
CLRmodel13 <-  MicMix13[[2]]
CLRmodel16 <-  MicMix16[[2]]

model13Taxa <- MicMix13[[3]]
model16Taxa <- MicMix16[[3]]

PropVar13<-MCMCRep(model13, scale="link") %>% mutate(year="2013")
PropVar16<-MCMCRep(model16, scale="link") %>%  mutate(year="2016")

PropVarCLR13<-MCMCRep(CLRmodel13, scale="link") %>% mutate(year="2013") %>% 
  mutate(Component=str_replace(Component, "age_class", "age"))
PropVarCLR13[3,1] <- "mev"
PropVarCLR16<-MCMCRep(CLRmodel16, scale="link") %>%  mutate(year="2016") %>% 
  mutate(Component=str_replace(Component, "asv:id", "host:asv"))
PropVarCLR16[4,1] <- "mev"
PropVarCLR16 <- PropVarCLR16 %>% mutate(Component=as.factor(Component) %>% fct_relevel(., "season:asv"))

PropVarTaxa13<-MCMCRep(model13Taxa, scale="link") %>% mutate(year="2013")
PropVarTaxa16<-MCMCRep(model16Taxa, scale="link") %>%  mutate(year="2016")


PropVarMain <- bind_rows(PropVar13, PropVar16) %>% 
  mutate(
         Component=str_replace(Component, "age_class", "age"),
         Mode = as.numeric(as.character(Mode)))

PropVarMain[10,]<-c("poisson", "0.12", NA, NA, "2013")
PropVarMain[11,]<-c("poisson", "0.10", NA, NA, "2016")

PropVarMain %>%
  mutate(Component=fct_relevel(Component,
                               "age:asv", 
                               "season:asv", "id:asv", 
                               "sample", 
                               "asv",
                               "units", "poisson") %>% 
           str_replace(., "units", "units (residual)") %>% 
           str_replace(., "id:asv", "host:asv"))-> PropVarMain


PropVarMainTaxa <- bind_rows(PropVarTaxa13, PropVarTaxa16) %>% 
  mutate(Component=str_replace(Component, "fam$", "family"),
         Component=str_replace(Component, "family2", "family"), 
         Component=str_replace(Component, "family:sample", "sample:family"),
         Component=str_replace(Component, "fam:sample", "sample:family"),
         Component=str_replace(Component, "phylum:sample", "sample:phylum"),
         Component=str_replace(Component, "age_class", "age"),
         Component=str_replace(Component, "sample_id", "sample"),
         Mode = as.numeric(as.character(Mode)))

PropVarMainTaxa[24,]<-c("poisson", "0.14", NA, NA, "2013")
PropVarMainTaxa[25,]<-c("poisson", "0.13", NA, NA, "2016")


PropVarMainTaxa %>% 
  mutate(Component=str_replace(Component, "asv|family|phylum", "taxonomy")) %>% 
  as.data.frame() %>% 
  group_by(year, Component) %>% mutate(Mode=as.numeric(Mode)) %>% 
  summarise(Mode=sum(Mode)) %>% as.data.frame() %>% 
  mutate(Component=fct_relevel(Component,
                               "age:taxonomy", 
                               #"season:taxonomy", "id:taxonomy", 
                               "sample:taxonomy", "sample", 
                               "taxonomy",
                               "units", "poisson"))-> PropVarTaxaHiLevel

### Plots 

PropVarMain %>% 
  ggplot(aes(x=year, y=as.numeric(Mode), fill=Component)) + 
  geom_bar(stat="identity", colour="transparent") + 
  scale_fill_manual(values=MicroColours) +
  theme_bw(base_size = 14) + 
  labs(x='Dataset', y='Proportion Variance') -> PropVarMainPlot 


PropVarMain13 <- PropVarMain %>% filter(year==2013) %>% droplevels() %>% 
  mutate(Component = fct_relevel(Component, "age:asv", "asv", "sample", "poisson"))

PropVarMain16 <- PropVarMain %>% filter(year==2016) %>% droplevels() %>% 
  mutate(Component = fct_relevel(Component, "season:asv", "asv", "sample", "host:asv", "poisson"))

Components13 <- PropVarMain13 %>% pull(Component) %>% as.character()
Components16 <- PropVarMain16  %>% pull(Component) %>% as.character()

Components <- union(Components13, Components16) %>% unique() 
Components <- Components[c(7,3,6,2,1,5,4)]
names(MicroColours) <- Components


PropVarMain13%>% 
  ggplot(aes(x=1, y=as.numeric(Mode), fill=Component)) + 
  geom_bar(stat="identity", colour="transparent") + 
  scale_x_continuous(labels=NULL) + 
  scale_fill_manual(values=MicroColours[Components13]) +
  theme_bw(base_size = 14) + 
  #theme(legend.position = "top") + 
  theme(axis.ticks.x = element_blank(), 
        legend.text = element_text(size=11), 
        legend.title = element_text(size=12)) + 
  labs(x='2013', y='Proportion Variance') -> PropVarMainPlot13


PropVarMain16 %>% 
  ggplot(aes(x=1, y=as.numeric(Mode), fill=Component)) + 
  geom_bar(stat="identity", colour="transparent") + 
  scale_x_continuous(labels=NULL) + 
  scale_fill_manual(values=MicroColours[Components16]) +
  theme_bw(base_size = 14) + 
  #theme(legend.position = "top") + 
  theme(axis.ticks.x = element_blank(), 
        legend.text = element_text(size=11), 
        legend.title = element_text(size=12)) + 
  labs(x='2016', y='Proportion Variance') -> PropVarMainPlot16

require(patchwork)

PropVarMainPlot13 + PropVarMainPlot16 + plot_annotation(tag_levels = "A") -> PropVarMainPlotSplit


ComponentsCLR13 <- PropVarCLR13 %>% pull(Component) %>% as.character()
ComponentsCLR16 <- PropVarCLR16  %>% pull(Component) %>% as.character()

ComponentsCLR <- union(ComponentsCLR13, ComponentsCLR16) %>% unique() 
ComponentsCLR <- ComponentsCLR[c(6,1,2,5,3,4)]

MicroColours <- tol(n=12)[c(1:5, 8:12)]
MicroCLR <- MicroColours[c(1,2,4,5,7,9)]
names(MicroCLR) <- ComponentsCLR


PropVarCLR13%>% 
  ggplot(aes(x=1, y=as.numeric(as.character(Mode)), fill=Component)) + 
  geom_bar(stat="identity", colour="transparent") + 
  scale_x_continuous(labels=NULL) + 
  scale_fill_manual(values=MicroCLR[ComponentsCLR13]) +
  theme_bw(base_size = 14) + 
  #theme(legend.position = "top") + 
  theme(axis.ticks.x = element_blank(), 
        legend.text = element_text(size=11), 
        legend.title = element_text(size=12)) + 
  labs(x='2013', y='Proportion Variance') -> PropVarCLRPlot13


PropVarCLR16 %>% 
  ggplot(aes(x=1, y=as.numeric(as.character(Mode)), fill=Component)) + 
  geom_bar(stat="identity", colour="transparent") + 
  scale_x_continuous(labels=NULL) + 
  scale_fill_manual(values=MicroCLR[ComponentsCLR16]) +
  theme_bw(base_size = 14) + 
  #theme(legend.position = "top") + 
  theme(axis.ticks.x = element_blank(), 
        legend.text = element_text(size=11), 
        legend.title = element_text(size=12)) + 
  labs(x='2016', y='Proportion Variance') -> PropVarCLRPlot16

PropVarCLRPlot13 + PropVarCLRPlot16 + plot_annotation(tag_levels = "A") -> PropVarCLRPlotSplit


PropVarTaxaHiLevel %>% 
  ggplot(aes(x=year, y=as.numeric(Mode), fill=Component)) + 
  geom_bar(stat="identity", colour="transparent") + 
  scale_fill_manual(values=MicroColours) +
  theme_bw(base_size = 14) + 
  labs(x='Dataset', y='Proportion Variance') -> PropVarMainTaxaPlot 

PropVarMainTaxa %>% 
  filter(str_detect(Component, "age:")) %>% 
  mutate(Mode=as.numeric(Mode)) %>% 
  mutate(Mode=Mode/0.20) -> PropVarTaxaAge

PropVarMainTaxa %>% 
  filter(str_detect(Component, "season:")) %>% 
  mutate(Mode=as.numeric(Mode)) %>% 
  mutate(Mode=Mode/0.01) -> PropVarTaxaSeason


PropVarTaxaAge %>% 
  ggplot(aes(x=year, y=as.numeric(Mode), fill=Component)) + 
  geom_bar(stat="identity", colour="transparent") + 
  scale_fill_brewer(palette="Purples") +
  theme_minimal(base_size = 14) + 
  theme(axis.ticks = element_blank(),
        axis.text.x = element_blank()) + 
  labs(x='', y='Proportion Variance') -> PropVarAgePlot 



PropVarMainPlot + ggsave("PropVarMainPlot.png", units="mm", height=140, width=160)
PropVarMainPlotSplit + ggsave("PropVarMainPlotSplit.png", units="mm", height=140, width=220)
PropVarCLRPlotSplit + ggsave("PropVarCLRPlotSplit.png", units="mm", height=140, width=220)
PropVarMainTaxaPlot + ggsave("PropVarMainTaxaPlot.png", units="mm", height=140, width=160)
PropVarAgePlot  + ggsave("PropVarAgePlot.png", units="mm", height=90, width=110)

write.csv(PropVarMain, "PropVarMCMC_Main.csv", row.names = FALSE)




# Model Tables -------------------------------------------------------------


# main mods 
tidyMCMC13 <- clean.MCMC.GLMM(model13) %>% 
  mutate(variable=as.factor(variable) %>%  
           recode(`age_class:asv`="age:asv", 
                 # units="units (id:asv)", 
                  `(Intercept)`= "Intercept", 
                  age_classA="ageAdult") %>% 
           as.character()) %>% 
  mutate(variable=str_replace(variable, "age_class", "age"))

tidyMCMC16 <- clean.MCMC.GLMM(model16) %>% 
  mutate(variable=as.factor(variable) %>%  
           recode(sample_id= "sample",
                  #`season:asv`="asv:season", 
                  `(Intercept)`= "Intercept", 
                  seasonsummer="seasonSummer", 
                #  units="units (season:id:asv)"
                ) 
         %>% as.character()) %>% 
  mutate(variable=str_replace(variable, "fam", "family"))

tidyMCMC13fin <- left_join(tidyMCMC13, PropVarMain %>% filter(year==2013), 
                        by=c("variable"="Component")) %>% 
  mutate_at(c("Mode", "lHPD", "uHPD"), as.character)
tidyMCMC16fin <- left_join(tidyMCMC16, PropVarMain %>% filter(year==2016),
                           by=c("variable"="Component")) %>% 
  mutate_at(c("Mode", "lHPD", "uHPD"), as.character)

# taxa mods 
tidyMCMC13taxa <- clean.MCMC.GLMM(model13Taxa) %>% 
  mutate(variable=as.factor(variable) %>%  
           recode(`age_class:asv`="age:asv", 
                #  units="units (id:asv)", 
                  `(Intercept)`= "Intercept", 
                  age_classA="ageAdult", 
                `phylum:sample`="sample:phylum", 
                `family:sample`="sample:family") %>% 
           as.character()) %>% 
  mutate(variable=str_replace(variable, "age_class", "age")) %>%
  mutate(variable=str_replace(variable, "family2", "family")) %>% 
  mutate(variable=str_replace(variable, "family:sample", "sample:family"))
  
tidyMCMC16taxa <- clean.MCMC.GLMM(model16Taxa) %>% 
  mutate(variable=as.factor(variable) %>%  
           recode(sample_id= "sample",
                  #`season:asv`="asv:season", 
                  `(Intercept)`= "Intercept", 
                  seasonsummer="seasonSummer", 
                #  units="units (season:id:asv)", 
                `phylum:sample`="sample:phylum", 
                `family:sample`="sample:family") 
         %>% as.character()) %>% 
  mutate(variable=str_replace(variable, "fam", "family")) %>% 
  mutate(variable=str_replace(variable, "family:sample", "sample:family"))

tidyMCMC13taxafin <- left_join(tidyMCMC13taxa, PropVarMainTaxa %>% filter(year==2013), 
                        by=c("variable"="Component")) %>% 
  mutate_at(c("Mode", "lHPD", "uHPD"), as.character)
tidyMCMC16taxafin <- left_join(tidyMCMC16taxa, PropVarMainTaxa%>% filter(year==2016),
                               by=c("variable"="Component")) %>% 
  mutate_at(c("Mode", "lHPD", "uHPD"), as.character)


## clr models 

tidyCLR13 <- clean.MCMC.GLMM(CLRmodel13) %>% 
  mutate(variable=as.factor(variable) %>%  
           recode(`sqrt(mev):sqrt(mev).meta`= "mev",
                  `(Intercept)`= "Intercept", 
                  `age_class:asv`="age:asv", 
                 # units="asv:id"
                 ) %>% as.character())

tidyCLR16 <- clean.MCMC.GLMM(CLRmodel16) %>% 
  mutate(variable=as.factor(variable) %>%  
           recode(`sqrt(mev):sqrt(mev).meta`= "mev",
                  `season:asv`="season:asv", 
                  `(Intercept)`= "Intercept", 
                 # units="asv:season:id"
                 ) %>% as.character())

tidyCLR13fin <- left_join(tidyCLR13, PropVarCLR13, by=c("variable"="Component")) %>% 
  mutate_at(c("Mode", "lHPD", "uHPD"), as.character)
tidyCLR16fin <- left_join(tidyCLR16, PropVarCLR16, by=c("variable"="Component")) %>% 
  mutate_at(c("Mode", "lHPD", "uHPD"), as.character)

ModelOutBase <- bind_rows(tidyMCMC13fin, tidyMCMC16fin) 
ModelOutTaxa <- bind_rows(tidyMCMC13taxafin, tidyMCMC16taxafin) 
ModelOutCLR <- bind_rows(tidyCLR13fin, tidyCLR16fin) 


write.csv(ModelOutBase, "PilotAnalysisOutputs/ModelOutputMain.csv", row.names = FALSE)
write.csv(ModelOutTaxa, "PilotAnalysisOutputs/ModelOutputTaxa.csv", row.names = FALSE)
write.csv(ModelOutCLR, "PilotAnalysisOutputs/ModelOutputCLR.csv", row.names = FALSE)


# Differential Abundances -------------------------------------------------

Solutions13<-model13$Sol[, 1:dim(model13$Z)[2]]
Solutions16<-model16$Sol[, 1:dim(model16$Z)[2]]

saveRDS(Solutions13, file="PilotAnalysisOutputs/Poisson_SolutionsMCMC2013_filtered.rds")
saveRDS(Solutions16, file="PilotAnalysisOutputs/Poisson_SolutionsMCMC2016_filtered.rds") 

Solutions13 <- readRDS(file="Poisson_SolutionsMCMC2013_filtered.rds")
Solutions16 <- readRDS(file="Poisson_SolutionsMCMC2016_filtered.rds")

Solutions13DF<- Solutions13 %>% as.data.frame()
Solutions16DF<- Solutions16 %>% as.data.frame()

SolutionsTaxa13<-model13Taxa$Sol[, 1:dim(model13Taxa$Z)[2]]
SolutionsTaxa16<-model16Taxa$Sol[, 1:dim(model16Taxa$Z)[2]]

SolutionsTaxa13DF<- SolutionsTaxa13 %>% as.data.frame()
SolutionsTaxa16DF<- SolutionsTaxa16 %>% as.data.frame()

SolKeepAge <- Solutions13DF %>% select(contains("age_class:asv")) %>% colnames()
SolKeepSeason <-Solutions16DF %>% select(contains("season:asv")) %>% colnames()

SolKeepAgeASV <- SolutionsTaxa13DF %>% select(contains("age_class:asv")) %>% colnames()
SolKeepSeasonASV <- SolutionsTaxa16DF %>% select(contains("season:asv")) %>% colnames()

SolKeepAgePhy <- SolutionsTaxa13DF %>% select(contains("age_class:phylum")) %>% colnames()
SolKeepSeasonPhy <-SolutionsTaxa16DF %>% select(contains("season:phylum")) %>% colnames()

SolKeepAgeFam <- SolutionsTaxa13DF %>% select(contains("age_class:family2")) %>% colnames()
SolKeepSeasonFam <-SolutionsTaxa16DF %>% select(contains("season:fam")) %>% colnames()

SolutionsAge<- Solutions13DF %>% select(SolKeepAge) 
SolutionsSeason<- Solutions16DF %>% select(SolKeepSeason) 

SolutionsAgeASV<- SolutionsTaxa13DF %>% select(SolKeepAgeASV) 
SolutionsSeasonASV<- SolutionsTaxa16DF %>% select(SolKeepSeasonASV) 

SolutionsAgePhy<- SolutionsTaxa13DF %>% select(SolKeepAgePhy) 
SolutionsSeasonPhy<- SolutionsTaxa16DF %>% select(SolKeepSeasonPhy) 

SolutionsAgeFam<- SolutionsTaxa13DF %>% select(SolKeepAgeFam) 
SolutionsSeasonFam<- SolutionsTaxa16DF %>% select(SolKeepSeasonFam) 


SolutionsAdult <- SolutionsAge %>% 
  gather(key="level", value="posterior") %>% as.data.frame() %>% 
  separate(level, into=c("term", "age_class", "asv"), sep="[.]")
SolutionsAdultASV <- SolutionsAgeASV %>% 
  gather(key="level", value="posterior") %>% as.data.frame() %>% 
  separate(level, into=c("term", "age_class", "asv"), sep="[.]")
SolutionsAdultPhy <- SolutionsAgePhy %>% 
  gather(key="level", value="posterior") %>% as.data.frame() %>% 
  separate(level, into=c("term", "age_class", "phylum"), sep="[.]")
SolutionsAdultFam <- SolutionsAgeFam %>% 
  gather(key="level", value="posterior") %>% as.data.frame() %>% 
  separate(level, into=c("term", "age_class", "family"), sep="[.]")


SolutionsSummer <- SolutionsSeason %>% 
  gather(key="level", value="posterior") %>% as.data.frame() %>% 
  separate(level, into=c("term", "season", "asv"), sep="[.]")
SolutionsSummerASV <- SolutionsSeasonASV %>% 
  gather(key="level", value="posterior") %>% as.data.frame() %>% 
  separate(level, into=c("term", "season", "asv"), sep="[.]")
SolutionsSummerPhy <- SolutionsSeasonPhy %>% 
  gather(key="level", value="posterior") %>% as.data.frame() %>% 
  separate(level, into=c("term", "season", "phylum"), sep="[.]")
SolutionsSummerFam <- SolutionsSeasonFam %>% 
  gather(key="level", value="posterior") %>% as.data.frame() %>% 
  separate(level, into=c("term", "season", "family"), sep="[.]")



ASVByAgeBase <- RanefPostComp(SolutionsAdult, comps=c("L", "A"), 
                              taxa="asv", var="age_class") # create a df of the comparisons 
ASVBySeasonBase<-RanefPostComp(SolutionsSummer, comps=c("spring", "summer"), 
                           taxa="asv", var="season") # create a df of the comparisons 

ASVByAge<-RanefPostComp(SolutionsAdultASV, comps=c("L", "A"), 
                        taxa="asv", var="age_class") # create a df of the comparisons 
ASVBySeason<-RanefPostComp(SolutionsSummerASV, comps=c("spring", "summer"), 
                           taxa="asv", var="season") # create a df of the comparisons 

PhylumByAge <- RanefPostComp(SolutionsAdultPhy, comps=c("L", "A"), 
                           taxa="phylum", var="age_class") # create a df of the comparisons 
PhylumBySeason <-RanefPostComp(SolutionsSummerPhy, comps=c("spring", "summer"), 
                           taxa="phylum", var="season") # create a df of the comparisons 

FamilyByAge <- RanefPostComp(SolutionsAdultFam, comps=c("L", "A"), 
                           taxa="family", var="age_class") # create a df of the comparisons 
FamilyBySeason <-RanefPostComp(SolutionsSummerFam, comps=c("spring", "summer"), 
                             taxa="family", var="season") # create a df of the comparisons 


write.csv(ASVByAgeBase, 
          "PilotAnalysisOutputs/ASVByAgeBase.csv", row.names = FALSE)
write.csv(ASVBySeasonBase, 
          "PilotAnalysisOutputs/ASVBySeasonBase.csv", row.names = FALSE)

write.csv(ASVByAge, 
          "PilotAnalysisOutputs/ASVByAge.csv", row.names = FALSE)
write.csv(ASVBySeason, 
          "PilotAnalysisOutputs/ASVBySeason.csv", row.names = FALSE)

write.csv(PhylumByAge, 
          "PilotAnalysisOutputs/PhylumByAge.csv", row.names = FALSE)
write.csv(PhylumBySeason, 
          "PilotAnalysisOutputs/PhylumBySeason.csv", row.names = FALSE)

write.csv(FamilyByAge, 
          "PilotAnalysisOutputs/FamilyByAge.csv", row.names = FALSE)
write.csv(FamilyBySeason, 
          "PilotAnalysisOutputs/FamilyBySeason.csv", row.names = FALSE)


ASVByAgeBase <- read.csv("PilotAnalysisOutputs/ASVByAgeBase.csv") 
ASVBySeasonBase <- read.csv("PilotAnalysisOutputs/ASVBySeasonBase.csv") 

taxo13 <- dada13 %>% filter(!duplicated(asv)) %>% 
  select(asv, phylum, family) %>% as.data.frame() 
  #rename(family=family2)

prev13 <- microbiome::prevalence(ps13) %>% as.data.frame() %>% 
  arrange(desc(.)) %>% 
  rename(Prevalence=1) %>% 
  mutate(Prevalence=Prevalence*100) %>% 
  rownames_to_column(var="asv")

taxo13 <- left_join(prev13, taxo13) %>% 
  left_join(., asv_prev13)

taxo16 <- dada16 %>% filter(!duplicated(asv)) %>% 
  select(asv, phylum, family) %>% as.data.frame() 

dfEndTaxAge <- merge(ASVByAgeBase, taxo13, by="asv", all.x=TRUE)
dfEndTaxSeason <- merge(ASVBySeasonBase, taxo16, by="asv", all.x=TRUE)


dfEndTaxAgeASV <- merge(ASVByAge, taxo13, by="asv", all.x=TRUE)
dfEndTaxSeasonASV <- merge(ASVBySeason, taxo16, by="asv", all.x=TRUE)


dfEndTaxAgeFam <- merge(FamilyByAge, taxo13 %>% select(family, phylum) %>% 
                        filter(!duplicated(family)), 
                      by="family", all.x=TRUE, all.y=FALSE)
dfEndTaxSeasonFam <- merge(FamilyBySeason, taxo16 %>% select(family, phylum) %>% filter(!duplicated(family)), 
                     by="family", all.x=TRUE)


dfEndTaxAge %>% group_by(phylum) %>% na.omit() %>% 
  summarise(nASV=length(asv)) %>% filter(nASV>5) %>% 
  pull(phylum) %>% as.character() -> keepPhyBase13

dfEndTaxAgeASV %>% group_by(phylum) %>% na.omit() %>% 
  summarise(nASV=length(asv)) %>% filter(nASV>5) %>% 
  pull(phylum) %>% as.character() -> keepPhy13

dfEndTaxAgeASV %>% group_by(family) %>% na.omit() %>% 
  summarise(nASV=length(asv)) %>% filter(nASV>10) %>% 
  pull(family) %>% as.character() -> keepFam13

dfEndTaxSeason %>% group_by(phylum) %>% na.omit() %>% 
  summarise(nASV=length(asv)) %>% filter(nASV>5) %>% 
  pull(phylum) %>% as.character() -> keepPhyBase16

dfEndTaxSeasonASV %>% group_by(phylum) %>% na.omit() %>% 
  summarise(nASV=length(asv)) %>% filter(nASV>5) %>% 
  pull(phylum) %>% as.character() -> keepPhy16

dfEndTaxSeasonASV %>% group_by(family) %>% na.omit() %>% 
  summarise(nASV=length(asv)) %>% filter(nASV>5) %>% 
  mutate(family=factor(family)) %>%  
  pull(family) %>% as.character()  -> keepFam16


union(keepPhyBase13, keepPhyBase16) -> keepPhyBaseAll
union(keepPhy13, keepPhy16) ->  keepPhyAll 
union(keepFam13, keepFam16) -> keepFamAll


#paste0("ASV", seq(1:500))-> keepASV

Prev13<- dada13 %>% group_by(asv) %>% 
  summarise(prevASV=length(abundance[abundance>0])/ 58) # n samples 
Prev16 <- dada16 %>% group_by(asv) %>% 
  summarise(prevASV=length(abundance[abundance>0])/ 72) # n samples 

sigAgeASV<- dfEndTaxAge %>% #filter(!is.na(phylum)) %>% 
  filter(mean>0 &lower>0| mean<0 &upper<0) 

sigSeasonASV<- dfEndTaxSeason %>% #filter(!is.na(phylum)) %>% 
  filter(mean>0 &lower>0| mean<0 &upper<0) 

sigASVPos13 <- sigAgeASV %>% 
  filter(mean>0) %>% 
  count(phylum) %>% as.data.frame() %>% 
  rename(nPos=n) %>% 
  mutate(percentPOS=nPos/sum(nPos)*100)

sigASVNeg13 <- sigAgeASV %>% 
  filter(mean<0) %>% 
  count(phylum) %>% as.data.frame() %>% 
  rename(nNeg=n) %>% 
  mutate(percentNEG=nNeg/sum(nNeg)*100)

sigASVPos16 <- sigSeasonASV %>% 
  filter(mean>0) %>% 
  count(phylum) %>% as.data.frame() %>% 
  rename(nPos=n) %>% 
  mutate(percentPOS=nPos/sum(nPos)*100)

sigASVNeg16 <- sigSeasonASV %>% 
  filter(mean<0) %>% 
  count(phylum) %>% as.data.frame() %>% 
  rename(nNeg=n) %>% 
  mutate(percentNEG=nNeg/sum(nNeg)*100)

sigAgeAll <- 
  full_join(sigASVPos13, sigASVNeg13) %>% 
  mutate(effect="Age") %>% 
  arrange(desc(nPos))

sum(sigAgeAll[, 2] %>% na.omit()) -> nPosTotal  #347 post 
sum(sigAgeAll[, 4] %>% na.omit()) -> nNegTotal  #336 

percentPos <- nPosTotal/ (nPosTotal+ nNegTotal) #50.805
percentNeg <- nNegTotal/ (nPosTotal+ nNegTotal) #49.195 


sigSeasonAll <- 
  full_join(sigASVPos16, sigASVNeg16) %>% 
  mutate(effect="Season") %>% 
  arrange(desc(nPos))

sum(sigSeasonAll[, 2] %>% na.omit()) -> nPosTotal  #16 
sum(sigSeasonAll[, 4] %>% na.omit()) -> nNegTotal # 8 

percentPos <- nPosTotal/ (nPosTotal+ nNegTotal) #33.33 
percentNeg <- nNegTotal/ (nPosTotal+ nNegTotal) #66.66 


write.csv(sigAgeAll, "PilotAnalysisOutputs/PhylaLevelSigShiftsAgeFiltered100Base.csv", row.names=FALSE)
write.csv(sigSeasonAll, "PilotAnalysisOutputs/PhylaLevelSigShiftsSeasonFiltered100Base.csv", row.names=FALSE)


top10 <- dfEndTaxSeason %>% filter(phylum %in% keepPhyBaseAll) %>% 
  arrange(-mean) %>% select(asv, phylum, lower, upper, mean) %>% 
  filter(lower>0) %>%  #make sure errors don't span zero 
  top_n(10) %>% 
  mutate(group="real")

bottom10 <- dfEndTaxSeason %>% filter(phylum %in% keepPhyBaseAll) %>% 
  arrange(-mean) %>% select(asv, phylum, lower, upper, mean) %>%  
  filter(upper<0) %>% 
  top_n(-10) %>% 
  mutate(group="real")

top100 <- dfEndTaxAge %>% arrange(-mean) %>% select(asv, phylum, lower, upper, mean) %>% 
  #filter(lower>0) %>%  #make sure errors don't span zero 
  top_n(100) %>% 
  mutate(group="real")

bottom100 <- dfEndTaxAge %>% arrange(-mean) %>% select(asv, phylum, lower, upper, mean) %>%  
  #filter(upper<0.1) %>% 
  top_n(-100) %>% 
  mutate(group="real")

top50 <- dfEndTaxAge %>% arrange(-mean) %>% filter(phylum %in% keepPhyBaseAll) %>% 
  select(asv, phylum, lower, upper, mean) %>% 
  #filter(lower>0) %>%  #make sure errors don't span zero 
  top_n(50) %>% 
  mutate(group="real")

bottom50 <- dfEndTaxAge %>% arrange(-mean) %>% filter(phylum %in% keepPhyBaseAll) %>% 
  select(asv, phylum, lower, upper, mean) %>%  
  #filter(upper<0.1) %>% 
  top_n(-50) %>% 
  mutate(group="real")


topAgePos <- top50 %>% 
  count(phylum) %>% as.data.frame() %>% 
  rename(nPos=n) 
topAgeNeg <- bottom50 %>% 
  count(phylum) %>% as.data.frame() %>% 
  rename(nNeg=n) 


dummiesAge<-paste0("dummy", rep(1:20))
groupAge<-rep("dummy", 20)

dummiesSzn <- paste0("dummy", rep(1:5))
groupSzn<-rep("dummy", 5)

dummyAge<- cbind(dummiesAge, groupAge) %>% as.data.frame()
colnames(dummyAge)<- c("asv", "group")

dummySzn<- cbind(dummiesSzn, groupSzn) %>% as.data.frame()
colnames(dummySzn)<- c("asv", "group")

select10<- bind_rows(top10, dummySzn, bottom10) 
select50<- bind_rows(top50, dummyAge, bottom50) 

select_x_age<- select50 %>% pull(asv)
select_x_szn<- select10 %>% pull(asv)

MicroColoursPhy <- MicroColoursExpand[c(TRUE, FALSE)][3:14]

PlotPal <- data.frame(keepPhyBaseAll, MicroColoursPhy)
colnames(PlotPal) <- c("phylum", "colour")


SeasonTopBottom<-
  select10 %>% filter(!is.na(phylum)) %>% 
  ggplot(aes(x=reorder(asv, mean), y=mean, colour=phylum))+
  geom_point(position=position_dodge(w=0.4), size=1.5)+
  geom_errorbar(position=position_dodge(w=0.5),
                aes(ymin=lower,ymax=upper),size=0.3,width=0.2)+
  geom_segment(aes(y=0, yend=0, x=25, xend=15), colour="black", size=0.7)+
  geom_segment(aes(y=0, yend=0, x=15, xend=10), linetype="dashed", colour="black", size=0.7)+
  geom_segment(aes(y=0, yend=0, x=10, xend=0), colour="black", size=0.7)+
  #geom_hline(aes(yintercept=0, linetype=group))+
  theme_classic(base_size = 14)+labs(x=NULL)+
  scale_colour_manual(values=PlotPal %>% 
                        filter(phylum %in% select10$phylum) %>% 
                        mutate(colour=as.character(colour)) %>% 
                        pull(colour), 
                      breaks=unique(as.character(select10$phylum)) %>% 
                        sort()) +
  scale_x_discrete(limits=rev(select_x_szn))+
  #geom_vline(aes(xintercept=50.5, colour="grey"), size=0.2) + 
  coord_flip()+
  theme(strip.text.x = element_text(size = 0.8), 
        strip.background =element_rect(fill="white"),
        axis.text.y = element_blank(), 
        axis.ticks.y=element_blank())+
  labs(subtitle = "D")
#labs(title="ASV-by-Season Estimate, Spring-Summer Shift", subtitle = "Top & Bottom 10") 

SeasonTopBottom+ggsave('SeasonTopBottom.tiff', units="mm", height=160, width=130)
SeasonTopBottom+ggsave('SeasonTopBottom.png', units="mm", height=160, width=130)



AgeTopBottom50<-
  select50 %>%filter(!is.na(phylum)) %>% 
  ggplot(aes(x=reorder(asv, mean), y=mean, colour=phylum))+
  geom_point(position=position_dodge(w=0.4), size=1.5)+
  geom_errorbar(position=position_dodge(w=0.5),
                aes(ymin=lower,ymax=upper),size=0.3,width=0.2)+
  geom_segment(aes(y=0, yend=0, x=120, xend=70), colour="black", size=0.7)+
  geom_segment(aes(y=0, yend=0, x=70, xend=50), linetype="dashed", colour="black", size=0.7)+
  geom_segment(aes(y=0, yend=0, x=50, xend=0), colour="black", size=0.7)+
  #geom_hline(aes(yintercept=0, linetype=group))+
  theme_classic(base_size = 14)+labs(x=NULL)+
  scale_colour_manual(values=PlotPal %>% 
                        filter(phylum %in% select50$phylum) %>% 
                        mutate(colour=as.character(colour)) %>% 
                        pull(colour), 
                      breaks=unique(as.character(select50$phylum)) %>% 
                        sort()) + 
  scale_x_discrete(limits=rev(select_x_age))+
  #geom_vline(aes(xintercept=50.5, colour="grey"), size=0.2) + 
  coord_flip()+
  theme(strip.text.x = element_text(size = 0.8), 
        strip.background =element_rect(fill="white"),
        axis.text.y = element_blank(), 
        axis.ticks.y=element_blank())+
  labs(subtitle = "B")
#labs(title="ASV-by-Age Estimate, Lamb-Adult Shift", subtitle = "Top & Bottom 50") 

AgeTopBottom100<-
  select100 %>%filter(!is.na(phylum)) %>% 
  ggplot(aes(x=reorder(asv, mean), y=mean, colour=phylum))+
  geom_point(position=position_dodge(w=0.4), size=1.5)+
  geom_errorbar(position=position_dodge(w=0.5),
                aes(ymin=lower,ymax=upper),size=0.3,width=0.2)+
  geom_segment(aes(y=0, yend=0, x=240, xend=140), colour="black", size=0.7)+
  geom_segment(aes(y=0, yend=0, x=140, xend=100), linetype="dashed", colour="black", size=0.7)+
  geom_segment(aes(y=0, yend=0, x=100, xend=0), colour="black", size=0.7)+
  #geom_hline(aes(yintercept=0, linetype=group))+
  theme_classic(base_size = 14)+labs(x=NULL)+
  scale_colour_viridis_d(option="plasma")+
  scale_x_discrete(limits=rev(select_x))+
  #geom_vline(aes(xintercept=50.5, colour="grey"), size=0.2) + 
  coord_flip()+
  theme(strip.text.x = element_text(size = 0.8), 
        strip.background =element_rect(fill="white"),
        axis.text.y = element_blank(), 
        axis.ticks.y=element_blank())+
  labs(title="ASV-by-Age Estimate, Lamb-Adult Shift", subtitle = "Top & Bottom 100") 


AgeTopBottom50+ggsave('AgeTopBottom.tiff', units="mm", height=180, width=160)
AgeTopBottom100+ggsave('AgeTopBottom100.tiff', units="mm", height=180, width=160)
AgeTopBottom+ggsave('AgeTopBottom.tiff', units="mm", height=180, width=160)
AgeTopBottom50+ggsave('AgeTopBottom.png', units="mm", height=180, width=160)



ForestAge<-list()

Phyla<-dada13 %>% filter(phylum %in% keepPhy) %>% 
  group_by(phylum) %>% 
  summarise(nPhy=length(abundance)) %>% 
  arrange(desc(nPhy)) %>% 
  pull(phylum) 

BiCols <-c("#1874CD", "#FF1493")
Viridis<- viridis_pal(alpha=1, option="D")

for(x in 1:length(Phyla)){
  require(ggExtra)
  p<-
    dfEndTax2 %>% filter(phylum == Phyla[x]) %>% 
    ggplot(aes(x=reorder(asv, mean), y=mean, colour=phylum))+
    geom_point(position=position_dodge(w=0.5), size=0.3)+
    geom_errorbar(position=position_dodge(w=0.5),
                  aes(ymin=lower,ymax=upper),size=0.3,width=0.2)+
    scale_color_gradient(high=BiCols[2], low=BiCols[1], limits=c(0,1),
                         breaks=c(seq(0,1, by=0.2)), name="ASV Prevalence")+
    geom_hline(aes(yintercept=0),lty=2, size=0.9)+
    theme_classic(base_size = 14)+labs(x=NULL)+coord_flip()+
    labs(title=Phyla[x], y="estimate")+
    theme(strip.text.x = element_text(size = 12),
          axis.text.y = element_blank(), 
          axis.title.x = element_text(size=10),
          plot.title = element_text(hjust=0.5, size=14, face="bold"),
          legend.position = "none",
          legend.key.size = unit(0.5, "cm"), 
          legend.text = element_text(size=8, angle = 45, hjust=1), 
          plot.margin = unit(c(2,2,2,2), "lines"))
  
  ForestAge[[x]]<- ggMarginal(p=p, type="density", margins = "x", size=4) 
}

require(ggpubr)
require(patchwork)
require(cowplot)

# ggarrange option 
AgePhyla<-ggarrange(ForestAge[[1]], ForestAge[[2]], ForestAge[[3]], ForestAge[[4]], 
                    ForestAge[[5]], ForestAge[[6]], ForestAge[[7]], ForestAge[[8]], 
                    ForestAge[[9]], ForestAge[[10]], ForestAge[[11]], ForestAge[[12]], ncol=4, nrow=3)

legend<- get_legend(ForestAge[[1]])  #+theme(legend.position="bottom"))

title_gg <- ggplot() + 
  labs(title = "ASV-by-Age Effects", subtitle = "Lamb-Adult Shift") +
  theme(plot.title=element_text(size=18),
        plot.subtitle=element_text(size=14, face="italic"))

# patchwork option 
AgePhylaFin <- plot_grid(title_gg, AgePhyla, legend, ncol=1, rel_heights = c(0.1, 1, 0.1)) + 
  ggsave("ForestAgeASVFin.tiff", units="mm", height=250, width=280, dpi=250)


dfEndTax  %>% na.omit() %>% subset(phylum %in% keepPhy) %>% 
  ggplot(aes(y=mean, x=effect, size=1/var, colour=phylum)) +
  geom_jitter(alpha=0.35)+
  scale_x_discrete(limits=c("season"), labels=c("season:asv")) +
  scale_size(range=c(0.06, 7), name="Inverse Variance") +
  scale_colour_viridis_d(guide=FALSE, option="d")+ geom_hline(yintercept=0, linetype="dashed")+
  coord_flip() +theme_ipsum(base_size = 10) + theme(axis.text.x = element_text())+facet_wrap(~phylum) # +
ggsave("FruitCloudCLR.tiff", units="mm", height=230, width=340, dpi=150)

library(ggExtra)
library(hrbrthemes)
library(ggpointdensity)

Cloud13<-
  dfEndTax2  %>% na.omit() %>% subset(phylum %in% keepPhy) %>% 
  ggplot(aes(y=mean, x=prevASV, size=1/var, colour=phylum)) +
  #geom_pointdensity() +  
  geom_jitter(alpha=0.35)+
  #scale_x_discrete(limits=c("age"), labels=c("age:asv")) +
  scale_size(range=c(0.3, 5), name="Inverse Variance") +
  scale_colour_viridis_d()+ 
  geom_hline(yintercept=0, linetype="dashed")+
  coord_flip() +theme_ipsum(base_size = 10) + 
  theme(axis.text.x = element_text()) + facet_wrap(~phylum) +
  ggsave("FruitCloudCLR_2013.tiff", units="mm", height=230, width=340, dpi=200)


dfEndTaxAge  %>% na.omit() %>% subset(phylum %in% keepPhy13) -> SinaDFAgeBase
dfEndTaxAgeASV  %>% na.omit() %>% subset(phylum %in% keepPhy13) -> SinaDFAgeASV
dfEndTaxAgeFam  %>% na.omit() %>% subset(family %in% keepFam13) %>% 
  mutate(phylum=factor(phylum)) -> SinaDFAgeFam
PhylumByAge %>% na.omit() -> SinaDFAgePhy  #%>% subset(phylum %in% keepPhy13) 


dfEndTaxSeason  %>% na.omit() %>% subset(phylum %in% keepPhy16) -> SinaDFSeasonBase
dfEndTaxSeasonASV  %>% na.omit() %>% subset(phylum %in% keepPhy16) -> SinaDFSeason
dfEndTaxSeasonFam  %>% na.omit() %>% subset(family %in% keepFam16) %>% 
  mutate(phylum=factor(phylum), family=factor(family)) -> SinaDFSeasonFam 
PhylumBySeason %>% na.omit() -> SinaDFSeasonPhy  #%>% subset(phylum %in% keepPhy13) 

Sina13<-
  ggplot() +
  geom_violin(data=SinaDFAgeBase,aes(y=mean, x=phylum, fill=phylum), colour="transparent", alpha=0.2) + 
  geom_sina(data=SinaDFAgeBase,aes(y=mean, x=phylum, size=1/var, colour=phylum), alpha=0.3) + 
  scale_size(range=c(0.01, 3), name="Inverse Variance") +
  scale_colour_manual(values=PlotPal %>% 
                        filter(phylum %in% SinaDFAgeBase$phylum) %>% 
                        mutate(colour=as.character(colour)) %>% 
                        pull(colour), 
                      breaks=unique(as.character(SinaDFAgeBase$phylum)) %>% 
                        sort()) + 
  scale_fill_manual(values=PlotPal %>% 
                        filter(phylum %in% SinaDFAgeBase$phylum) %>% 
                        mutate(colour=as.character(colour)) %>% 
                        pull(colour), 
                      breaks=unique(as.character(SinaDFAgeBase$phylum)) %>% 
                        sort()) + 
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw(base_size = 14) + 
  #coord_flip() +theme_ipsum(base_size = 10) + 
  theme(axis.text.x = element_text(angle=45, hjust=1), 
        axis.title.x = element_blank(),
        legend.position = "top")+ #+facet_wrap(~phylum) +
  guides(colour=FALSE, fill=FALSE)

Sina13Rect<-Sina13+
  geom_rect(data=top50, 
            aes(xmin=0.75, xmax=12.25, ymin=min(mean-0.10), ymax=max(mean+0.10)),
            fill="transparent", color="grey60",size=0.5) +
  geom_rect(data=bottom50,
            aes(xmin=0.75, xmax=12.25, ymin=min(mean-0.30), ymax=max(mean+0.25)),
            fill="transparent", color="grey60",size=0.5)+
  labs(subtitle = "A")


Sina13Rect+ggsave("SinaCLR_2013_Rect.tiff", units="mm", height=150, width=180, dpi=400)
Sina13Rect+ggsave("SinaCLR_2013_Rect.png", units="mm", height=150, width=180, dpi=400)
Sina13+ggsave("SinaCLR_2013.tiff", units="mm", height=150, width=180, dpi=400)


fam_ordered <- SinaDFAgeFam %>% arrange(phylum) %>% pull(family) %>% unique

Sina13Fam<-
  ggplot() +
  geom_sina(data=SinaDFAgeASV %>% filter(family %in% keepFam13), 
              aes(y=mean, x=family), size=2, alpha=0.08) + 
  geom_point(data=SinaDFAgeFam,aes(y=mean, x=family, colour=phylum), 
             size=3) + 
  geom_errorbar(data=SinaDFAgeFam,
                aes(x=family, ymax=upper, ymin=lower, colour=phylum), width=0.4)+
  scale_colour_manual(values=PlotPal %>% 
                        filter(phylum %in% SinaDFAgeASV$phylum) %>% 
                        mutate(colour=as.character(colour)) %>% 
                        pull(colour), 
                      breaks=unique(as.character(SinaDFAgeASV$phylum)) %>% 
                        sort()) + 
  #scale_fill_viridis_d(option="plasma")+ 
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw(base_size = 14) + 
  scale_x_discrete(limits=fam_ordered) + 
  #coord_flip() +theme_ipsum(base_size = 10) + 
  theme(axis.text.x = element_text(angle=45, hjust=1), 
        axis.title.x = element_blank(),
        legend.position = "none")+ #+facet_wrap(~phylum) +
  guides(fill=FALSE)

Sina16Fam<-
  ggplot() +
  geom_sina(data=SinaDFSeason[SinaDFSeason$family %in% SinaDFSeasonFam$family,] %>% 
              mutate(family=factor(family)),
            aes(y=mean, x=family), size=2, alpha=0.08) + 
  geom_point(data=SinaDFSeasonFam,aes(y=mean, x=family, colour=phylum), 
             size=3) + 
  geom_errorbar(data=SinaDFSeasonFam,
                aes(x=family, ymax=upper, ymin=lower, colour=phylum), width=0.4)+
  #scale_size(range=c(0.01, 3), name="Inverse Variance") +
  scale_colour_manual(values=MicroColours) + 
  #scale_fill_viridis_d(option="plasma")+ 
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw(base_size = 14) + 
  scale_x_discrete(limits=fam_ordered) + 
  #coord_flip() +theme_ipsum(base_size = 10) + 
  theme(axis.text.x = element_text(angle=45, hjust=1), 
        axis.title.x = element_blank(),
        legend.position = "right")+ #+facet_wrap(~phylum) +
  guides(fill=FALSE)


require(patchwork)
  
(Sina13Fam / Sina16Fam) + plot_annotation(tag_levels = "A") -> SinaNew 

Sina13Fam + ggsave("SinaPlotTaxa13.png", units="mm", height=150, width=300, dpi=500)

SinaNew + ggsave("SinaPlotTaxa.png", units="mm", height=300, width=300, dpi=600)

Sina13Phy<-
  ggplot() +
  geom_point(data=SinaDFAgePhy,aes(y=mean, x=phylum, colour=phylum), 
             size=3) + 
  geom_errorbar(data=SinaDFAgePhy,
                aes(x=phylum, ymax=upper, ymin=lower, colour=phylum), width=0.4) + 
  #geom_pointdensity() +  
  #geom_jitter(alpha=0.2)+
  #scale_x_discrete(limits=c("age"), labels=c("age:asv")) +
  scale_size(range=c(0.01, 3), name="Inverse Variance") +
  scale_colour_viridis_d(option="plasma")+ 
  scale_fill_viridis_d(option="plasma")+ 
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw(base_size = 14) + 
  #coord_flip() +theme_ipsum(base_size = 10) + 
  theme(axis.text.x = element_text(angle=45, hjust=1), 
        axis.title.x = element_blank(),
        legend.position = "top")+ #+facet_wrap(~phylum) +
  guides(colour=FALSE, fill=FALSE)



dfEndTaxSeason%>% na.omit() %>% subset(phylum %in% keepPhy16) ->SinaDF16 
Sina16<-
  ggplot() +
  geom_violin(data=SinaDF16,aes(y=mean, x=phylum, fill=phylum), colour="transparent", alpha=0.2) + 
  geom_sina(data=SinaDF16,aes(y=mean, x=phylum, size=1/(var*2), colour=phylum), alpha=0.3) + 
  #geom_pointdensity() +  
  #geom_jitter(alpha=0.2)+
  #scale_x_discrete(limits=c("age"), labels=c("age:asv")) +
  scale_size(range=c(0.01, 3), name="Inverse Variance") +
  scale_colour_manual(values=PlotPal %>% 
                        filter(phylum %in% SinaDFSeason$phylum) %>% 
                        mutate(colour=as.character(colour)) %>% 
                        pull(colour), 
                      breaks=unique(as.character(SinaDFSeasonBase$phylum)) %>% 
                        sort()) + 
  scale_fill_manual(values=PlotPal %>% 
                      filter(phylum %in% SinaDFSeasonBase$phylum) %>% 
                      mutate(colour=as.character(colour)) %>% 
                      pull(colour), 
                    breaks=unique(as.character(SinaDFSeasonBase$phylum)) %>% 
                      sort()) + 
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw(base_size = 14) + 
  #coord_flip() +theme_ipsum(base_size = 10) + 
  theme(axis.text.x = element_text(angle=45, hjust=1), 
        axis.title.x = element_blank(),
        legend.position = "top")+ #+facet_wrap(~phylum) +
  guides(colour=FALSE, fill=FALSE)

Sina16Rect<-Sina16+
  geom_rect(data=top10, 
            aes(xmin=0.75, xmax=10.25, ymin=min(mean-0.05), ymax=max(mean+0.05)),
            fill="transparent", color="grey60",size=0.5) +
  geom_rect(data=bottom10,
            aes(xmin=0.75, xmax=10.25, ymin=min(mean-0.05), ymax=max(mean+0.05)),
            fill="transparent", color="grey60",size=0.5)+
  labs(subtitle = "C")


Sina16Rect+ggsave("SinaCLR_2016_Rect.tiff", units="mm", height=150, width=180, dpi=400)
Sina16Rect+ggsave("SinaCLR_2016_Rect.png", units="mm", height=150, width=180, dpi=400)
Sina16+ggsave("SinaCLR_2016.tiff", units="mm", height=150, width=180, dpi=400)


require(patchwork)

TopBottoms<- Sina13Rect + AgeTopBottom50 + Sina16Rect + SeasonTopBottom  +
  plot_layout(widths = c(2.2, 1))

TopBottoms13<- Sina13Rect + AgeTopBottom50 +
  plot_layout(widths = c(2.2, 1))

TopBottoms16<- Sina16Rect + SeasonTopBottom  +
  plot_layout(widths = c(2.2, 1))



TopBottoms + ggsave('TopsNBottoms.tiff', units="mm", height=250, width=250, dpi=500)
TopBottoms + ggsave('TopsNBottoms.png', units="mm", height=250, width=250, dpi=500)

TopBottoms13 + ggsave('TopsNBottoms13.tiff', units="mm", height=160, width=250, dpi=500)
TopBottoms16 + ggsave('TopsNBottoms16.tiff', units="mm", height=160, width=250, dpi=500)

