

# admin -------------------------------------------------------------------

## packages 
library(microbiome)
library(knitr)
library(ggpubr)
library(reshape2)
library(RColorBrewer)
library(microbiomeutilities)
library(viridis)
library(tibble)
library(metagMisc)
library(tidyverse)
library(janitor)
library(phyloseq)
library(pals)


MicroColours <- tol(n=12)[c(1:5, 8:10)]
MicroColoursRamp <- colorRampPalette(MicroColours)
MicroColoursExpand <- MicroColoursRamp(n=20)
MicroColoursCat <- list(
  Age = tol(12)[c(1,4)],
  Season= tol(12)[c(8,10)]
) %>% unlist()


# data  -------------------------------------------------------------------
ps16<-readRDS("Manuscript files/Data/phyloseq_2016_models.rds")
ps13<-readRDS("Manuscript files/Data/phyloseq_2013_models.rds")
ps16_2 <- prune_taxa(taxa_sums(ps16) >100, ps16) # total abundance filter 100 
ps13_2 <- prune_taxa(taxa_sums(ps13) >100, ps13) # total abundance filter 100


dada16<- read.csv("Manuscript files/Data/output_dada_2016.csv") %>% clean_names() %>% 
  mutate(season=case_when(sample_month==4 ~ "spring",
                          sample_month==8 ~ "summer"),
         asv_id=paste(asv, id, sep="_")) %>% 
  filter(!grepl("star", sample)) %>%  #take out the repeats  
  mutate(sample=factor(sample))  #618696 obs 

dada13<- read_csv("Manuscript files/Data/output_dada_2013.csv") %>% clean_names()  %>% 
  mutate(asv_id = paste(asv, id, sep="_")) %>% 
  filter(!is.na(age_class), !sample_id %in% c("S2", "S4")) %>% # take out duplicated samples 
  mutate(age_class=factor(age_class, levels=c("L", "A"))) # 299048 obs 


samdf16<-read_csv("Manuscript files/Data/sample_meta_2016.csv") %>% as.data.frame() %>% 
  mutate(sample_month=as.factor(as.character(sample_month))) %>% 
  filter(!grepl("star", sample_id)) %>%  #take out the repeats  
  mutate(sample_id=factor(sample_id), 
         Season=case_when(sample_month==4 ~ "Spring",
                          sample_month==8 ~ "Summer")) #72 levels 

samdf13<-read_csv("Manuscript files/Data/sample_meta_2013.csv") %>% clean_names() %>% 
  as.data.frame() %>% 
  filter(!is.na(age_class), !sample_id %in% c("S2", "S4")) %>% # take out duplicated samples
  mutate(sample_id=factor(sample_id)) %>% #58 levels 
  mutate(age_class=factor(age_class, levels=c("L", "A")),
         Age=case_when(age_class=="L"~ "Lamb",
                       age_class=="A"~"Adult"), 
         Age=factor(Age, levels=c("Lamb", "Adult"))) 



# functions ---------------------------------------------------------------

source("Manuscript files/Functions/CleanMCMC.R")
source("Manuscript files/Functions/MicrobiomeFunctions.R")
