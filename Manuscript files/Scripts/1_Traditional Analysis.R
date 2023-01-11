### Code for alpha & beta diversity analysis 


# set-up  ----------------------------------------------------------------

source('Manuscript files/Scripts/0_Set-up.R')

# summary stats  ----------------------------------------------------------

## sample numbers by groups 
samdf16 %>% count(Season, sex) -> counts16 
samdf16 %>% pull(age) %>% summary() 

samdf13 %>% count(Age, sex) -> counts13 
samdf13 %>% filter(Age=="Adult") %>%  pull(age) %>% summary()

## read summaries 
SeqDepth13 = rowSums(otu_table(ps13))
sample_data(ps13_2)$SeqDepth = SeqDepth13

SeqDepth16 = rowSums(otu_table(ps16))
sample_data(ps16_2)$SeqDepth = SeqDepth16

# community summary -------------------------------------------------------

## Aggregate the reaeds by phlum 
psdat.phy13 <- tax_glom(ps13_2, taxrank = "phylum") 
ps.melt.phy13 <- psmelt(psdat.phy13) 

psdat.phy16 <- tax_glom(ps16_2, taxrank = "phylum") 
ps.melt.phy16 <- psmelt(psdat.phy16) 

## create plot objects for relative abundance 
Rel13<-ggplot(ps.melt.phy13, aes(x = Sample, y = Abundance, fill = phylum)) + 
  geom_bar(stat = "identity", position = "fill", colour="transparent") + 
  scale_fill_manual(values=MicroColoursExpand) + 
  facet_grid(~Age, scales="free_x") + 
  theme_bw(base_size = 14) + # labs(subtitle = "A") + 
  guides(fill=guide_legend(ncol=2)) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=6), 
        strip.background = element_rect(fill="white"))

Rel16<-ggplot(ps.melt.phy16, aes(x = Sample, y = Abundance, fill = phylum)) + 
  geom_bar(stat = "identity", position = "fill", colour="transparent") + 
  scale_fill_manual(values=MicroColoursExpand) + 
  facet_grid(~Season, scales="free_x") + 
  theme_bw(base_size = 14) + # labs(subtitle = "B") + 
  guides(fill=guide_legend(ncol=2)) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=6),
        strip.background = element_rect(fill="white"))

## create summaries for text 

SampleDepth16 <- ps.melt.phy16 %>% 
  group_by(Sample) %>% 
  summarise(Depth=sum(Abundance)) %>% as.data.frame() 

SampleDepth13 <- ps.melt.phy13 %>% 
  group_by(Sample) %>% 
  summarise(Depth=sum(Abundance)) %>% as.data.frame() 

RelAbundance13 <-full_join(ps.melt.phy13, SampleDepth13) %>% 
  group_by(Sample, phylum) %>% 
  summarise(RelAbundance=Abundance/Depth) %>% as.data.frame()

RelAbundance16 <-full_join(ps.melt.phy16, SampleDepth16) %>% 
  group_by(Sample, phylum) %>% 
  summarise(RelAbundance=Abundance/Depth) %>% as.data.frame()

require(plotrix)
RelAbundanceMean13 <- RelAbundance13 %>% 
  select(-Sample) %>%  group_by(phylum) %>% 
  summarise_each(funs(mean, sd, std.error)) %>% 
  arrange(desc(mean))
RelAbundanceMean16 <- RelAbundance16 %>% 
  select(-Sample) %>%  group_by(phylum) %>% 
  summarise_each(funs(mean, sd, std.error)) %>% 
  arrange(desc(mean))

Top5_13 <- RelAbundanceMean13 %>% slice(1:5) 
Top5_16 <- RelAbundanceMean16 %>% slice(1:5) 

RelAbundance13 %>% filter(phylum %in% Top5_13$phylum) %>% 
  SinaGraph("phylum", "RelAbundance")+ 
  scale_colour_viridis_d() +theme_bw(base_size=14) -> TopPhylaSina13

RelAbundance16 %>% filter(phylum %in% Top5_16$phylum) %>% 
  SinaGraph("phylum", "RelAbundance")+ 
  scale_colour_viridis_d() + theme_bw(base_size=14) -> TopPhylaSina16


# alpha diversity ---------------------------------------------------------

## Alpha 
ps.even13 <- evenness(ps13_2, index = "all") 
ps.even16 <- evenness(ps16_2, index = "all") 

ps13.meta <- meta(ps13_2)
ps13.meta$simpson <- ps.even13$simpson 
hist(ps13.meta$simpson) #distribution of diversity roughly normal 

ps16.meta <- meta(ps16_2)
ps16.meta$simpson <- ps.even16$simpson 
hist(ps16.meta$simpson) #distribution of diversity roughly normal 

# test for normality 
shapiro.test(ps13.meta$simpson) # p=0.09002 - data normal 
qqnorm(ps13.meta$simpson) # okay 

shapiro.test(ps16.meta$simpson) # p=0.22 - normal
qqnorm(ps16.meta$simpson) # good 

# conditioning on variable type 
# create a list of pairwise comaprisons
smtype13 <- levels(as.factor(ps13.meta$Age)) # get the variables
smtype16 <- levels(as.factor(ps16.meta$Season)) # get the variables


Simpson13<- t.test(ps13.meta$simpson[ps13.meta$age_class=="L"], 
                   ps13.meta$simpson[ps13.meta$age_class=="A"])
Simpson16<- t.test(ps16.meta$simpson[ps16.meta$sample_month=="4"], 
                   ps16.meta$simpson[ps16.meta$sample_month=="8"])


# simpson plot objects 
alpha13 <- ps13.meta %>% 
  ggviolin(., x="Age", y="simpson", 
           add = "boxplot", fill = "Age", 
           #palette = c("#FFA54F", "#AB82FF"),
           palette = c(MicroColoursCat[[1]], MicroColoursCat[[2]]), 
           alpha=0.6, add.params = list(alpha=0.6), 
           legend = "top") +
  stat_compare_means(method = "t.test", 
                     aes(label=..p.signif..), 
                     label.x=1.5, label.y=0.4,
                     size=7)

alpha16 <- ps16.meta %>% 
  ggviolin(., x = "Season", y = "simpson",
           add = "boxplot", fill = "Season", 
           #palette = c("#3A5FCD", "#66CDAA"),
           palette = c(MicroColoursCat[[3]], MicroColoursCat[[4]]),
           alpha=0.75, add.params = list(alpha=0.75), 
           legend = "top") +
  stat_compare_means(method = "t.test", 
                     aes(label=..p.signif..), 
                     label.x=1.5, label.y=0.35,
                     size=7)

require(patchwork) 
alpha13+alpha16


# beta diversity ----------------------------------------------------------

ps13.rel <- microbiome::transform(ps13_2, "compositional")
bx.ord_pcoa_bray13 <- ordinate(ps13.rel, "PCoA", "bray")
bx.ord_nmds_bray13 <- ordinate(ps13.rel, "NMDS", "bray")

ps16.rel <- microbiome::transform(ps16_2, "compositional")
bx.ord_pcoa_bray16 <- ordinate(ps16.rel, "PCoA", "bray")
bx.ord_nmds_bray16 <- ordinate(ps16.rel, "NMDS", "bray")


# Axis 1 and 2 are of interest.
beta.ps13 <- plot_ordination(ps13.rel, 
                             bx.ord_pcoa_bray13, 
                             color="Age",
                             #label = "sample_id",
) + 
  geom_point(aes(shape = age_class), size= 4) + 
  scale_colour_manual(values=c(MicroColoursCat[[1]], MicroColoursCat[[2]]))+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0, size = 12),
        legend.position = "top")+
  stat_ellipse()+
  guides(shape=FALSE)

beta.ps13


beta.ps16 <- plot_ordination(ps16.rel, 
                             bx.ord_pcoa_bray16, 
                             color="Season",
                             #label = "sample_id",
) + 
  geom_point(aes(shape = sample_month), size= 4) + 
  scale_colour_manual(values=c(MicroColoursCat[[3]], MicroColoursCat[[4]]))+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0, size = 12),
        legend.position = "top")+
  stat_ellipse()+
  guides(shape=FALSE)

beta.ps16 


## PERMANOVA
metadf.bx13 <- data.frame(sample_data(ps13.rel))
bray_ps.bxn13 <- phyloseq::distance(physeq = ps13.rel, method = "bray")

metadf.bx16 <- data.frame(sample_data(ps16.rel))
bray_ps.bxn16 <- phyloseq::distance(physeq = ps16.rel, method = "bray")

set.seed(999)
# Adonis test
library(vegan)
adonis.test13 <- adonis(bray_ps.bxn13 ~ age_class, data = metadf.bx13)
adonis.test16 <- adonis(bray_ps.bxn16 ~ sample_month, data = metadf.bx16)

adonis.test13<-clean_names(adonis.test13$aov.tab)
adonis.test16<-clean_names(adonis.test16$aov.tab)

adonis.test13
adonis.test16

# divergence --------------------------------------------------------------

## dissimilarity stuff 
b.lamb <- as.data.frame(divergence(subset_samples(ps13_2, Age == "Lamb"))) %>% 
  mutate(variable="Lamb") %>% rename(value=1)
b.adult <- as.data.frame(divergence(subset_samples(ps13_2, Age == "Adult"))) %>% 
  mutate(variable="Adult") %>% rename(value=1)
dif13 <- bind_rows(b.lamb, b.adult)

b.apr <- as.data.frame(divergence(subset_samples(ps16_2, Season == "Spring")))
b.aug <- as.data.frame(divergence(subset_samples(ps16_2, Season == "Summer")))
div_df16 <- data.frame(b.apr, b.aug)
colnames(div_df16) <- c("Spring", "Summer")
dif16<- reshape2::melt(div_df16)

Divergence13 <- t.test(b.lamb$value, b.adult$value)
Divergence16 <- t.test(div_df16$Spring, div_df16$Summer)
