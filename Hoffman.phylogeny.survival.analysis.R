## survival analysis in Jessicas Phylogenetic Data:
# analysis of traits (not metabolites) in Jessica's multi-species data:
# load and analyze longevity data:

library(survival)
library(survminer)
library(lattice)
library(lme4)
library(ggplot2)
library(splitstackshape)
library(lmerTest)
library(sjstats)
library(phytools)
library(viridis)
library(MCMCpack)
library(MCMCglmm)
library(tidyverse)
library(gplots) 
library(plyr)
require(graphics) 
library(ggpubr)

setwd("/Volumes/GoogleDrive/My Drive/Documents/Hoffman_Metabolite_Phylogeny_Paper/Harrison_Hoffman_GitHub") 

# load longevity data
ls <- read.csv('Expanded comparative longevity data.csv', stringsAsFactors = T, header=T)
head(ls)
str(ls)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### examine sample sizes:
table(ls$Sex, ls$Genotype)
range(table(ls$Sex, ls$Genotype))
mean(table(ls$Sex, ls$Genotype))
sd(table(ls$Sex, ls$Genotype))
table(ls$Species, ls$Sex)
range(table(ls$Species, ls$Sex))

ls$Event <- 1 ## there are no censored flies in these data, so we could analyze the lifespans without survival modeling.  
############################################################################################################

survfit(Surv(ls$Hour, ls$Event) ~ Species + Sex, data = ls)
hist(ls$Hours/24, col='grey', 30, border=0, las=1, xlab='individual lifespan (days)', main='study-wide lifespans')

head(ls)
ls$Genotype <- as.character(ls$Genotype)
ls$Sex <- as.character(ls$Sex)
ls$Genotype
ls$Genotype[ls$Genotype == 'mjci1002'] <- 'moj1002'
ls$Genotype[ls$Genotype == 'mjci2008'] <- 'moj2008'
ls$Sex[ls$Sex=='female'] <- 'F'
ls$Sex[ls$Sex=='male'] <- 'M'
head(ls)

## lifespan summary by strain (a.k.a. Genotype):
strain_kap <- survfit(Surv(ls$Hour, ls$Event) ~ Genotype + Sex, data = ls)
strain_kap <- summary(strain_kap, print.rmean=TRUE)
strain.stats <- data.frame(strain_kap$table)
head(strain.stats)
strain.stats$info <- rownames(strain.stats)
cSplit(strain.stats, 'info', sep=',') -> strain.stats
cSplit(strain.stats, c('info_1', 'info_2'), sep='=') -> strain.stats
head(strain.stats)
colnames(strain.stats)[c(1,5,6,7,11,13)] <- c('n.flies.ls', 'mean.ls', 'se.mean.ls', 'median.ls', 'strain', 'sex')
head(strain.stats)
strain.stats$mean.ls <- strain.stats$mean.ls/24 
strain.stats$se.mean.ls <- strain.stats$se.mean.ls/24
strain.stats$median.ls <- strain.stats$median.ls/24
strain.stats <- strain.stats[ ,c(11,13,1,5:7)]

head(strain.stats)

write.table(strain.stats, file='mean.lifespan.by.strain.txt', row.names=F, quote=F) 

## lifespan summary by species:
head(strain.stats) 
head(ls)
tmp <- ls[ ,c(1,3)]
tmp <- tmp[!duplicated(tmp), ]
tmp
ls.stats <- merge(strain.stats, tmp, by.x='strain', by.y='Genotype')
ls.stats <- aggregate(ls.stats, by=list(ls.stats$Species, ls.stats$sex), mean, na.rm=T)
head(ls.stats)
ls.stats$species <- ls.stats$Group.1
ls.stats$sex <- ls.stats$Group.2
head(ls.stats)
ls.stats <- ls.stats[ ,c(10, 4:8)]

write.table(ls.stats, "lifespan_summary_stats", sep=',', quote=F, row.names=F)
ls.stats <- read.table("lifespan_summary_stats", sep=',', header=T)

# violin plots:
ggplot(ls, aes(x=Genotype, y=Hours/24, fill=Sex)) + 
geom_violin(trim=T) + 
theme_minimal() + 
theme(axis.text.x=element_blank())+
labs(x="Species:Strain", y = "Lifespan (days)") +
theme_classic() +
theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(subset(ls, Sex == 'female'), aes(x=Genotype, y=Hours/24, fill=Species)) + 
  geom_violin(trim=T) + 
  theme_minimal() + 
  theme(axis.text.x=element_blank())+
  labs(x="Species:Strain", y = "Lifespan (days)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip() +
  annotate("text", x = 3, y = 140, label = "female", size=5) +
  ylim(0, 170)

ggplot(subset(ls, Sex == 'male'), aes(x=Genotype, y=Hours/24, fill=Species)) + 
  geom_violin(trim=T) + 
  theme_minimal() + 
  theme(axis.text.x=element_blank())+
  labs(x="Species:Strain", y = "Lifespan (days)") +
  theme_classic() +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  annotate("text", x = 3, y = 140, label = "male", size=5) +
  ylim(0, 170)


#### plot survival curves:
unique(ls$Species) -> pull.species.from.here
unique(ls$Species) -> species
paste('D.', tolower(species)) -> species
length(species)

par(mar=c(5,4,5,1))

par(mfrow=c(3,4))
for(i in 1:length(species)) {
ls[ls$Species == pull.species.from.here[i], ] -> x
as.factor(x$Sex) -> x$Sex
plot(survfit(Surv(Hours/24, Event) ~ Genotype + Sex, data=x), las=1, col=c(1,1,2,2,8,8), lwd=2, lty=c(1,2), cex.lab=1, cex.axis=1, main=species[i], cex.main=1, xlab='time (day)', ylab='survivorship')
}
hist(ls$Hours/24, col='grey', 30, border=0, las=1, xlab='individual lifespan (days)', main='study-wide lifespans')
