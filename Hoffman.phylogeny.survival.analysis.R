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
flytree <- read.tree('fly.species.list.nwk')
plot(flytree)
axisPhylo(side = 1, root.time = NULL, backward = TRUE)

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


mean(table(ls$Species, ls$Genotype, ls$Sex)[table(ls$Species, ls$Genotype, ls$Sex) >0])
sd(table(ls$Species, ls$Genotype, ls$Sex)[table(ls$Species, ls$Genotype, ls$Sex) >0])

ls$Event <- 1 ## there are no censored flies in these data, so we could analyze the lifespans without survival modeling.  
############################################################################################################

survfit(Surv(ls$Hour, ls$Event) ~ Species + Sex, data = ls)
hist(ls$Hours/24, col='grey', 30, border=0, las=1, xlab='individual lifespan (days)', main='study-wide lifespans')

## lifespan summary by strain (a.k.a. Genotype):
strain_kap <- survfit(Surv(ls$Hour, ls$Event) ~ Genotype + Sex, data = ls)
strain_kap <- summary(strain_kap, print.rmean=TRUE)
strain.stats <- data.frame(strain_kap$table)
head(strain.stats)
strain.stats$info <- rownames(strain.stats)

cSplit(strain.stats, 'info', sep=',') -> strain.stats
cSplit(strain.stats, c('info_1', 'info_2'), sep='=') -> strain.stats
head(strain.stats)
strain.stats[ ,c(11,13,1,5:9)] -> strain.stats
colnames(strain.stats)[1:2] <- c('strain', 'sex')
head(strain.stats)

write.table(strain.stats, file='mean.lifespan.by.strain.txt', row.names=F, quote=F) 


## lifespan summary by species:
kaplan <- survfit(Surv(ls$Hour, ls$Event) ~ Species + Sex, data = ls)
summary(kaplan, print.rmean=TRUE) -> kap.sum
data.frame(kap.sum$table) -> ls.stats
head(ls.stats)
rownames(ls.stats) -> ls.stats$info

cSplit(ls.stats, 'info', sep=',') -> ls.stats
cSplit(ls.stats, c('info_1', 'info_2'), sep='=') -> ls.stats
head(ls.stats)
ls.stats[ ,c(11,13,1,5:9)] -> ls.stats
colnames(ls.stats)[1:2] <- c('species', 'sex')
head(ls.stats)

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



