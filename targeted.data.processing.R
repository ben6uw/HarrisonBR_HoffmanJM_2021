## Working with Jessica Hoffman's phylogenetic metabolome data:

library(ape) # Ananlysis of Phylogenetics and Evolution - R phylogeny tools
library(phytools) # Phylogenetic Tools for Comparative Biology (and Other Things)
library(corrplot) # visualizing correlation matricies
library(geiger) # Analysis of Evolutionary Diversification
library(nlme) # Linear and Nonlinear Mixed Effects Models
library(plyr) # data wrangling
library(caper) # Comparative Analyses of Phylogenetics and Evolution in R
require(graphics) # for plotting trees more easily
library(sjstats) # a package for investigating the results of mixed models
library(caret)
library(corrr)
library(lme4)
library(lmerTest)
library(gplots)
library(gridExtra)
library(MuMIn)
library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)
library(tidyverse)
library(viridis)# color palettes

gc()

setwd("/Users/ben/Google Drive/My Drive/Documents/Hoffman_Metabolite_Phylogeny_Paper/Harrison_Hoffman_GitHub") 

###################################
# SKIP to line 93 to load data (d.f. dat), everything prior to that is processing the raw data.
###################################

dat <- read.table('Comparative metabolomics targeted raw data.csv', header=T, stringsAsFactors = F, sep=',')
targeted_mz_info <- as.data.frame(read_csv("targeted.mz.info.csv"))

head(dat[ ,1:20])
dat[ ,1:14] -> sample.info
head(sample.info)

names(dat)[15:ncol(dat)] <- targeted_mz_info$metabolite[match(names(dat)[15:ncol(dat)], targeted_mz_info$mz.name)] # replace metabolite codes with metabolite names

# identify internal standards (13C compounds added to the samples by the Raftery lab):
stds <- dat[ ,colnames(dat) == 'internal_standard']
names(stds) <- c('internal_std_A', 'internal_std_B')
plot(stds$internal_std_A , stds$internal_std_B, pch=19, main='2 internal standards', log='xy')

# remove the standards, and later add them back as special variables 
dat[ ,colnames(dat) != 'internal_standard'] -> dat

sample.info <- cbind(sample.info, stds)

# remove several of the metabolites based on Danijel's post-hoc analysis of highly-correlated metabolite pairs (Danijel called these particular data 'ghost peaks' of other metabolites, meaning that enough of their signal comes from other features in the LC that they shoudl be excluded):  
ghost_peaks <- c('4-Trimethylammoniobutanoate', 'Sorbitol', 'Xanthosine', 'Oxypurinol', 'Uridine', 'Trimethylamine', 'Valine', 'Adenine', 'Pipecolic acid')

table(ghost_peaks %in% colnames(dat))
dat <- dat[ ,!colnames(dat) %in% ghost_peaks]

# normalize the data
head(dat)
as.matrix(dat[ ,15:ncol(dat)]) -> m2 # pull the mz data as a matrix
apply(m2, 2, log) -> logmz # log each metabolite
pairs(t(logmz)[ ,1:5]) # metabolite data across samples are highly correlated (probably driven by mz features that are either high in every sample, or low in every sample)
pairs(logmz[ ,1:5]) # metabolites are not (generally) highly correlated with each other; but see below

apply(logmz, 1, function(x) scale(x, scale=F)) -> norm.dat # center (but don't also scale) by sample
norm.dat[1:5,1:5]
t(norm.dat) -> norm.dat
norm.dat[1:5,1:4]
head(dat)
colnames(norm.dat) <- colnames(dat[ ,15:ncol(dat)])
pairs(norm.dat[,1:5])
pairs(t(norm.dat)[,1:5])

par(mfrow=c(1,3))
par(mar=c(2,4,2,2))
boxplot(t(m2), ylab='raw')
boxplot(t(logmz), ylab='log.mz')
boxplot(t(norm.dat), ylab='log.mz centered.by.sample')

par(mfrow=c(1,1))
par(mar=c(5,4,4,2))

norm.dat[1:5,1:5]
cbind(sample.info, norm.dat) -> dat
write.table(dat, 'comparative.mz.targeted.data.log_centered', row.names=F, quote=F, sep='\t')

cbind(sample.info, logmz) -> logged.dat
write.table(logged.dat, 'comparative.mz.targeted.data.log_only', row.names=F, quote=F, sep='\t')

rm(list = ls())

################# load data and continue
mycolors = c("#FF0000", '#FF6347', "#AFEEEE", "#87CEFA")
mycol4 = c('#F4A460', "#CD6600", "#87CEFA", "#104E8B")

read.table('comparative.mz.targeted.data.log_only', header=T, sep='\t') -> dat
head(dat) # a leading X is added to the metabolite names that begin with a number

names(dat)[1:18]
dat[, 17:ncol(dat)] -> mzdat
mzdat[ ,colSums(is.na(mzdat))==0] -> dat.no.miss
mzs <- colnames(dat.no.miss)
cbind(dat[ ,1:16], dat.no.miss) -> dat.no.miss

mzdat[1:4,1:4]
mets <- as.matrix(mzdat)
mets <- mets[ ,colSums(is.na(mets)) == 0] # rid matrix of mzs with missing data
centered_mets <- t(scale(t(mets), scale=F)) # center, but don't scale metabolites by sample
cbind(dat[ ,1:16], centered_mets) -> centered_dat

####################################################################################################################################
# examine the potential relationship between the samples, the metabolome and/or individual metabolites and the two internal 13C standards that were added to the samples.  Note, the internal standards are not log-normalized in this analysis, as they are fairly normal to begin with:
####################################################################################################################################

# standards vs. metabolome (PCA):
pca <- prcomp(centered_mets, scale=F)

par(mfrow=c(1,1))
plot(pca$x, col=mycol4[as.numeric(as.factor(paste(dat$sex, dat$age.days)))], pch=16)
legend('bottomleft', legend=c(levels(as.factor(paste(dat$sex, dat$age.days)))), pch=16, col=mycol4)

pca <- pca$x

par(mfrow=c(2,2))
plot(pca[ ,1] ~ dat$internal_std_A, pch=16, main='standard A')
plot(pca[ ,2] ~ dat$internal_std_A, pch=16, main='standard A')
plot(pca[ ,1] ~ dat$internal_std_B, pch=16, main='standard B') 
summary(lm(pca[ ,2] ~ dat$internal_std_A)) # standard_A has an inadvertent relationship with the samples by species and sex (see below), which may drive this relationship with PC2
abline(lm(pca[ ,2] ~ dat$internal_std_A)) 
summary(lm(pca[ ,1] ~ dat$internal_std_B)) # standard_B has an inadvertent relationship with the samples by species and sex (see below), which may drive this relationship with PC1
abline(lm(pca[ ,1] ~ dat$internal_std_B)) 
plot(pca[ ,2] ~ dat$internal_std_B, pch=16, main='standard B')
summary(lm(pca[ ,2] ~ dat$internal_std_B)) # standard_B has an inadvertent relationship with the samples by species and sex (see below), which may drive this relationship with PC2
abline(lm(pca[ ,2] ~ dat$internal_std_B))


# individual metabolites
tmp <- matrix(nrow=nrow(mets), ncol=ncol(mets))
for(i in 1:ncol(mets)){
  tmp[ ,i] <- lm(centered_mets[ ,i] ~ dat$species + dat$strain * dat$sex + dat$age.days)$residuals # analyze residuals after species, strain * sex, and age
}
colnames(tmp) <- colnames(mets)

par(mfrow=c(3,2))
hist(mets, 30)
hist(moments::kurtosis(mets), 30)
hist(centered_mets, 30)
hist(moments::kurtosis(centered_mets), 30)
hist(tmp, 30) # residuals from species, strain * sex and age reduced sample variation substantially, they appear more gausian, however I wonder if the
hist(moments::kurtosis(tmp), 30) # kurtosis of gausian data = 3, and higher values indicate underdispersion so, as seen in the histograms, the residuals reduce varaition (too much?)

internal_standards_effects <- matrix(nrow=ncol(mets), ncol=2)
for(i in 1:ncol(mets)){
  m <- anova(lm(tmp[ ,i] ~ dat$internal_std_A + dat$internal_std_B)) # regress each standard over the species, strain * sex, age residuals
  internal_standards_effects[i, ] <- m$`Pr(>F)`[1:2]
}
rownames(internal_standards_effects) <- colnames(mets)


# relationships with standard_A:
names(which.min(internal_standards_effects[ ,1])) # metabolite at Pmin

par(mfrow=c(2,2))
plot(tmp[ ,names(which.min(internal_standards_effects[ ,1]))] ~ dat$internal_std_A, pch=16, ylab=paste('residual', names(which.min(internal_standards_effects[ ,1])))) # relationship for mz at min P value is driven by outlier
plot(centered_mets[ ,names(which.min(internal_standards_effects[ ,1]))] ~ dat$internal_std_A, pch=16, ylab=paste('centered', names(which.min(internal_standards_effects[ ,1])))) 

# relationships with standard_B:
plot(tmp[ ,names(which.min(internal_standards_effects[ ,2]))] ~ dat$internal_std_B, pch=16, ylab=paste('residual', names(which.min(internal_standards_effects[ ,2])))) # 
plot(centered_mets[ ,names(which.min(internal_standards_effects[ ,2]))] ~ dat$internal_std_B, pch=16, ylab=paste('centered', names(which.min(internal_standards_effects[ ,2])))) # relationship for the centered mz (prior to taking residuals of sp, strain * sex, age) and standard B is 'tight'.  Interestingly, the residuals remove this effect.   

anova(lm(dat$internal_std_B ~ dat$species + dat$strain * dat$sex + dat$age.day)) # effects of sex and species on the internal standard_B

par(mfrow=c(2,2))    
plot(dat$internal_std_B ~ as.factor(dat$sex))
stripchart(dat$internal_std_B ~ as.factor(dat$sex), add=T, vertical=T, pch=19)
plot(dat$internal_std_B ~ as.factor(dat$species))
stripchart(dat$internal_std_B ~ as.factor(dat$species), add=T, vertical=T, pch=19)

y <- centered_mets[ ,names(which.min(internal_standards_effects[ ,2]))]
plot(y ~ dat$internal_std_B, pch=16, ylab=paste('centered', names(which.min(internal_standards_effects[ ,2]))))
abline(lm(y ~ dat$internal_std_B))
plot(lm(y ~ dat$sex * dat$species)$residuals ~ dat$internal_std_B, pch=16 )

summary(lm(lm(y ~ dat$sex * dat$species)$residuals ~ dat$internal_std_B)) 
plot(lm(y ~ dat$sex * dat$species)$residuals ~ dat$internal_std_B, pch=16) # residuals of species and sex reduced, but not eliminated, the relationship between std_B and these variables, and at this point, the relationship looks driven by an outlier

rownames(internal_standards_effects)[order(internal_standards_effects[ ,2])] [1:10] # the most significant 10 

par(mfrow=c(2,2))
for(i in 1:4){
  tmp_mz <- rownames(internal_standards_effects)[order(internal_standards_effects[ ,2])] [i]
  plot(tmp[ ,tmp_mz] ~ dat$internal_std_B, pch=16, ylab=tmp_mz)
}
####################################################################################################################################
# conclusions: 
# internal standards A does not have a relationship with the experimental variables or any metabolites, except for relationships that are driven by outliers.  
# internal standard B has a relationship with species and sex, however after correcting for their effect, none of the metabolites has a significant relationship with std B
####################################################################################################################################
