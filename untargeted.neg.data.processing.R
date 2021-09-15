library(splitstackshape)
library(ape) # Ananlysis of Phylogenetics and Evolution - R phylogeny tools
library(phytools) # Phylogenetic Tools for Comparative Biology (and Other Things)
library(corrplot) # visualizing correlation matricies
library(geiger) # Analysis of Evolutionary Diversification
library(nlme) # Linear and Nonlinear Mixed Effects Models
library(tidyr) # data wrangling
library(plyr) # data wrangling
library(caper) # Comparative Analyses of Phylogenetics and Evolution in R
library(impute) # knn imputation
library(DMwR) # has knn imputation function
require(graphics) # for plotting trees more easily
library(sjstats) # a package for investigating the results of mixed models
library(lme4)
library(lmerTest)
library(gplots)


setwd("/Volumes/GoogleDrive/My Drive/Documents/Hoffman_Metabolite_Phylogeny_Paper/Harrison_Hoffman_GitHub") 
dir()

raw <- read.table('Comparative metabolomics global negative raw data.csv', sep=',' , stringsAsFactors = F, header=T)
raw[1:4,1:6]
names(raw)[1] <- 'Met.sample.number'
raw$Met.sample.number <- as.numeric(gsub('X', '', raw$Met.sample.number))
raw[1:4,1:6]

sampID <- read.table('globalSampleInfo.csv', sep=',', header=T, stringsAsFactors = T) # column 'Raftery' has the metabolomics sample identification to match the LC-TOR-MS data 
head(sampID)
str(sampID)

table(raw$Met.sample.number %in% sampID$Raftery) # only 168 samples in BOTH data, 4 were omitted due to ambiguous IDs in the sample processing steps.
TR <-raw[!raw$Met.sample.number %in% sampID$Raftery, ] 
TR[ ,1:16] # samples 41, 42, 107 and 134 were omitted
rm(TR)

sampID[1:4, ]
neg <- merge(sampID, raw, by.x='Raftery', by.y='Met.sample.number', all.y = F)
neg[1:4, 1:20]
neg <- neg[ ,-c(7:10)] # remove the 'sppLevel_' columns as they refer to nodes in an outdated phylogeny 

# NA is coded as '1' in these data.
table(neg[ ,10:ncol(neg)] == 1)
table(is.na(neg[ ,10:ncol(neg)])) # currently no NA values are coded in the LC-TOR-MS data
neg[ ,10:ncol(neg)][neg[ ,10:ncol(neg)] == 1] <- NA # convert to NA
table(is.na(neg[ ,10:ncol(neg)]))

# number of replicates per: species, genotype, sex, etc:
table(neg$Line, neg$Sex) # 1 to 3 replicates per genotype and sex

neg[1:4, 1:20]

write.table(neg, file='Comparative metabolomics.untargeted neg mode.csv', sep=',') 

######## transform mz data:
neg[ , 10:ncol(neg)] -> neg.mz
apply(neg.mz, 2, log) -> log.neg.mz # log transform by metabolite
apply(log.neg.mz, 1, scale) -> sc.neg.mz # scale by sample
data.frame(t(sc.neg.mz)) -> sc.neg.mz 
pairs(sc.neg.mz[ ,c(1:5)], pch=16)
names(sc.neg.mz) <- names(neg.mz)
neg[1:4, 1:10]
# replace mz data in neg with sc.neg.mz
neg[ ,10:ncol(neg)] <- sc.neg.mz
write.table(neg, file='Comparative metabolomics.untargeted neg mode logged_scaled.csv', sep=',') 

dim(log.neg.mz)
log.neg <- neg
log.neg[ ,10:ncol(log.neg)] <- log.neg.mz
pairs(log.neg.mz[ ,c(1:5)], pch=16)
# save the logged-only data
write.table(log.neg, file='Comparative metabolomics.untargeted neg mode logged-only.csv', sep=',') 
########################################

rm(list = ls())

flytree <- read.tree('fly.species.list.nwk') 
plot(flytree)

read.table('Comparative metabolomics.untargeted neg mode logged-only.csv', sep=',', header = T, stringsAsFactors = F) -> neg  

# load info about each feature:
read.table('Negative ion mode metabolite info.csv', sep=',', header = T, stringsAsFactors = F) -> neg.feature.dat
head(neg.feature.dat)

neg[1:4, 1:10]
names(neg)[10:ncol(neg)] -> neg.features
gsub("X", "",  neg.features) -> neg.features

### the code below fixes some of the un-matching feature names (a common problem is that some of them are given two decimals in their names):
strsplit(neg.features, split='\\.') -> x
head(x)
numeric() -> fixed.masses
for(i in 1:length(x)){
  paste(x[[i]][1], x[[i]][2], sep='.') -> fixed.masses[i]
}
fixed.masses
table(fixed.masses %in% neg.feature.dat$Mass)
fixed.masses[!(fixed.masses %in% neg.feature.dat$Mass)]
fixed.masses[!(fixed.masses %in% neg.feature.dat$Mass)] <- '322.0'
fixed.masses -> neg.feature.masses

neg[1:5, 1:10]
write.table(neg, file='Comparative metabolomics.untargeted neg mode logged-only.csv', sep=',') 

neg <- read.table('Comparative metabolomics.untargeted neg mode logged-only.csv', sep=',', header = T, stringsAsFactors = F)

########################################
# plot signal over missingness:
table(is.na(neg))
colSums(is.na(neg[ ,9:ncol(neg)])) -> nas
hist(nas)
1-(nas/168) -> completeness
hist(completeness, col='grey', 20, las=1)
length(completeness) - sum(completeness != 1) # 15 complete features
sd(completeness)
hist(nas, 20, col='grey', main='missing values by mz feature (n=168 samples)\nthere are a lot of missing data\n~1058 of 2000 features are missing from more than half of the samples', las=1)
#####################################################

as.matrix(neg[ ,10:ncol(neg)]) -> neg.mat
dim(neg.mat)
table(colSums(is.na(neg.mat))/nrow(neg.mat) <= 0.1) # 102 metabolites with less than 10% missing values

## using strain or species-level means to evaluate the level of missing data within/between species and strains:
aggregate(neg.mat, by=list(neg$ID), mean, na.rm=T) -> x
x[1:4,1:4]
as.matrix(x[ ,-c(1)]) -> x
table(colSums(is.na(x))==0)  # 164 features measured at least once in each strain x sex combo.

x <- aggregate(neg.mat, by=list(species = neg$Species, sex = neg$Sex), mean, na.rm=T)
x[1:4,1:4]
table(colSums(is.na(x[ ,-c(1, 2)]))==0) #693 features detected at least once in both sexes of all species 

