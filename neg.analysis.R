library(splitstackshape)
library(Rphylip) # an R interface for Joseph Felsenstein's PHYLIP phylogeny methods program package [you will also need the non-R application PHYLIP loaded onto your computer]
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


setwd("/Volumes/GoogleDrive/My Drive/Documents/Hoffman_Metabolite_Phylogeny_Paper") 

flytree <- read.tree('fly.species.list.nwk')
plot(flytree)

raw <- read.table('Comparative metabolomics global negative raw data.csv', sep=',' , stringsAsFactors = F, header=T)
raw[1:4,1:6]
names(raw)[1] <- 'Met.sample.number'

neg.samp.ID <- read.table('Comparative metabolomics negative global centered and scaled.csv', sep=',', header=T)
as.character(neg.samp.ID$Met.sample.number) -> neg.samp.ID$Met.sample.number

table(raw$Met.sample.number %in% neg.samp.ID$Met.sample.number) # only 168 samples in BOTH data, 4 were omitted due to ambiguous IDs in the sample processing steps.
TR <-raw[!raw$Met.sample.number %in% neg.samp.ID$Met.sample.number, ] 
TR[ ,1:16] # samples X41, X42, X107 and X134 were omitted
rm(TR)

neg <- raw[raw$Met.sample.number %in% neg.samp.ID$Met.sample.number, ] 
neg.samp.ID[1:4, 1:11]
neg.samp.ID[ , 2:9] -> neg.sample.info
head(neg.sample.info)
merge(neg.sample.info, neg, by.x='Met.sample.number') -> neg
neg[1:4, 1:20]
neg[ ,-c(8)] -> neg

# NA is coded as '1' in these data.
table(neg == 1)
table(is.na(neg))
neg[neg == 1] <- NA # convert to NA
neg$Day.frozen[is.na(neg$Day.frozen)] = 1 # some '1's in Day.frozen, which have been converted to NA, so switch these back

# number of replicates per: species, genotype, sex, etc:
table(neg$Line, neg$Sex) # 2 to 3 replicates per genotype and sex

# add species labels to untargeted data that match the tip labels in the phylogeny:
separate(neg, 'Line', into=c('species', 'strain'), sep='-') -> neg
neg[1:4, 1:10]
table(neg$species)
neg$species[neg$species == 'Ana'] <- 'ananassae'
neg$species[neg$species == 'Ere'] <- 'erecta'
neg$species[neg$species == 'Mel'] <- 'melanogaster'
neg$species[neg$species == 'MJ'] <- 'mojavensis'
neg$species[neg$species == 'MJCI'] <- 'mojavensis' # Jessica confimed that MJ and MJCI are the same, 'mojavensis'
neg$species[neg$species == 'Per'] <- 'persimilis'
neg$species[neg$species == 'Pse'] <- 'pseudoobscura'
neg$species[neg$species == 'Sec'] <- 'sechellia'
neg$species[neg$species == 'Sim'] <- 'simulans'
neg$species[neg$species == 'Vir'] <- 'virilis'
neg$species[neg$species == 'Wil'] <- 'willistoni'
neg$species[neg$species == 'Yak'] <- 'yakuba'
table(neg$species)
##############

write.table(neg, file='Comparative metabolomics.untargeted neg mode.csv', sep=',') 

######## transform mz data:
neg[ , 9:ncol(neg)] -> neg.mz
apply(neg.mz, 2, log) -> log.neg.mz # log transform by metabolite
apply(log.neg.mz, 1, scale) -> sc.neg.mz # scale by sample
data.frame(t(sc.neg.mz)) -> sc.neg.mz 
pairs(sc.neg.mz[ ,c(1:5)], pch=16)
names(sc.neg.mz) <- names(neg.mz)
neg[1:4, 1:10]
# replace mz data in neg with sc.neg.mz
neg[ ,9:ncol(neg)] <- sc.neg.mz
write.table(neg, file='Comparative metabolomics.untargeted neg mode logged_scaled.csv', sep=',') 

dim(log.neg.mz)
log.neg <- neg
log.neg[ ,9:ncol(log.neg)] <- log.neg.mz
pairs(log.neg.mz[ ,c(1:5)], pch=16)
# save the logged-only data
write.table(log.neg, file='Comparative metabolomics.untargeted neg mode logged-only.csv', sep=',') 
########################################

rm(raw)
rm(log.neg)
rm(log.neg.mz)
rm(sc.neg.mz)
rm(neg.mz)
rm(neg)
rm(neg.samp.ID)
rm(neg.sample.info)

read.table('Comparative metabolomics.untargeted neg mode logged-only.csv', sep=',', header = T, stringsAsFactors = F) -> neg  

# load info about each feature:
read.table('Negative ion mode metabolite info.csv', sep=',', header = T, stringsAsFactors = F) -> neg.feature.dat
head(neg.feature.dat)

neg[1:4, 1:10]
names(neg)[9:ncol(neg)] -> neg.features
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

########################################
# plot signal over missingness:
table(is.na(neg))
colSums(is.na(neg[ ,9:ncol(neg)])) -> nas
hist(nas)
1-(nas/168) -> completeness
hist(completeness, col='grey', 20, las=1)
length(completeness) - sum(completeness != 1) # only 5 complete features
sd(completeness)
hist(nas, 20, col='grey', main='missing values by mz feature (n=168 samples)\nthere are a lot of missing data\n~1058 of 2000 features are missing from more than half of the samples', las=1)
#####################################################

as.matrix(neg[ ,9:ncol(neg)]) -> neg.mat
dim(neg.mat)
table(colSums(is.na(neg.mat))/nrow(neg.mat) <= 0.1) # 102 metabolites with less than 10% missing values


## using strain or species-level means to evaluate the level of missing data within/between species and strains:
neg[1:4,1:10]
rownames(neg.mat) <- neg$ID
aggregate(neg.mat, by=list(rownames(neg.mat)), mean, na.rm=T) -> x
rownames(x) <- x$Group.1
as.matrix(x[ ,-c(1)]) -> x
table(colSums(is.na(x))==0)  # 164 features measured at least once in each strain x sex combo.

rownames(neg.mat) <- paste(neg$species, neg$Sex)
aggregate(neg.mat, by=list(rownames(neg.mat)), mean, na.rm=T) -> x
rownames(x) <- x$Group.1
as.matrix(x[ ,-c(1)]) -> x
table(colSums(is.na(x))==0) #693 features detected at least once in both sexes of all species 


