
library(splitstackshape)
library(ape) # Analysis of Phylogenetics and Evolution - R phylogeny tools
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

flytree <- read.tree('fly.species.list.nwk') 
plot(flytree)

raw <- read.table('Comparative metabolomics global positive raw data.csv', sep=',' , stringsAsFactors = F, header=T)
raw[1:4,1:6]
raw[ ,1] # last row is a QC sample, remove it:
raw <- raw[-(172), ]
names(raw)[1] <- 'Met.sample.number'
raw$Met.sample.number <- as.numeric(gsub('X', '', raw$Met.sample.number))
raw[1:4,1:6]

sampID <- read.table('globalSampleInfo.csv', sep=',', header=T, stringsAsFactors = T) # column 'Raftery' has the metabolomics sample identification to match the LC-TOR-MS data 
head(sampID)
str(sampID)

table(raw$Met.sample.number %in% sampID$Raftery) # only 167 samples in BOTH data, 4 were omitted due to ambiguous IDs in the sample processing steps.
TR <-raw[!raw$Met.sample.number %in% sampID$Raftery, ] 
TR[ ,1:16] # samples 41, 42, 107 and 134 were omitted
rm(TR)

pos <- merge(sampID, raw, by.x='Raftery', by.y='Met.sample.number', all.x = F)
pos[1:4,1:14]
pos <- pos[ ,-c(7:10)] # remove the 'sppLevel_' columns as they refer to nodes in an outdated phylogeny 
pos[1:4,1:10]

# NA is coded as '1' in these data.
table(pos[ ,10:ncol(pos)] == 1) 
table(is.na(pos[ ,10:ncol(pos)])) # conformed, no NA in raw data
pos[ ,10:ncol(pos)][pos[ ,10:ncol(pos)] == 1] <- NA # convert 1 to NA

# number of replicates per: species, genotype, sex, etc:
table(pos$Line, pos$Sex) # 2 to 3 replicates per genotype and sex


######## transform mz data:
pos[ , 10:ncol(pos)] -> pos.mz
apply(pos.mz, 2, log) -> log.pos.mz # log transform by metabolite
apply(log.pos.mz, 1, scale) -> sc.pos.mz # scale by sample
data.frame(t(sc.pos.mz)) -> sc.pos.mz 
pairs(sc.pos.mz[ ,c(1:5)], pch=16)
names(sc.pos.mz) <- names(pos.mz)
pos[1:4, 1:10]
# replace mz data in pos with sc.pos.mz
pos[ ,10:ncol(pos)] <- sc.pos.mz
write.table(pos, file='Comparative metabolomics.untargeted pos mode logged_scaled.csv', sep=',') 

dim(log.pos.mz)
log.pos <- pos
log.pos[ ,10:ncol(log.pos)] <- log.pos.mz
pairs(log.pos.mz[ ,c(1:5)], pch=16)
# save the logged-only data
write.table(log.pos, file='Comparative metabolomics.untargeted pos mode logged-only.csv', sep=',') 
########################################

###################################################################################################################
pos <- read.table('Comparative metabolomics.untargeted pos mode logged-only.csv', sep=',' , stringsAsFactors = F, header=T)

# load info about each feature:
pos.feature.dat <- read.table('Positive ion mode metabolite info.csv', sep=',', header = T, stringsAsFactors = F) 
head(pos.feature.dat)

pos[1:4, 1:10]
names(pos)[10:ncol(pos)] -> pos.features
gsub("X", "",  pos.features) -> pos.features

### the code below fixes some of the un-matching feature names (a common problem is that some of them are given two decimals in their names):
strsplit(pos.features, split='\\.') -> x
head(x)
numeric() -> fixed.masses
for(i in 1:length(x)){
  paste(x[[i]][1], x[[i]][2], sep='.') -> fixed.masses[i]
}
fixed.masses
table(fixed.masses %in% pos.feature.dat$Mass)
fixed.masses -> pos.feature.masses

pos[1:5, 1:10]
write.table(pos, file='Comparative metabolomics.untargeted pos mode logged-only.csv', sep=',') 
#################################################################################################
