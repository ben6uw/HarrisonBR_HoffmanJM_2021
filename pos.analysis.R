
library(ape) # Analysis of Phylogenetics and Evolution - R phylogeny tools
library(phytools) # Phylogenetic Tools for Comparative Biology (and Other Things)
library(caper) # Comparative Analyses of Phylogenetics and Evolution in R
library(geiger) # Analysis of Evolutionary Diversification

library(nlme) # Linear and Nonlinear Mixed Effects Models

library(tidyr) # data wrangling
library(plyr) # data wrangling
library(splitstackshape) # data wrangling
library(impute) # knn imputation
library(DMwR) # has knn imputation function

require(graphics) # for plotting trees more easily
library(corrplot) # visualizing correlation matrices
library(gplots)

setwd("/Volumes/GoogleDrive/My Drive/Documents/Hoffman_Metabolite_Phylogeny_Paper/Harrison_Hoffman_GitHub") 

###################################
# SKIP to line 100 to load data (d.f. pos), everything prior to that is massaging the raw data.
###################################

read.table('Comparative metabolomics global positive raw data.csv', sep=',', header=T) ->pos
pos$X <- as.character(pos$X)

read.table('Comparative metabolomics positive global centered and scaled.csv', sep=',', header=T) -> pos.samp.ID
as.character(pos.samp.ID$Met.sample.number) -> pos.samp.ID$Met.sample.number

par(mar=c(5,5,5,5))

table(pos$X %in% pos.samp.ID$Met.sample.number)# only 167 samples in BOTH data, 5 were omitted due to ambiguous IDs in the sample processing steps.
pos[!pos$X %in% pos.samp.ID$Met.sample.number, ] -> TR
TR[ ,1:16] # samples X41, X42, X107 and X134 were omitted

pos[pos$X %in% pos.samp.ID$Met.sample.number, ] -> pos
pos.samp.ID[1:4, 1:11]
pos.samp.ID[ , 2:9] -> pos.sample.info
merge(pos.sample.info, pos,by.x='Met.sample.number', by.y = 'X') -> pos
pos[ ,-c(8)] -> pos
str(pos)

# NA is coded as '1' in these data.
table(pos == 1)
table(is.na(pos))
pos[pos == 1] <- NA # convert to NA
pos$Day.frozen[is.na(pos$Day.frozen)] = 1 # some '1's in Day.frozen, which have been converted to NA, so switch these back

# number of replicates per: species, genotype, sex, etc:
table(pos$Line, pos$Sex) # 2 to 3 replicates per genotype and sex

pos[1:6, 1:10]
pairs(pos[ ,10:13])
pairs(log(pos[ ,10:13])) 
pairs(log(pos[ ,10:13]), col=as.factor(pos$Sex), pch=19) 
pairs(log(pos[ ,10:13]), col=as.factor(pos$Line), pch=19) 

# add species labels to untargeted data that match the tip labels in the phyloeny:
separate(pos, 'Line', into=c('species', 'strain'), sep='-') -> pos
pos[1:4, 1:10]
table(pos$species)
pos$species[pos$species == 'Ana'] <- 'ananassae'
pos$species[pos$species == 'Ere'] <- 'erecta'
pos$species[pos$species == 'Mel'] <- 'melanogaster'
pos$species[pos$species == 'MJ'] <- 'mojavensis'
pos$species[pos$species == 'MJCI'] <- 'mojavensis' # Jessica confimed that MJ and MJCI are the same, 'mojavensis'
pos$species[pos$species == 'Per'] <- 'persimilis'
pos$species[pos$species == 'Pse'] <- 'pseudoobscura'
pos$species[pos$species == 'Sec'] <- 'sechellia'
pos$species[pos$species == 'Sim'] <- 'simulans'
pos$species[pos$species == 'Vir'] <- 'virilis'
pos$species[pos$species == 'Wil'] <- 'willistoni'
pos$species[pos$species == 'Yak'] <- 'yakuba'
table(pos$species)
pairs(log(pos[ ,10:13]), col=as.factor(pos$species), pch=19) 
##############

######## transform mz data:
pos[ , 9:ncol(pos)] -> pos.mz
apply(pos.mz, 2, log) -> log.pos.mz # log transform by metabolite
apply(log.pos.mz, 1, scale) -> sc.pos.mz # scale by sample
data.frame(t(sc.pos.mz)) -> sc.pos.mz 
pairs(sc.pos.mz[ ,c(1:5)], pch=16)
names(sc.pos.mz) <- names(pos.mz)
pos[1:4, 1:10]
# replace mz data in pos with sc.pos.mz
pos[ ,9:ncol(pos)] <- sc.pos.mz
write.table(pos, file='Comparative metabolomics.untargeted pos mode logged_scaled.csv', sep=',') 

dim(log.pos.mz)
log.pos <- pos
log.pos[ ,9:ncol(log.pos)] <- log.pos.mz
pairs(log.pos.mz[ ,c(1:5)], pch=16)
# save the logged-only data
write.table(log.pos, file='Comparative metabolomics.untargeted pos mode logged-only.csv', sep=',') 
########################################


###################################################################################################################
##### SKIP to HERE to load data:
read.table('Comparative metabolomics.untargeted pos mode logged-only.csv', sep=',' , stringsAsFactors = F, header=T) -> pos

# load info about each feature:
read.table('Positive ion mode metabolite info.csv', sep=',', header = T, stringsAsFactors = F) -> pos.feature.dat
head(pos.feature.dat)

pos[1:4, 1:10]
names(pos)[9:ncol(pos)] -> pos.features
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
