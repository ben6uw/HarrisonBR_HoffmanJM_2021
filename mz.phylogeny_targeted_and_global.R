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
library(lme4)
library(lmerTest)
library(gplots)
library(gridExtra)
library(MuMIn)
library(ggplot2)

mycolors = c("#FF0000", '#FF6347', "#AFEEEE", "#87CEFA")
mycol4 = c('#F4A460', "#CD6600", "#87CEFA", "#104E8B")

setwd("~/Google Drive/My Drive/Documents/Hoffman_Metabolite_Phylogeny_Paper")
flytree <- read.tree('fly.species.list.nwk')

#####################################################################
# targeted data
read.table('comparative.mz.targeted.data.log_only', header=T, sep='\t') -> dat
targeted_mz_info <- as.data.frame(read.csv("targeted.mz.info.csv"))

head(dat) # a leading X is added to the metabolite names that begin with a number
names(dat)[1:18]
dat[, 18:ncol(dat)] -> mzdat
mzdat[ ,colSums(is.na(mzdat))==0] -> dat.no.miss
mzs <- colnames(dat.no.miss)
cbind(dat[ ,1:17], dat.no.miss) -> dat.no.miss
dat.no.miss$age_sex <- paste(dat.no.miss$age, dat.no.miss$sex)

mzdat[1:4,1:4]
mets <- as.matrix(mzdat)
mets <- mets[ ,colSums(is.na(mets)) == 0] # rid matrix of mzs with missing data
centered_mets <- t(scale(t(mets), scale=F)) # center, but don't scale metabolites by sample
cbind(dat[ ,1:17], centered_mets) -> centered_dat
centered_dat$age_sex <- paste(centered_dat$age, centered_dat$sex)

dat$species -> Species
dat$strain -> Strain
dat$sex -> Sex
dat$age.days -> Age
#####################################################################

#####################################################################
# untargeted data
read.table('combined.untargeted.data.csv', sep=',', stringsAsFactors = F, header=T) -> glob
glob[1:4,1:10]
glob[ ,colnames(glob) != 'pos_mode'] -> glob # gets rid of mode name
glob[ ,colnames(glob) != 'neg_mode'] -> glob # gets rid of mode name
glob$strain <- paste(glob$species, glob$strain, sep='_')

feats <- glob[ ,9:ncol(glob)]
centered_feats <- t(scale(t(feats), scale=F)) # center,but not scale, by sample
centered_feats[1:4,1:4]
glob <- cbind(glob[ ,4:6], centered_feats)

## analyze untargeted data by sex and species
glob[1:4, 1:8]
aggregate(glob[, 4:ncol(glob)], list(glob$sex, glob$species), mean, na.rm=T) -> glob.species.means 
glob.species.means[1:5, 1:5]
names(glob.species.means)[1:2] <- c('sex', 'species')

table(colSums(is.na(glob.species.means))[3:ncol(glob.species.means)]) # 1558 features present at least once in each species AND in BOTH sexes
table(colSums(is.na(glob.species.means[glob.species.means$sex == 'females', ]))[3:ncol(glob.species.means)]) # 2366 features present at least once in each species in females
table(colSums(is.na(glob.species.means[glob.species.means$sex == 'males', ]))[3:ncol(glob.species.means)]) # 2083 features present at least once in each species in males

glob.species.means[1:6, 1:4]

glob.species.means[glob.species.means$sex == 'males', ] -> m
m[1:4,1:4]
m[ ,colSums(is.na(m)) ==0] -> m
as.matrix(m[ ,3:ncol(m)]) -> m2
rownames(m2) <- m$species
m2[1:4,1:4]

glob.species.means[glob.species.means$sex == 'females', ] -> f
f[,colSums(is.na(f)) ==0] -> f
f[1:4,1:5]
as.matrix(f[ ,3:ncol(f)]) -> f2
rownames(f2) <- f$species
m2[rownames(f2), ] -> m2
f2[1:4,1:4]

# a simple tree plot:
?stats::hclust
ftree <- hclust(dist(f2), method='ave')
plot(as.phylo(ftree), cex = 1, label.offset = 0.5)

## bootstrapping the tree
tmptree <- as.phylo(hclust(dist(f2)), method='ave')
bs <- boot.phylo(tmptree, f2, FUN = function(xx) as.phylo(hclust(dist(xx))), B=1000)
tmptree$node.label <- bs

par(mfrow=c(1,1))
par(mar=c(3,3,3,3))
par(xpd=TRUE)
plot(tmptree, main='untargeted, female species means of 2366 features')
nodelabels(pie=(tmptree$node.label/1000), cex=1, piecol=c('grey',0))



tmptree <- as.phylo(hclust(dist(m2)), method='ave')
bs <- boot.phylo(tmptree, m2, FUN = function(xx) as.phylo(hclust(dist(xx))), B=1000)
tmptree$node.label <- bs

par(mfrow=c(1,1))
par(mar=c(3,3,3,3))
par(xpd=TRUE)
plot(tmptree, main='untargeted, male species means of 2083 features')
nodelabels(pie=(tmptree$node.label/1000), cex=1, piecol=c('grey',0))

plot(flytree)


## tree-tree distance:
# scale trees relative to their total edge length to make comparable:
flytree$edge.length <- flytree$edge.length/sum(flytree$edge.length)

# permute random trees to assess significance of tree-tree distance:
distance <- numeric()
set.seed(123456)
for(i in 1:10000){
  rtree(11) -> r # make random tree with 11 tips
  r$edge.length <- r$edge.length/sum(r$edge.length) # standardize tree length
  dist.topo(flytree, r, method='score') -> distance[i]# get branch score
}

hist(distance, 20, border = 0, las=1, xlim=c(0,0.4))

# scale tree relative to its total edge length to make comparable to the scaled flytree:
tmptree <- as.phylo(hclust(dist(m2)), method='ave') # specify data here (e.g. male or female data)
tmptree$edge.length <- tmptree$edge.length/sum(tmptree$edge.length)

dist.topo(flytree, tmptree, method='score')
sum(dist.topo(flytree, tmptree, method='score') > distance) / 10000  # empirical P value for tree-tree distance

hist(distance, border=0, las=1, xlim=c(0, 0.4), xlab='branch score', main='distance flytree to randomized tree', 30)

xpd=F
tmptree <- as.phylo(hclust(dist(f2)))
tmptree$edge.length <- tmptree$edge.length/sum(tmptree$edge.length)
abline(v=dist.topo(flytree, tmptree, method='score'), col='red', lty=2)
tmptree <- as.phylo(hclust(dist(m2)))
tmptree$edge.length <- tmptree$edge.length/sum(tmptree$edge.length)
abline(v=dist.topo(flytree, tmptree, method='score'), col='blue', lty=2)

tmptree <- as.phylo(hclust(dist(cbind(f2, m2))))
tmptree$edge.length <- tmptree$edge.length/sum(tmptree$edge.length)
abline(v=dist.topo(flytree, tmptree, method='score'), col='blue', lty=2)




## analyze untargeted data by strain
glob[1:4, 1:8]
aggregate(glob[, 4:ncol(glob)], list(glob$sex, glob$strain), mean, na.rm=T) -> glob.st.means
glob.st.means[1:5, 1:5]
names(glob.st.means)[1:2] <- c('sex', 'strain')

glob.st.means$species <- sapply(strsplit(glob.st.means$strain, split='_'), "[[", 1)
glob.st.means <- glob.st.means[ ,c(4421, 1:4420)]
glob.st.means[1:5, 1:5]

table(colSums(is.na(glob.st.means))[4:ncol(glob.st.means)]) # 361 features present at least once in each strain AND in BOTH sexes
table(colSums(is.na(glob.st.means[glob.st.means$sex == 'females', ]))[4:ncol(glob.st.means)]) # 553 features present at least once in each strain in females
table(colSums(is.na(glob.st.means[glob.st.means$sex == 'males', ]))[4:ncol(glob.st.means)]) # 608 features present at least once in each strain in males

glob.st.means[1:5, 1:5]
glob.st.means[glob.st.means$sex == 'males', ] -> m
m[1:4,1:4]
m <- m[ ,colSums(is.na(m)) ==0]
m <- aggregate(m[, 4:ncol(m)], list(m$species), mean, na.rm=T)
rownames(m) <- m$Group.1
m <- as.matrix(m[ ,-1])


glob.st.means[glob.st.means$sex == 'females', ] -> f
f[1:4,1:4]
f <- f[ ,colSums(is.na(f)) ==0]
f <- aggregate(f[, 4:ncol(f)], list(f$species), mean, na.rm=T)
rownames(f) <- f$Group.1
f <- as.matrix(f[ ,-1])



# a simple tree plot:
ftree <- hclust(dist(f), method='ave')
plot(as.phylo(ftree), cex = 1, label.offset = 0.5)

## bootstrapping the tree
bs <- boot.phylo(ftree, f, FUN = function(xx) as.phylo(hclust(dist(xx))), B=1000)
tmptree$node.label <- bs

par(mfrow=c(1,1))
par(mar=c(3,3,3,3))
par(xpd=TRUE)
plot(tmptree, main='untargeted, female n=553 features', label.offset=0.5)
nodelabels(pie=(tmptree$node.label/1000), cex=1, piecol=c('grey',0))


## measure the distance between trees:
tmptree <- as.phylo(hclust(dist(f)), method='ave') # specify data here (e.g. male or female data)
tmptree$edge.length <- tmptree$edge.length/sum(tmptree$edge.length)

dist.topo(flytree, tmptree, method='score')
sum(dist.topo(flytree, tmptree, method='score') > distance) / 10000  # empirical P value for tree-tree distance