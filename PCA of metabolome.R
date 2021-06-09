# PCA of the multi-species metabolome:

library(DMwR)
library(vegan)
library(splitstackshape)
library(phytools)
library(plyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(RColorBrewer)
library(AssocTests)
library(dendextend)

mycolors = c("#FF0000", '#FF6347', "#AFEEEE", "#87CEFA")
mycol4 = c('#F4A460', "#CD6600", "#87CEFA", "#104E8B")

setwd("/Volumes/GoogleDrive/My Drive/Documents/Hoffman_Metabolite_Phylogeny_Paper/Harrison_Hoffman_GitHub") 

flytree <- read.tree('fly.species.list.nwk')
drop.tip(flytree, 'persimilis') ->flytree_no.persimilis # the day 31 data has no D.persimilis (n=10 species)

read.table('comparative.mz.targeted.data.log_only', header=T, sep='\t') -> targ
head(targ) # a leading X is added to the metabolite names that begin with a number

names(targ)[1:18]
targ[, 18:ncol(targ)] -> mzdat
mzdat[1:4,1:4]
mets <- as.matrix(mzdat)
mets <- mets[ ,colSums(is.na(mets)) == 0] # rid matrix of mzs with missing data
centered_mets <- t(scale(t(mets), scale=F)) # center, but don't scale metabolites by sample
cbind(targ[ ,1:17], centered_mets) -> centered_dat

# PCA of targeted data
##########################
names(centered_dat)[1:19]
targ.pca <- prcomp(centered_dat[ ,18:ncol(centered_dat)], scale=F)
plot(targ.pca)
summary(targ.pca)
eigs <- targ.pca$sdev^2

barplot(eigs/sum(eigs), space=0, border=NA, xlab='Principal Components', las=1, col=c(rep('grey20', 15), rep(8, length(targ.pca$sdev))), ylab='% variance explained')

pairs(targ.pca$x[ ,1:4], pch=16, col=as.numeric(as.factor(targ$species)))
pairs(targ.pca$x[ ,1:4], pch=16, col=as.numeric(as.factor(targ$sex)))
pairs(targ.pca$x[ ,1:4], pch=16, col=mycol4[as.numeric(as.factor(targ$age))*2])

pairs(targ.pca$x[ ,1:4], pch=16, col=as.numeric(as.factor(targ$sex)), upper.panel = NULL, diag.panel = NULL)
pairs(targ.pca$x[ ,1:4], pch=16, col=as.numeric(as.factor(targ$age))+2, lower.panel = NULL, diag.panel = NULL)

targ.pca <- as.data.frame(targ.pca$x[,1:6])
targ.pca <- merge(targ[ ,c(2,4,5,8)], targ.pca, by.x="row.names", by.y="row.names")
head(targ.pca)
targ.pca$species <- as.factor(targ.pca$species)
targ.pca$strain <- as.factor(targ.pca$strain)
targ.pca$age <- as.factor(targ.pca$age.days)
targ.pca$sex <- as.factor(targ.pca$sex)

str(targ.pca)
targ.pca <- targ.pca[ ,c(2:4,12,6:11)]

age.centroids <- merge(targ.pca, aggregate(cbind(mean.PC1=PC1, mean.PC2=PC2, mean.PC3=PC3, mean.PC4=PC4) ~ age, targ.pca, mean), by='age') # calculate centroids for PCs 1 to 4


ggplot(age.centroids, aes(PC1 , PC2, color=factor(age))) +
  stat_ellipse(aes(x=PC1,y=PC2, fill=as.factor(age)), geom="polygon", level=0.5, alpha=0.1) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_point(size=1) +
  geom_point(aes(x=mean.PC1, mean.PC2), size=2) +
  geom_segment(aes(x=mean.PC1, mean.PC2, xend=PC1, yend=PC2))


as.factor(paste (targ.pca$sex, targ.pca$age)) -> targ.pca$sex_age
levels(targ.pca$sex_age) 
levels(targ.pca$sex_age) <- c('old female', 'young female', 'old male', 'young male')

sex_age.centroids <- merge(targ.pca, aggregate(cbind(mean.PC1=PC1, mean.PC2=PC2, mean.PC3=PC3, mean.PC4=PC4) ~ sex_age, targ.pca, mean), by='sex_age') # calculate centroids for PCs 1 to 4

mean.PCS <- aggregate(cbind(mean.PC1=PC1, mean.PC2=PC2, mean.PC3=PC3, mean.PC4=PC4) ~ sex_age, targ.pca, mean) 
sex_age.centroids <- merge(targ.pca, mean.PCS, by= 'sex_age')

ggplot(sex_age.centroids, aes(PC1 , PC2, color=factor(sex_age))) +
  scale_color_manual(values=c(mycol4)) +
  xlab("PC1 (17.7%)")+
  ylab("PC2 (12.4%)")+
  stat_ellipse(aes(x=PC1,y=PC2, fill=as.factor(sex)), geom="polygon", level=0.5, alpha=0.2, linetype = 0) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_point(size=1) 

## a simple plot to steal the legend from to pair with the previous plot
ggplot(sex_age.centroids, aes(PC1 , PC2, color=factor(sex_age))) +
  stat_ellipse(aes(x=PC1,y=PC2, fill=as.factor(sex_age)), geom="polygon", level=0.5, alpha=0.2, linetype = 0) +
  scale_fill_manual(values=c(mycol4)) +
  scale_color_manual(values=c(mycol4)) +
  xlab("PC1 (17.7%)")+
  ylab("PC2 (12.4%)")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_point(size=1) +
  geom_point(aes(mean.PC1, mean.PC2), size=0) +
  geom_segment(aes(mean.PC1, mean.PC2, xend=PC1, yend=PC2), size=0.2)

# plot PCs on phylogeny:
as.factor(paste (targ.pca$species, targ.pca$sex, targ.pca$age)) -> targ.pca$species_sex_age
levels(targ.pca$species_sex_age) 

targ.pca[1:4,1:12]
ssa.means <- aggregate(targ.pca[ ,c(5:10)], by=list(targ.pca$species_sex_age), mean) # calculate centroids for PCs 1 to 4

head(ssa.means)
ssa.means <- cSplit(ssa.means, 'Group.1', sep=' ')
colnames(ssa.means)[7:9] <- c('species', 'sex', 'age')
ssa.means$age <- as.factor(ssa.means$age)
str(ssa.means)
head(ssa.means)

x <- as.data.frame(print(pivot_wider(ssa.means[ ,c(1:2, 7:9)], names_from= c(sex, age), values_from=c(PC1, PC2))))
xx <- print(as.matrix(x[ ,2:ncol(x)]))
rownames(xx) <- x$species
xx <- xx[flytree$tip.label, ]

plotTree.barplot(flytree, xx[ ,1:4])


plotTree.barplot(flytree, xx[ ,1:4], args.barplot=list(beside=TRUE, col=mycol4, space=c(0,2), legend.text=F)) # plot PC1
mtext("PC1 (17.5%)",1,at=0.5,line=2.5)
plotTree.barplot(flytree, xx[ ,5:8], args.barplot=list(beside=TRUE, col=mycol4, space=c(0,2), legend.text=F)) # plot PC2
mtext("PC2 (11.3%)",1,at=0.5,line=2.5)



################################################
# get phylosignal (lambda and K) values for PCs
# for each age and sex combination, need species-level means, and the se of the species means from the strain values
targ.pca[1:5,1:10]

# persimilis has no day 31 data, create a tree to use with the day 31 data:
flytree_no_per <- drop.tip(flytree, 'persimilis') 
plot(flytree_no_per)
##

unique(targ.pca$age) -> ages
as.character(levels(targ.pca$sex)) -> sexes

sem <- function(x) {sd(x)/sqrt(length(x))}

output <- matrix(nrow=4, ncol=4)
unique(paste0(targ.pca$sex, '_', targ.pca$age)) -> colnames(output)
rownames(output) <- c('lambda', 'lambda_P', 'K', 'K_P')
output

PC <- 3 # enter PC to test here

for(i in 1:2) {
  sexes[i] -> SEX
  targ.pca[targ.pca$sex == SEX, ] -> tmp
  for(j in 1:2) {
      ages[j] -> AGE
    tmp[tmp$age == AGE, ] -> tmp2
    paste0(SEX, '_', AGE) -> tested.group

    head(tmp2)
means <- aggregate(tmp2[ ,5:10], by=list(tmp2$species), mean)
rownames(means) <- means$Group.1
as.matrix(means[ ,2:ncol(means)]) -> means
ifelse(AGE == '5', means[flytree$tip.label, ] -> means, means[flytree_no_per$tip.label, ] -> means) # aligns phenotype and phylogeny

sems <- aggregate(tmp2[ ,5:10], by=list(tmp2$species), sem)
rownames(sems) <- sems$Group.1
as.matrix(sems[ ,2:ncol(sems)]) -> sems
ifelse(AGE == '5', sems[flytree$tip.label, ] -> sems, sems[flytree_no_per$tip.label, ] -> sems) # aligns phenotype and phylogeny

# fill in NA (from a species with a single strain) with max se
apply(sems, 2, max, na.rm=T) -> max.se
rownames(sems)[is.na(sems[ ,1])] -> species.missing.error
for(k in 1:ncol(sems)) {
colnames(sems)[k] -> name
sems[species.missing.error, name] <- max.se[name]} 
#

# PC is specified in these lines: means[ ,_]
ifelse(AGE == '5', phylosig(flytree, means[ ,PC], method='lambda', test=T, se=sems[ ,3]) -> l, phylosig(flytree_no_per, means[ ,PC], method='lambda', test=T, se=sems[ ,3]) -> l)

ifelse(AGE == '5', phylosig(flytree, means[ ,PC], method='K', test=T, se=sems[ ,3]) -> k, phylosig( flytree_no_per, means[ ,PC], method='K', test=T, se=sems[ ,3]) -> k)

c(l$lambda, l$P, k$K, k$P) -> out
out -> output[ ,tested.group] }}

t(output)



#########################################################################################
# PCA of the untargeted metabolome
#########################################################################################
read.table('Comparative metabolomics.untargeted pos mode logged-only.csv', sep=',' , stringsAsFactors = F, header=T) -> pos

read.table('Comparative metabolomics.untargeted neg mode logged-only.csv', sep=',', header = T, stringsAsFactors = F) -> neg  


### combining pos and neg data:
pos[1:4,1:10]
neg[1:4, 1:10]
table(neg$Met.sample.number %in% pos$Met.sample.number)
pos[ ,c(1, 9:ncol(pos))] -> pos.data
merge(neg, pos.data, by='Met.sample.number') -> glob
glob[1:4,1:10]
colnames(glob)[6] <- 'sex'
table(glob$ID)
paste(glob$species, glob$strain, sep='_') -> glob$strain
str(glob)

table(glob$strain, glob$sex) # sample size per sex per strain
## pca of individual samples (will do strain-level means later, which increases the number of non-missing features)

globmat <- as.matrix(glob[ ,9:ncol(glob)])

# impute features
hist(colSums(is.na(globmat)), las=1, 40)

plot(colSums(is.na(globmat))[order(colSums(is.na(globmat)))])

impute.these <- globmat[ ,colSums(is.na(globmat)) < 40]
dim(impute.these)
imputed <- bnstruct::knn.impute(impute.these, k=10)
glob.imputed <- cbind(glob[ ,4:6], imputed)
glob.imputed[1:4,1:6]

pca <- prcomp(glob.imputed[ ,4:ncol(glob.imputed)], scale=T)

# tracy widom test to determine the number of significant PCs:
par(mfrow=c(1,2))
plot(pca$sdev^2, las=1, ylab='eigenvalue', xlab='principal component', pch=16)
ev <- pca$sdev^2
tw(ev, length(ev), criticalpoint = 0.9793) -> h
sigPCs <- h$SigntEigenL # number of significant PCs
plot(h$statistic, las=1, pch=16, ylab='Tracy Widom statistic', xlab='principal component')
abline(h=0)

pca_all <- as.data.frame(pca$x[ ,1:sigPCs])
pca_all <- cbind(glob[ ,4:6], pca_all)

head(pca_all)
pca_all$species <- as.factor(pca_all$species)
pca_all$strain <- as.factor(pca_all$strain)
pca_all$sex <- as.factor(pca_all$sex )

ggplot(pca_all, aes(x=PC1,y=PC2))+
  geom_point(aes(color=species, size=3, shape=sex))+ 
  scale_shape_manual(values=c(1, 19))+
  stat_ellipse(aes(x=PC1,y=PC2, fill=sex), geom="polygon", level=0.6, alpha=0.1)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

all_centroids <- merge(pca_all, aggregate(cbind(mean.PC1=PC1, mean.PC2=PC2, mean.PC3=PC3, mean.PC4=PC4) ~ species + sex, pca_all, mean), by='species') # calculate centeroids for PCs 1 to 4

all_centroids[1:8, ]
all_centroids <- all_centroids[all_centroids$sex.x == all_centroids$sex.y, ]
all_centroids <- all_centroids[ ,colnames(all_centroids) != 'sex.y']
colnames(all_centroids)[3] <- 'sex'
all_centroids[1:8, ]

ggplot(all_centroids, aes(PC1 , PC2, color=species)) +
  xlab("PC1 (___%)")+
  ylab("PC2 (___%)")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_point(size=4, aes(shape=sex)) +
  scale_shape_manual(values=c(1, 19))+
  geom_point(aes(mean.PC1, mean.PC2), size=0) +
  geom_segment(aes(mean.PC1, mean.PC2, xend=PC1, yend=PC2), size=0.2)


##############################
# PCA on each sex separately:
##############################
glob.imputed[1:4,1:6]
fpca <- prcomp(glob.imputed[glob.imputed$sex == 'females', 4:ncol(glob.imputed)])

# tracy widom test to determine the number of significant PCs:
par(mfrow=c(1,2))
plot(fpca$sdev^2, las=1, ylab='eigenvalue', xlab='principal component', pch=16)
ev <- fpca$sdev^2
tw(ev, length(ev), criticalpoint = 0.9793) -> h
sigPCs <- h$SigntEigenL # number of significant PCs
plot(h$statistic, las=1, pch=16, ylab='Tracy Widom statistic', xlab='principal component')
abline(h=0)

fpca <- as.data.frame(fpca$x[ ,1:sigPCs])
fpca <- cbind(glob[glob$sex == 'females',4:6], fpca)

head(fpca)
fpca$species <- as.factor(fpca$species)
fpca$strain <- as.factor(fpca$strain)
fpca$sex <- as.factor(fpca$sex )

ggplot(fpca, aes(x=PC1,y=PC2))+
  geom_point(aes(color=species, size=3))+ 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

ggplot(fpca, aes(x=PC1,y=PC2))+
  geom_point(aes(color=species))+ 
  stat_ellipse(aes(x=PC1,y=PC2, fill=species), geom="polygon", level=0.6, alpha=0.1)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


# mean of PC by species
fcents <- merge(fpca[ ,c(1,3)], aggregate(cbind(mean.PC1=PC1, 
                                                  mean.PC2=PC2, 
                                                  mean.PC3=PC3, 
                                                  mean.PC4=PC4,
                                                  mean.PC5=PC5,
                                                  mean.PC6=PC6,
                                                  mean.PC7=PC7,
                                                  mean.PC8=PC8,
                                                  mean.PC9=PC9) ~ species, fpca, mean), by='species') 

fcents <- fcents[!duplicated(fcents), ]
head(fcents)

# neighbor-joining tree of distance between species-mean PCs 1-N, where N is the largest significant PC (see tracy-widom step):
fmat <- as.matrix(fcents[ ,3:ncol(fcents)])
rownames(fmat) <- fcents$species
df <- dist(fmat)
fpctree <- nj(df) # 

table(glob$species, glob$sex)
table(glob$species, glob$sex)[ ,1] # sample size for females
table(glob$species, glob$sex)[ ,2] # sample size for males

fpctree$tip.label <- paste0(fpctree$tip.label, ' (', table(glob$species, glob$sex)[ ,1], ')') # add sample size info to tree labels

# inspect tree, perform manual node rotation and rooting to compare to flytree
par(mfrow=c(1,3)) 
plot(flytree, cex=1.5)
plot(fpctree, cex=1.5)
nodelabels()
plot(reroot(fpctree, 18), cex=1.5)
plot(rotateNodes(reroot(fpctree, 18), c(12, 19)), cex=1.5)
plot(rotateNodes(rotateNodes(reroot(fpctree, 18), c(19)), 'all'), cex=1.5)

fpctree$edge.length <- fpctree$edge.length/sum(fpctree$edge.length)
dist.topo(flytree, fpctree, 'score')




# PCA for male data:
mpca <- prcomp(glob.imputed[glob.imputed$sex == 'males', 4:ncol(glob.imputed)])

# tracy widom test to determine the number of significant PCs:
par(mfrow=c(1,2))
plot(mpca$sdev^2, las=1, ylab='eigenvalue', xlab='principal component', pch=16)
ev <- mpca$sdev^2
tw(ev, length(ev), criticalpoint = 0.9793) -> h
sigPCs <- h$SigntEigenL # number of significant PCs
plot(h$statistic, las=1, pch=16, ylab='Tracy Widom statistic', xlab='principal component')
abline(h=0)

mpca <- as.data.frame(mpca$x[ ,1:sigPCs])
mpca <- cbind(glob[glob$sex == 'males',4:6], mpca)

head(mpca)
mpca$species <- as.factor(mpca$species)
mpca$strain <- as.factor(mpca$strain)
mpca$sex <- as.factor(mpca$sex )

ggplot(mpca, aes(x=PC1,y=PC2))+
  geom_point(aes(color=species))+ 
  stat_ellipse(aes(x=PC1,y=PC2, fill=species), geom="polygon", level=0.6, alpha=0.1)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


mcents <- merge(mpca[ ,c(1,3)], aggregate(cbind(mean.PC1=PC1, 
                                                mean.PC2=PC2, 
                                                mean.PC3=PC3, 
                                                mean.PC4=PC4,
                                                mean.PC5=PC5,
                                                mean.PC6=PC6,
                                                mean.PC7=PC7,
                                                mean.PC8=PC8,
                                                mean.PC9=PC9) ~ species, mpca, mean), by='species') 

mcents <- mcents[!duplicated(mcents), ]
head(mcents)

# neighbor-joining tree of distance between species-mean PCs 1-N, where N is the largest significant PC (see tracy-widom step):
fmat <- as.matrix(mcents[ ,3:ncol(mcents)])
rownames(fmat) <- mcents$species
dm <- dist(fmat)
mpctree <- nj(dm) # 

table(glob$species, glob$sex)
table(glob$species, glob$sex)[ ,1] # sample size for females
table(glob$species, glob$sex)[ ,2] # sample size for males

mpctree$tip.label <- paste0(mpctree$tip.label, ' (', table(glob$species, glob$sex)[ ,2], ')') # add sample size info to tree labels

# inspect tree, perform manual node rotation and rooting to compare to flytree
par(mfrow=c(1,3)) 
plot(flytree, cex=1.5)
plot(mpctree, cex=1.5)
nodelabels()
plot(rotateNodes(rotateNodes(reroot(mpctree, 15), c(19)), 'all'), cex=1.5)

mpctree$edge.length <- mpctree$edge.length/sum(mpctree$edge.length)
dist.topo(flytree, mpctree, 'score')


####################################################################
## strain-level data
####################################################################

## get means for each strain.sex:
aggregate(glob[, 9:ncol(glob)], list(glob$sex, glob$strain), mean, na.rm=T) -> glob.strain.means
glob.strain.means[1:9, 1:5]
names(glob.strain.means)[1:2] <- c('sex', 'strain')
separate(glob.strain.means, 'strain', into=c('species', 'strain')) -> glob.strain.means
glob.strain.means[1:4, 1:10]
glob.strain.means$strain <- paste(glob.strain.means$species, glob.strain.means$strain, sep='_')
glob.strain.means[1:9, 1:5]

table(colSums(is.na(glob.strain.means[ ,4:ncol(glob.strain.means)]))==0)
table(colSums(is.na(glob.strain.means[ ,4:ncol(glob.strain.means)]))==1)

# impute metabolites with only a single missing sample
impute.these <- glob.strain.means[ ,colSums(is.na(glob.strain.means[ ,4:ncol(glob.strain.means)]))<2]
impute.these[1:4,1:4]
impute.these <- as.matrix(impute.these[ ,4:ncol(impute.these)])
imputed <- bnstruct::knn.impute(impute.these, k=10)
table(colSums(is.na(imputed)))    
glob.strain.means <- cbind(glob.strain.means[ ,1:3], imputed)

pca <- prcomp(glob.strain.means[ ,4:ncol(glob.strain.means)])
summary(pca)

# tracy widom test to determine the number of significant PCs:
par(mfrow=c(1,2))
plot(pca$sdev^2, las=1, ylab='eigenvalue', xlab='principal component', pch=16)
ev <- pca$sdev^2
tw(ev, length(ev), criticalpoint = 0.9793) -> h
sigPCs <- h$SigntEigenL # number of significant PCs
plot(h$statistic, las=1, pch=16, ylab='Tracy Widom statistic', xlab='principal component')
abline(h=0)

glob.pca <- as.data.frame(pca$x[ ,1:sigPCs])
glob.pca <- cbind(glob.strain.means[ ,1:3], glob.pca)

head(glob.pca)

centroids <- merge(glob.pca, aggregate(cbind(mean.PC1=PC1, mean.PC2=PC2, mean.PC3=PC3, mean.PC4=PC4) ~ species + sex, glob.pca, mean), by='species') # calculate centeroids for PCs 1 to 4

centroids[1:8, ]
centroids <- centroids[centroids$sex.x == centroids$sex.y, ]
centroids <- centroids[ ,colnames(centroids) != 'sex.y']
colnames(centroids)[colnames(centroids) == 'sex.x'] <- 'sex'

ggplot(centroids, aes(x=PC1,y=PC2))+
  geom_point(aes(color=species))+ 
  stat_ellipse(aes(x=PC1,y=PC2, fill=sex), geom="polygon", level=0.6, alpha=0.1) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

ggplot(centroids, aes(PC1 , PC2, color=species)) +
  xlab("PC1 (17.7%)")+
  ylab("PC2 (12.4%)")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_point(size=1) +
  geom_point(aes(mean.PC1, mean.PC2), size=0) +
  geom_segment(aes(mean.PC1, mean.PC2, xend=PC1, yend=PC2), size=0.2)




