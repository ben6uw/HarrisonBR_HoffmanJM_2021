# Analyze distance between samples over species divergence time

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.11")
BiocManager::install("limma")

library(plyr) # data wrangling
library(splitstackshape)
library(dplyr)
library(tidyr)
library(reshape2)
library(readr)
library(DMwR) 
library(impute)

library(sjstats) # modeling tools
library(lme4)
library(lmerTest)
library(nlme) 
library(MuMIn) 
library(MASS)
library(smatr)
library(lmodel2)
library(mvtnorm)
library(nlreg)
library(lsmeans)
library(igraph)
library(AssocTests)

library(graphics) # plotting
library(ggplot2) 
library(gplots) 
library(corrplot)
library(RColorBrewer)
library(limma)
library(psych)
require(gridExtra)

library(ape) # phylogenetic analysis
library(phytools) 
library(geiger) 
library(caper) 

gc()
dev.off()

mycol7 <- brewer.pal(7, 'RdYlBu')
mycol4 = c("#EEC900", "#CD6600", "#87CEFA", "#104E8B")

setwd("/Volumes/GoogleDrive/My Drive/Documents/Hoffman_Metabolite_Phylogeny_Paper/Harrison_Hoffman_GitHub") 

flytree <- read.tree('fly.species.list.nwk')
plot(flytree)
axisPhylo(side = 1, root.time = NULL, backward = TRUE)

dat <- read.table('comparative.mz.targeted.data.log_only', header=T, sep='\t')
head(dat) # a leading X is added to the metabolite names that begin with a number
table(dat$species) * 97

names(dat)[1:18]
dat[, 18:ncol(dat)] -> mets
mets[1:4,1:4]
mets <- as.matrix(mets)
mets <- mets[ ,colSums(is.na(mets)) == 0] # rid matrix of mzs with missing data
# remove metabolites identified by the good doctor Danijel Djukovic as artifacts
ghost_peaks <- c('4-Trimethylammoniobutanoate', 'X4.Trimethylammoniobutanoate', 'Sorbitol', 'Xanthosine', 'Oxypurinol', 'Uridine', 'Trimethylamine', 'Valine', 'Adenine', 'Pipecolic acid', 'Pipecolic.acid') # these were one member of a few pair of highly-correlated metabolties (r>0.98 in another analysis),Danijel suspected that these metabolties are actualy LC-MS artifacts
mets <- mets[ ,!colnames(mets) %in% ghost_peaks] 
centered_mets <- t(scale(t(mets), scale=F)) # center, but don't scale metabolites by sample

# a look at the data
par(mfrow=c(2,2))
boxplot(mets[order(apply(mets, 1, median)), ], use.cols = F, main='log only', xlab='sample', las=1)
boxplot(centered_mets[order(apply(mets, 1, median)), ], use.cols = F, main='log, mean-centered', xlab='sample', las=1)
boxplot(mets[ ,order(colMeans(mets))], use.cols = T, main='log only', xlab='metabolite', las=1)
boxplot(centered_mets[ ,order(colMeans(mets))], use.cols = T, main='log, mean-centered', xlab='metabolite', las=1)

dat <- cbind(dat[ ,1:17], centered_mets)
dat[1:5, 1:20]

# sample parameters and other useful indexes
dat$species -> Species 
dat$sex -> Sex
dat$strain -> Strain
dat$age.days -> Age

setLabels = c("young female", "old female", 'young male', 'old male')
shortLabels = c('F 5', 'F 31' ,'M 5', 'M 31')
nSets=4

#################################################
## metabolome distances between samples:
#################################################

# the (mean of) squared difference in log-metabolite levels between samples. 
# (logA - logB)^2 is the 'conventional' measure of trait divergence (Bedford and Hartl 2009, Lande 1976; Felsenstein 2004, Khaitovich et al., 2004)
rownames(centered_mets) <- paste(Species, Strain, Sex, Age, 1:nrow(centered_mets))
centered_mets[1:10,1:6]

unique(Age) -> ages
unique(Sex) -> sexes

for(i in 1:length(ages)){
  for(j in 1:2) {
    ages[i] -> day
    sexes[j] -> sex.type
    
    mat <- centered_mets[Age == day & Sex == sex.type, ] 
     p <- outer(1:nrow(mat),1:nrow(mat), FUN = Vectorize( function(i,j) mean((mat[i,] - mat[j,])^2 )) ) # the mean squared difference in each log(mz)
    p[lower.tri(p ,diag=TRUE)] <- 0
  rownames(p) <- rownames(mat)
  colnames(p) <- rownames(mat)
  d <- reshape2::melt(p, varnames = c("sample_A", "sample_B"))   
   d <- d[d$value != 0, ]  
   d$sex <- sex.type
    d$age <- day
    head(d)
    cSplit(d, 'sample_A', sep=' ') -> d
    cSplit(d, 'sample_B', sep=' ') -> d
    as.character(d$`sample_A_1`) -> d$`sample_A_1`
    as.character(d$`sample_B_1`) -> d$`sample_B_1`
    as.factor(ifelse(d$`sample_A_1` == d$`sample_B_1`, 'within_species', 'between_species')) -> d$species.comparison

    if(i == 1 & j == 1){ d -> Euc.dist  }  
    if(i > 1 | j > 1) {rbind(Euc.dist, d) -> Euc.dist}
  }}

head(Euc.dist) 
Euc.dist$age <- as.factor(Euc.dist$age)
Euc.dist$sex <- as.factor(Euc.dist$sex)
table(Euc.dist$species.comparison)
table(Euc.dist$age, Euc.dist$species.comparison)
table(Euc.dist$sex, Euc.dist$species.comparison)
names(Euc.dist)[1] <- 'distance'


xyplot(distance ~ age | species.comparison + sex,  data = Euc.dist, type=c('p', 'r'), xlab = "age (days)", ylab = "metabolome distance")
xyplot(distance ~ species.comparison | age + sex,  data = Euc.dist, type=c('p', 'r'), xlab = "age (days)", ylab = "metabolome distance")

# relate metabolome distance to divergence times:
cophenetic.phylo(flytree)/2 -> tip_tips # pairwise divergence times

divergence_times <- numeric()
for(i in 1:nrow(Euc.dist)) {
  tip_tips[Euc.dist$'sample_A_1'[i], Euc.dist$'sample_B_1'[i]] -> divergence_times[i]
}
divergence_times -> Euc.dist$divergence_time
Euc.dist$sex <- as.factor(Euc.dist$sex)

xyplot(distance ~ divergence_time | age + sex,  data = Euc.dist, type=c('p'), xlab = "divergence time (MY)", ylab = "metabolome distance")
xyplot(distance ~ divergence_time | age + sex, data = Euc.dist, type=c('p', 'r'), xlab = "divergence time (MY)", ylab = "metabolome distance")

par(mfrow=c(2,2))
par(mar=c(5,4,4,2))
plot(lm(distance ~ divergence_time * age * sex, data=Euc.dist)) 
summary(lm(distance ~ divergence_time * age * sex, data=Euc.dist)) # r-squared of 0.38

age.x.sex.x.mod <- lm(distance ~ divergence_time * age * sex, data=Euc.dist)
age.sex.x.mod <- lm(distance ~ divergence_time + age * sex, data=Euc.dist)
age.x.sex.mod <- lm(distance ~ divergence_time * age + sex, data=Euc.dist)
age.sex.mod <- lm(distance ~ divergence_time + age + sex, data=Euc.dist)
age.mod <- lm(distance ~ divergence_time + age, data=Euc.dist)
age.x.mod <- lm(distance ~ divergence_time * age, data=Euc.dist)
sex.mod <- lm(distance ~ divergence_time + sex, data=Euc.dist)
sex.x.mod <- lm(distance ~ divergence_time * sex, data=Euc.dist)
mod <- lm(distance ~ divergence_time, data=Euc.dist) 

anova(age.x.sex.x.mod , age.x.sex.mod, age.sex.x.mod, age.sex.mod, age.x.mod, age.mod, sex.x.mod, sex.mod, mod) # model comparison
round(anova(age.x.sex.x.mod), 5)# anova table for full model

summary(age.x.sex.x.mod)
anova(age.x.sex.x.mod)
trends <- lstrends(age.x.sex.x.mod, c("age", "sex"), var="divergence_time")
trends
pairs(trends)
lstrends(age.x.sex.x.mod, c("age"), var="divergence_time")

summary(aov(distance ~ divergence_time * age * sex, Euc.dist)) # the slope of distance ~ divergence is different between the ages, but does not differ by sex.  
pf(184.3584, 1, 1057, lower.tail = F) # calculate exact P given the F and degrees of freedom (above)

# plot the regression within each group:
Euc.dist$sex_age <- (as.factor(paste(Euc.dist$sex, Euc.dist$age)))
str(Euc.dist)
levels(Euc.dist$sex_age)

qplot(x = divergence_time, y = distance, data = Euc.dist, facets = ~ sex, colour=sex_age) +
  geom_smooth(method = "lm") +
  scale_color_manual(name='sex_age', values = mycol4) +
  theme_minimal() + 
  theme(axis.text.x = element_blank())+
  labs(x="divergence time (million years)", y = "metabolome distance", size=4) +
  scale_y_continuous(breaks=seq(0, 1.4, 0.2)) +
  theme_classic() +
  theme(panel.spacing = unit(2, "lines")) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.y = element_text(size=14)) +
  theme(axis.title.x = element_text(size=14)) +
    theme(axis.text.x = element_text(size=14))

dev.off()

# within-species variance: do our data support the theory of metabolic dysregulation?
head(Euc.dist)
summary(lm(distance ~ age * sex, Euc.dist[Euc.dist$divergence_time == 0, ]))


# Approach 2:
#################################################
# the change (reduction) in correlation between samples for all metabolites.  This does have the issue of co-linearity
# this measure has been used to look at the evolution of the transcriptome:  Carvunis et al., 2015 eLife, Coolon et al., Genome Research 2014

rownames(centered_mets) <- paste(Species, Strain, Sex, Age, 1:nrow(centered_mets))
centered_mets[1:10,1:6]

unique(Age) -> ages
unique(Sex) -> sexes

for(i in 1:length(ages)){
  for(j in 1:2) {
    ages[i] -> day
    sexes[j] -> sex.type
    mat <- centered_mets[Age == day & Sex == sex.type, ] 
    1-cor(t(mat)) -> p
    p[lower.tri(p ,diag=TRUE)] <- 0
    d <- reshape2::melt(p, varnames = c("sample_A", "sample_B"))   
    d[d$value != 0, ] -> d
    d$sex <- sex.type
    d$age <- day
    head(d)
    cSplit(d, 'sample_A', sep=' ') -> d
    cSplit(d, 'sample_B', sep=' ') -> d
    as.character(d$`sample_A_1`) -> d$`sample_A_1`
    as.character(d$`sample_B_1`) -> d$`sample_B_1`
    as.factor(ifelse(d$`sample_A_1` == d$`sample_B_1`, 'within_species', 'between_species')) -> d$species.comparison
    
    if(i == 1 & j == 1){ d -> cor.dist  }  
    if(i > 1 | j > 1) {rbind(cor.dist, d) -> cor.dist}
  }}

head(cor.dist) 
cor.dist$age <- as.factor(cor.dist$age)
cor.dist$sex <- as.factor(cor.dist$sex)
table(cor.dist$species.comparison)
table(cor.dist$age, cor.dist$species.comparison)
table(cor.dist$sex, cor.dist$species.comparison)
names(cor.dist)[1] <- 'distance'


xyplot(distance ~ age | species.comparison + sex,  data = cor.dist, type=c('p', 'r'), xlab = "age (days)", ylab = "metabolome distance")
xyplot(distance ~ species.comparison | age + sex,  data = cor.dist, type=c('p', 'r'), xlab = "age (days)", ylab = "metabolome distance")

summary(lm(distance ~ age * sex * species.comparison,  data = cor.dist))
round(anova(lm(distance ~ age * sex * species.comparison,  data = cor.dist)),3)
m1 <- lm(distance ~ age * species.comparison * sex,  data = cor.dist)

# relate metabolome distance to divergence times:
cophenetic.phylo(flytree)/2 -> tip_tips # pairwise divergence times

divergence_times <- numeric()
for(i in 1:nrow(cor.dist)) {
  tip_tips[cor.dist$'sample_A_1'[i], cor.dist$'sample_B_1'[i]] -> divergence_times[i]
}
divergence_times -> cor.dist$divergence_time
cor.dist$sex <- as.factor(cor.dist$sex)

xyplot(distance ~ divergence_time | age + sex,  data = cor.dist, type=c('p'), xlab = "divergence time (MY)", ylab = "metabolome distance")
xyplot(distance ~ divergence_time | age + sex, data = cor.dist, type=c('p', 'r'), xlab = "divergence time (MY)", ylab = "metabolome distance")

par(mfrow=c(2,2))
par(mar=c(5,4,4,2))
plot(lm(distance ~ divergence_time * age * sex, data=cor.dist)) 
summary(lm(distance ~ divergence_time * age * sex, data=cor.dist)) # r-squared of 0.28

lm(distance ~ divergence_time * age * sex, data=cor.dist) -> age.sex.mod
lm(distance ~ divergence_time * age, data=cor.dist) -> age.mod
lm(distance ~ divergence_time * sex, data=cor.dist) -> sex.mod
lm(distance ~ divergence_time, data=cor.dist) -> mod

anova(mod, sex.mod, age.mod, age.sex.mod) # model comparison
round(anova(age.sex.mod), 6)# anova table for full model

cor.dist$sex_age <- (as.factor(paste(cor.dist$sex, cor.dist$age)))
str(cor.dist)
levels(cor.dist$sex_age)

qplot(x = divergence_time, y = distance, data = cor.dist, facets = ~ sex, colour=sex_age) +
  geom_smooth(method = "lm") +
  scale_color_manual(name='sex_age', values = mycol4) +
  theme_minimal() + 
  theme(axis.text.x = element_blank())+
  labs(x="divergence time (million years)", y = "metabolome distance", size=4) +
  scale_y_continuous(breaks=seq(0, 1.4, 0.2)) +
  theme_classic() +
  theme(panel.spacing = unit(2, "lines")) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.y = element_text(size=14)) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.text.x = element_text(size=14))

dev.off()

# the results of these two approaches are very similar, both agree regarding the linearity divergence and the effect of age on divergence


######################################################################################
## comparing BM and OU models
# ML estimation of BM and OU parameters, with model comparison by AIC
######################################################################################
# specify models:
OU <- function(x) dft / (2 * sel)  * (1 - exp(-2 * sel * x))  
BM <- function(x) (slope * x) + intercept 

######################################################################################
# Summarize metabolome-wide divergence
# plot fits
# plot residuals
# plot fitted values over time
######################################################################################
## compare OU and BM but omit the intraspecies data (divergence_time==0).  The OU model as-is does not fit these data, rather it is forced to the origin (intraspecies distance == 0), so to give it a chance to outperform BM:

# ML estimation of BM and OU parameters, with model comparison by AIC
par(mfrow=c(2,4))

modCompOut <- matrix(nrow=4, ncol=3)
distance_data <- Euc.dist[Euc.dist$divergence_time !=0, ] # time=0 data removed here
groups <- unique(paste(distance_data$sex, distance_data$age))
  
  for(i in 1:length(groups)){
    group <- groups[i]
    tmp <- distance_data[paste(distance_data$sex, distance_data$age) == group, ] 
    
    OU.fit <- nlreg(distance ~ (drift / (2*selection)) * (1-exp(-2*selection*divergence_time)), start = c(drift = 2, selection = 0.1), data=tmp)
    
    dft <- coef(OU.fit)[1]
    sel <- coef(OU.fit)[2]
    
    BM.fit<- nlreg(distance ~ (slope*divergence_time) + intercept, start = c(slope = 0.006, intercept= 2.4), data=tmp)
    coef(BM.fit)[1] -> slope
    coef(BM.fit)[2] -> intercept
    
    plot(distance ~ divergence_time, tmp, las=1, pch=19, col=rgb(0,0,0,0.5), main=setLabels[i], xlab='time (myr)', ylab='metabolome distance')
    curve(OU, 1, 50, las=1, add=T, lty=2)
    curve(BM, 1, 50, las=1, add=T)
    
    modCompOut[i, 1] <- round(cor(tmp$distance, predict(BM.fit))^2, 4) # r-squared
    modCompOut[i, 2] <- round(cor(tmp$distance, predict(OU.fit))^2, 4) # r-squared
    
    plot(density(OU.fit$residuals), ylim=c(0, 0.5), lty=2, main='', xlab='residuals', las=1)
    lines(density(BM.fit$residuals))
    
    # compare models by AIC:
    f <- AIC(BM.fit, OU.fit)
    modCompOut[i, 3] <- round(f$AIC[1] - f$AIC[2], 3) 
  }
rm(f)

# model comparison table
rownames(modCompOut) <- setLabels
colnames(modCompOut) <- c('BM r-squared', 'OU r-squared', 'delta AIC')
modCompOut <- as.data.frame(modCompOut)
modCompOut$group <- rownames(modCompOut)
modCompOut <- modCompOut[ ,c(4,1:3)]
modCompOut

# make a table of these results, export as pdf as well:
ggsave("FigS2b_metabolome_divergence_AIC_table.pdf", tableGrob(modCompOut, rows=NULL), width=14, height=10) # export a table as pdf

modCompOut
######################################################################################



##############################################################################################################
# analyze the evolutionary models by WGCNA module (this will require you to first run the analysis file: "metabolome PIC network.R" to generate the module identifications)
##############################################################################################################
# we might expect that the OU model might be favored over BM if this analysis was not done metabolome-wide, but rather within metabolic modules, which might have particular evolutionary constraint:

# load module membership data from WGNCA analysis:
setLabels = c("young female", "old female", 'young male', 'old male')
shortLabels = c('F 5', 'F 31' ,'M 5', 'M 31')
nSets=4

Age # a vector of age classes from (dat, see above)
Sex # a vector of sex classes from (dat, see above)
group <- paste(Sex, Age) # an index of group ids

master.d <- list()

for(set in 1:nSets) {
  load(paste0('correlation among the metabolome/', shortLabels[set], '_networkConstruction.RData')) # load the module data
  moduleColors <- as.factor(moduleColors)
  names(moduleColors) <- colnames(centered_mets)
  table(moduleColors)
  mods.to.consider <- levels(moduleColors)[levels(moduleColors) != 'grey'] # omit the non-modular ('grey') metabolites
  
  for(i in 1:length(mods.to.consider)) {
    tmpod <- colnames(centered_mets)[moduleColors == mods.to.consider[i]]
    mat <- centered_mets[group == shortLabels[set], tmpod] 
    
    p <- outer(1:nrow(mat),1:nrow(mat), FUN = Vectorize( function(q,j) mean((mat[q,] - mat[j,])^2 )) )
    p[lower.tri(p ,diag=TRUE)] <- 0
    rownames(p) <- rownames(mat)
    colnames(p) <- rownames(mat)
    d <- reshape2::melt(p, varnames = c("sample_A", "sample_B"))   
    d <- d[d$value != 0, ]  
    d$module <- levels(moduleColors)[i]
    cSplit(d, 'sample_A', sep=' ') -> d
    cSplit(d, 'sample_B', sep=' ') -> d
    as.character(d$`sample_A_1`) -> d$`sample_A_1`
    as.character(d$`sample_B_1`) -> d$`sample_B_1`
    
    head(d) 
    names(d)[1] <- 'distance'
    
    tip_tips <- cophenetic.phylo(flytree)/2 # pairwise divergence times
    
    divergence_times <- numeric()
    for(k in 1:nrow(d)) {
      tip_tips[d$'sample_A_1'[k], d$'sample_B_1'[k]] -> divergence_times[k] }
    
    d$divergence_time <- divergence_times
    d <- d[ ,c(2,1,13)] 
    ifelse(i==1, master.d[[set]] <- d, master.d[[set]] <- rbind(master.d[[set]] , d)) } }


# ML estimation of BM and OU parameters, with model comparison by AIC 
# don't forget to remove the data at distance=0, the OU model assumes this is zero and so it wont fit

set=2
d <- master.d[[set]]
d <- d[d$divergence_time !=0, ]
d$module <- as.factor(d$module)

head(d)

by(d, d[ ,'module'], function(x) coef(nlreg(distance ~ (slope*divergence_time) + intercept, start = c(slope = 0.006, intercept= 2.4), data=x))[2])

by(d, d[ ,'module'], function(x) coef(nlreg(distance ~ (drift / (2*selection)) * (1-exp(-2*selection*divergence_time)), start = c(drift = 2, selection = 0.1), data=x))[2])


# remove intraspecies comparisons (divergence time == 0)
for(set in 1:nSets){
  master.d[[set]] <- master.d[[set]][master.d[[set]]$divergence_time !=0, ] }

# functions for BM and OU fits to be made over the list of distance data:
OUprefunc <- function(x) nlreg(distance ~ (drift / (2*selection)) * (1-exp(-2*selection*divergence_time)), start = c(drift = 2, selection = 0.1), data=x)
ouFunction <- function (x) {return(tryCatch(OUprefunc(x), error=function(e) return('error'))) }

ou <- lapply(master.d, function(x) by(x, x[ ,'module'], ouFunction))
ou2 <- lapply(ou, function(x)  x[x != 'error'])
lapply(ou2, function(x)  x != 'error') # these are the OU models for only those modules that didn't fail to fit.  failure to fit with the message 'singular gradient' means that the extra parameter on the OU model was not useful for the model (ML failed to get the model to 'sit down'), and so this could be interpreted as a failure of the OU model and the BM model would be favored (by default).


BMprefunc <- function(x) nlreg(distance ~ (slope*divergence_time) + intercept, start = c(slope = 0.006, intercept= 2.4), data=x)
bmFunction <- function (x) {return(tryCatch(BMprefunc(x), error=function(e) return('error'))) }

bm <- lapply(master.d, function(x) by(x, x[ ,'module'], bmFunction))
bm2 <- lapply(bm, function(x)  x[x != 'error'])
lapply(bm2, function(x)  x != 'error')



ouAICs <- list()
for(set in 1:nSets){
  tmp <- ou2[[set]]
  ouAICs[[set]] <- lapply(tmp, AIC) }

ouAICs

bmAICs <- list()
for(set in 1:nSets){
  tmp <- bm2[[set]]
  bmAICs[[set]] <- lapply(tmp, AIC) }


deltaAIC <- list()
topModel <- list() #list of model choice for each set and module

for(set in 1:nSets){
  b <- unlist(bmAICs[[set]] )
  o <- unlist(ouAICs[[set]] )
  deltaAIC[[set]] <- b[names(b) %in% names(o)] - o 
  
  bm.is.better <- c(names(b)[!names(b) %in% names(o)], names(b)[names(b) %in% names(o)][b[names(b) %in% names(o)] - o <0])
  topModel[[set]] <- ifelse(names(b) %in% bm.is.better, 'BM', 'OU') 
  names(topModel[[set]]) <- names(bmAICs[[set]])  }

topModel
deltaAIC

plyr::ldply(topModel, rbind)

########## we find evidence for module-level selection (the OU model is preferred over BM for some modules)


# make a useful summary table
# I would like to have columns for 'group', 'module', 'r2 BM', 'r2 OU', and 'delta AIC'

for(set in 1:nSets){
x <- rep(setLabels[set] , length(topModel[[set]]))
ifelse(set==1, grp <- x, grp <- c(grp, x)) }

model.table <- data.frame('group'=grp, 'module'=names(unlist(topModel)), 'model'=unlist(topModel))

unlist(deltaAIC)
deltaAIC

x <-topModel
a <- deltaAIC

for(set in 1:nSets) {
x[[set]][!names(x[[set]])  %in% names(a[[set]]) ] = NA
x[[set]][names(x[[set]])  %in% names(a[[set]]) ] = a[[set]] }
model.table$delta_AIC <- round(as.numeric(unlist(x)), 2)

model.table

write.table(model.table, file= "BM_vs_OU by module table.csv", sep=',', quote=F, row.names=F)

table(model.table$model)

dev.off()

par(mfrow=c(2,3)) 

# plot example(s)
# where OU is preferred:
set=2
d <- master.d[[set]]
module.of.interest <- 'F'
tmp <- d[d$module == module.of.interest ]

OU.fit <- nlreg(distance ~ (drift / (2*selection)) * (1-exp(-2*selection*divergence_time)), start = c(drift = 2, selection = 0.1), data=tmp)
coef(OU.fit)[1] -> dft
coef(OU.fit)[2] -> sel

BM.fit<- nlreg(distance ~ (slope*divergence_time) + intercept, start = c(slope = 0.006, intercept= 2.4), data=tmp)
coef(BM.fit)[1] -> slope
coef(BM.fit)[2] -> intercept

plot(distance ~ divergence_time, tmp, xlab = "divergence time (MY)", ylab = "metabolome distance", las=1, main=paste(setLabels[set], '\nModule', module.of.interest), col=4, pch=16)
curve(OU, 1, 50, las=1, add=T, lty=2)
curve(BM, 1, 50, las=1, add=T)

r2BM <- cor(tmp$distance, predict(BM.fit))^2 # r-squared 
r2OU <- cor(tmp$distance, predict(OU.fit))^2 # r-squared
delAIC <- AIC(BM.fit) - AIC(OU.fit) # delta AIC

legend('topleft', legend= c(paste('r2_BM', round(r2BM,3)), paste('r2_OU', round(r2OU,3)), paste('del_AIC', round(delAIC,3))), bty='n')

hist(BM.fit$residuals, 20, main='BM residuals', las=1, border=0)
hist(OU.fit$residuals, 20, main='OU residuals', las=1, border=0)

# where BM is preferred:
set=2
d <- master.d[[set]]
module.of.interest <- 'B'
tmp <- d[d$module == module.of.interest ]

OU.fit <- nlreg(distance ~ (drift / (2*selection)) * (1-exp(-2*selection*divergence_time)), start = c(drift = 2, selection = 0.1), data=tmp)
coef(OU.fit)[1] -> dft
coef(OU.fit)[2] -> sel

BM.fit<- nlreg(distance ~ (slope*divergence_time) + intercept, start = c(slope = 0.006, intercept= 2.4), data=tmp)
coef(BM.fit)[1] -> slope
coef(BM.fit)[2] -> intercept

plot(distance ~ divergence_time, tmp, xlab = "divergence time (MY)", ylab = "metabolome distance", las=1, main=paste(setLabels[set], '\nModule', module.of.interest), col=4, pch=16)
curve(OU, 1, 50, las=1, add=T, lty=2)
curve(BM, 1, 50, las=1, add=T)

r2BM <- cor(tmp$distance, predict(BM.fit))^2 # r-squared 
r2OU <- cor(tmp$distance, predict(OU.fit))^2 # r-squared
delAIC <- AIC(BM.fit) - AIC(OU.fit) # delta AIC

legend('topleft', legend= c(paste('r2_BM', round(r2BM,3)), paste('r2_OU', round(r2OU,3)), paste('del_AIC', round(delAIC,3))), bty='n')

hist(BM.fit$residuals, 20, main='BM residuals', las=1, border=0)
hist(OU.fit$residuals, 20, main='OU residuals', las=1, border=0)


###############################################################################
# for figure:

par(mfrow=c(1, 2)) 

# plot example(s)
# where OU is preferred:
set=2
d <- master.d[[set]]
module.of.interest <- 'F'
tmp <- d[d$module == module.of.interest ]

OU.fit <- nlreg(distance ~ (drift / (2*selection)) * (1-exp(-2*selection*divergence_time)), start = c(drift = 2, selection = 0.1), data=tmp)
coef(OU.fit)[1] -> dft
coef(OU.fit)[2] -> sel

BM.fit<- nlreg(distance ~ (slope*divergence_time) + intercept, start = c(slope = 0.006, intercept= 2.4), data=tmp)
coef(BM.fit)[1] -> slope
coef(BM.fit)[2] -> intercept

plot(distance ~ divergence_time, tmp, xlab = "divergence time (MY)", ylab = "metabolome distance", las=1, main=paste(setLabels[set], '\nModule', module.of.interest), col=4, pch=16)
curve(OU, 1, 50, las=1, add=T, lty=2)
curve(BM, 1, 50, las=1, add=T)

# where BM is preferred:
set=2
d <- master.d[[set]]
module.of.interest <- 'B'
tmp <- d[d$module == module.of.interest ]

OU.fit <- nlreg(distance ~ (drift / (2*selection)) * (1-exp(-2*selection*divergence_time)), start = c(drift = 2, selection = 0.1), data=tmp)
coef(OU.fit)[1] -> dft
coef(OU.fit)[2] -> sel

BM.fit<- nlreg(distance ~ (slope*divergence_time) + intercept, start = c(slope = 0.006, intercept= 2.4), data=tmp)
coef(BM.fit)[1] -> slope
coef(BM.fit)[2] -> intercept

plot(distance ~ divergence_time, tmp, xlab = "divergence time (MY)", ylab = "metabolome distance", las=1, main=paste(setLabels[set], '\nModule', module.of.interest), col=4, pch=16)
curve(OU, 1, 50, las=1, add=T, lty=2)
curve(BM, 1, 50, las=1, add=T)


