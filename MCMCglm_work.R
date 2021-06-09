library(MCMCpack)
# BiocManager::install('MCMCglmm') [remove hash to install MCMCglmm]
library(MCMCglmm)
library(phytools)
library(tidyverse)
library(lme4)
library(gplots) 
library(plyr)
require(graphics) 
library(viridis)
library(limma)
library(lattice)
library(igraph)
library(FELLA)
library(org.Mm.eg.db)
library(KEGGREST)
library(magrittr)

# using MCMCglmm to estimate the effect (betas) of sex, age and their interaction, while correcting for phylogenetic covariance.

# this method will also estimate the varince due to species and to strain (variance in the case of random effects, and not betas for the fixed effects)

# [NOTE] p-values from MCMC (pMCMC) are: two times the smaller of the two quantities: MCMC estimates of i) the probability that a<0 or ii) the probability that a>0, where a is the parameter value. Its not a p-value as such, and better ways of obtaining Bayesian p-values exist.  Jarrod Hadfield https://stat.ethz.ch/pipermail/r-sig-mixed-models/2011q3/013925.html

setwd("/Volumes/GoogleDrive/My Drive/Documents/Hoffman_Metabolite_Phylogeny_Paper/Harrison_Hoffman_GitHub") 

flytree <- read.tree('fly.species.list.nwk')
plot(flytree)
axisPhylo(side = 1, root.time = NULL, backward = TRUE)

# a phylogeny must be invertable to use it as a covariance matrix in a glmm:
inverseA(flytree, nodes="TIPS", scale=TRUE) # phylogeny inversion (can't have polytomes)
flytree2 <-flytree
flytree2$edge.length <- flytree2$edge.length + 0.0001 # addinG a trivial amount of branch fixes polytomes, also 
inverseA(flytree2, nodes="TIPS", scale=TRUE)
flytree3 <- force.ultrametric(flytree2, method='extend') # the tree also needs to be untrametric the distance of all tips to the root should be the same
inv.phylo <- inverseA(flytree3, nodes="TIPS", scale=TRUE) 

plot(flytree3)
axisPhylo(side = 1, root.time = NULL, backward = TRUE)

# read metabolome data
read.table('comparative.mz.targeted.data.log_only', header=T, sep='\t') -> dat
head(dat) # a leading X is added to the metabolite names that begin with a number
as.matrix(dat[ ,18:ncol(dat)]) -> norm.dat

apply(norm.dat, 1, function(x) scale(x, scale=F)) -> norm.dat # center (but don't also scale) by sample
norm.dat[1:5,1:5]
t(norm.dat) -> norm.dat
norm.dat[1:5,1:4]
colnames(norm.dat) <- colnames(dat[18:ncol(dat)])

table(colSums(is.na(norm.dat))) # missingness 
norm.dat <- norm.dat[ ,colSums(is.na(norm.dat)) == 0]

dat[1:3,1:16]
dat[ ,c(2,2,5,8)] -> tmp 
head(tmp)
colnames(tmp)[2] <- 'strain' # phylogenetic methods book (Mod. Phylogenetic Comparative Methods.. authors: http://www.mpcm-evolution.com/authors) includes a tutorial for multiple measures (http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11/chapter-11-2-multiple-measurements-model-mcmcglmm. here i will replicate species names to represent the strain level. 

cbind(tmp, norm.dat) -> tmp
str(tmp)
tmp$sex <- as.factor(tmp$sex)
tmp$age.days <- as.factor(tmp$age.days)
# fix the metabolite names that lead with an X:
gsub("X", "\\1", colnames(tmp)) -> tmp.names
tmp.names
tmp.names [tmp.names=='anthosine'] <- "Xanthosine"
tmp.names [tmp.names=='anthurenate'] <- "Xanthurenate"
tmp.names [tmp.names=='anthine'] <- "Xanthine"
colnames(tmp) <- tmp.names
head(tmp)

# specify priors for each of species and strain:
species_strain_priors <- list(G=list(species=list(V=1,nu=0.02), strain=list(V=1,nu=0.02)), R=list(V=1,nu=0.02)) # We have told MCMCglmmto pay little heed to our prior expectaion (V) by specifying a small degree of belief parameter (nu) of 0.02.

colnames(tmp)[5:ncol(tmp)] -> mzs

convergence.result <- numeric()
pvals <- matrix(nrow=length(mzs), ncol=4)
post.mode.fixed <- matrix(nrow=length(mzs), ncol=4)
post.mode.rand <- matrix(nrow=length(mzs), ncol=3)
# to run models without strain (to test its effect by DIC):
species_priors <- list(G=list(species=list(V=1,nu=0.02)), R=list(V=1,nu=0.02)) 
deltaDIC_strain <- numeric()


#############################################
## WARNING: the next steps are computationally heavy:
#############################################

iterations <- 500000 # chose a number of MCMC iterations, you will probably want to start 'small' i.e. 100 iterations, and then evaluate the posterior distributions before chosing a large number of iterations.  


## run MCMC glmm:
for(i in 1:length(mzs)){
  tmp[ ,colnames(tmp) == mzs[i]] -> y
  o <- summary(m <- MCMCglmm(y ~ sex * age.days, random = ~ strain + species, family="gaussian", ginverse=list(species= inv.phylo$Ainv), prior=species_strain_priors, data= tmp, nitt=iterations, burnin=1000, thin=500))
  
  o$solutions[ ,5] -> pvals[i, ]
  posterior.mode(m$VCV) -> post.mode.rand[i, ]
  posterior.mode(m$Sol) -> post.mode.fixed[i, ]
  heidel.diag(m$Sol) -> conv.test
  sum(conv.test[ ,1], conv.test[ ,4]) -> convergence.result[i] # 1=pass, 0=fail for stationary test (first row) and for for halfwidth test (2nd row), all pass: sum=8, sum<8 means at least on of the estimates failed
  
  m_minus_strain <- MCMCglmm(y ~ sex * age.days, random = ~ species, family="gaussian", ginverse=list(species= inv.phylo$Ainv), prior=species_priors, data= tmp, nitt=iterations, burnin=1000, thin=500)
  m$DIC - m_minus_strain$DIC-> deltaDIC_strain[i]
}

paste0('P.val_', names(o$solutions[ ,5])) -> colnames(pvals)
paste0('post.mode_', names(o$solutions[ ,1])) -> colnames(post.mode.fixed)
paste0('post.mode_', names(posterior.mode(m$VCV))) -> colnames(post.mode.rand)

out <- as.data.frame(cbind(convergence.result, post.mode.fixed, pvals, post.mode.rand, deltaDIC_strain))
rownames(out) <- mzs
head(out)
str(out)

copy_of_out <- out
copy_of_out -> out
# replace mz codes in 'out':
head(out)
targeted_mz_info <- as.data.frame(read_csv("targeted.mz.info.csv"))
head(targeted_mz_info)

out[rownames(out) != 'internal_standard', ] -> out # remove internal standards (13C compounds that the Raftery lab added to diagnose LC-MS run)

out$metabolite <- rownames(out)

## clean up the metabolites names to match (this gets messy):
out$metabolite [!out$metabolite %in% mzs] # show all the mismatches
mzs [!mzs %in% out$metabolite]

out$metabolite -> fix
gsub('-', '.', fix) -> fix
gsub(' ', '.', fix) -> fix
out$metabolite <- fix
## 
head(out)

write.table(out, 'mcmc.glmm.output')

out <- read.table('mcmc.glmm.output', stringsAsFactors = F) 
head(out)

# pull out models that didn't converge and re-run:
table(out$convergence.result) # summary of which models converged
table(is.na(out$convergence.result))
out$metabolite[is.na(out$convergence.result)] -> stragglers
out[!is.na(out$convergence.result), ] -> out

out$metabolite[out$convergence.result<8] -> more.stragglers
c(stragglers, more.stragglers) -> stragglers

stragglers

# re-run the stragglers, 10x iterations:
convergence.result <- numeric()
pvals <- matrix(nrow=length(stragglers), ncol=4)
post.mode.fixed <- matrix(nrow=length(stragglers), ncol=4)
post.mode.rand <- matrix(nrow=length(stragglers), ncol=3)
# to run models without strain (to test its effect by DIC):
species_priors <- list(G=list(species=list(V=1,nu=0.02)), R=list(V=1,nu=0.02)) 
deltaDIC_strain <- numeric()

stragglers[!stragglers %in% colnames(tmp)]
table(stragglers %in% colnames(tmp)) # all 63 stragglers are in tmp:

for(i in 1:length(stragglers)){
  tmp[ ,colnames(tmp) == stragglers[i]] -> y
  o <- summary(m <- MCMCglmm(y ~ sex * age.days, random = ~ strain + species, family="gaussian", ginverse=list(species= inv.phylo$Ainv), prior=species_strain_priors, data= tmp, nitt=5000000,burnin=10000,thin=10000))
  
  o$solutions[ ,5] -> pvals[i, ]
  posterior.mode(m$VCV) -> post.mode.rand[i, ]
  posterior.mode(m$Sol) -> post.mode.fixed[i, ]
  heidel.diag(m$Sol) -> conv.test
  sum(conv.test[ ,1], conv.test[ ,4]) -> convergence.result[i] # 1=pass, 0=fail for stationary test (first row) and for for halfwidth test (2nd row), all pass: sum=8, sum<8 means something failed
  
  m_minus_strain <- MCMCglmm(y ~ sex * age.days, random = ~ species, family="gaussian", ginverse=list(species= inv.phylo$Ainv), prior=species_priors, data= tmp, nitt=5000000,burnin=10000,thin=10000)
  m$DIC - m_minus_strain$DIC-> deltaDIC_strain[i]
}

paste0('P.val_', names(o$solutions[ ,5])) -> colnames(pvals)
paste0('post.mode_', names(o$solutions[ ,1])) -> colnames(post.mode.fixed)
paste0('post.mode_', names(posterior.mode(m$VCV))) -> colnames(post.mode.rand)

rerun <- as.data.frame(cbind(convergence.result, post.mode.fixed, pvals, post.mode.rand, deltaDIC_strain))
rownames(rerun) <- stragglers
rerun$metabolite <- rownames(rerun)
head(rerun)
str(rerun)

write.table(rerun, 'mcmc.glmm.on.the.stragglers') # save the results!
read.table('mcmc.glmm.on.the.stragglers', header=T) -> rerun

table(is.na(rerun$post.mode_sexM))
table(rerun$convergence.result)
table(is.na(rerun$convergence.result))

#### CAREFUL HERE< the straggler df may need the mz codes and the metabolite column added::::::::
# add the new re-run values to the output
head(out)
rownames(out) <- out$metabolite
out[!is.na(out$convergence.result), ] -> out
out[out$convergence.result==8, ] -> willing
names(willing)-> names(rerun)
rbind(willing, rerun) -> out
table(is.na(out$convergence.result))
table(out$convergence.result) # models coverged for at least 7 of 8 parameters for 105 of 106 metabolites.

write.table(out, 'mcmc.glmm.output_added_stragglers') #save the output

#########################################################################################################
out <- read.table('mcmc.glmm.output_added_stragglers') # read the output
head(out)
 
table(p.adjust(out$P.val_sexM, method='BH') < 0.05)
table(p.adjust(out$P.val_age.days31, method='BH') <0.05)
table(p.adjust(out$P.val_sexM.age.days31, method='BH') <0.05)

head(out)
# relationship between variance attributed to species, strain and residuals ('units').
pairs(out[ ,c(10:13)], log='xy')

# distribution of deltaDIC (the 'information about mz added' to the model by strain)I believe a negative value here means that adding strain makes a 'better' model, given the 'cost' of an additional parameter.  Warning: DIC is only valid when the posterior distribution is approximately multivariate normal. Some posterior distributions were not normal, when the parameter was near zero
par(mfrow=c(1,1))
h <- hist(out$deltaDIC_strain, 20)
cuts <- cut(h$breaks, c(-Inf, -0.1, +Inf)) # trick to color part of the histogram
plot(h, col=c('grey', 0)[cuts], las=1)

head(out)
p.adjust(out$P.val_sexM, method='BH') -> out$fdr.sex
p.adjust(out$P.val_age.days31, method='BH') -> out$fdr.age
p.adjust(out$P.val_sexM.age.days31, method='BH') -> out$fdr.sexbyage

out[out$fdr.sex <= 0.05 | out$fdr.age <= 0.05 | out$fdr.sexbyage <= 0.05, ] -> sigs

head(sigs)

# venn diagram
as.numeric(out$fdr.sex <= 0.05) -> sex
as.numeric(out$fdr.age <= 0.05) -> age
as.numeric(out$fdr.sexbyage <= 0.05) -> sex_by_age
cbind(sex, age, sex_by_age) -> fdr0.05
par(mfrow=c(1,1))
vennDiagram(fdr0.05, circle.col=viridis(3), names=c('sex', 'age', 'sex x age'))

head(sigs)
as.matrix(abs(sigs[ ,3:5])) -> fixeffmat
rownames(fixeffmat) <- sigs$metabolite
colnames(fixeffmat) <- c('sex', 'age', 'sex x age')

my_palette <- colorRampPalette(c("lightgrey", "blue"))(n = 10)

heatmap.2(fixeffmat, Colv=F, dendrogram='row', scale='none', trace='none',cexCol=1, cexRow = 0.5, margins=c(4,30), density.info = 'none', key.xlab=expression(beta), srtCol=45, col=my_palette)
#########################################################################################################



####################################################################################
### for the untargeted data:
####################################################################################
# read metabolome data
glob.dat <- read.table('combined_modes_untargeted_inputed_590.features', header=T, stringsAsFactors = T)
glob.dat[1:4,1:10]
str(glob.dat)

colnames(glob.dat)[4:ncol(glob.dat)] -> features


convergence.result <- numeric()
pvals <- numeric()
post.mode.fixed <- numeric()
post.mode.rand <- matrix(nrow=length(features), ncol=2)

for(i in 1:length(features)){
  glob.dat[ ,colnames(glob.dat) == features[i]] -> y
  o <- summary(m <- MCMCglmm(y ~ sex, random = ~ strain + species, family="gaussian", ginverse=list(species= inv.phylo$Ainv), prior=species_strain_priors, data= glob.dat, nitt=500000, burnin=1000, thin=500))
  o$solutions[2,5] -> pvals[i]
  posterior.mode(m$Sol)[2] -> post.mode.fixed[i]
  posterior.mode(m$VCV)[1:2] -> post.mode.rand[i, ]
  
  heidel.diag(m$Sol) -> conv.test
  sum(conv.test[ ,1], conv.test[ ,4]) -> convergence.result[i] # 1=pass, 0=fail for stationary test (first row) and for for halfwidth test (2nd row), all pass: sum=4, sum<4 means at least on of the estimates failed
}

globout <- as.data.frame(cbind(post.mode.rand, post.mode.fixed, pvals, convergence.result))
rownames(globout) <- features
colnames(globout) <- c('strain_effect', 'species_effect', 'sex_effect', 'sex_P', 'MCMC_convergence')
head(globout)


write.table(globout, 'imputed.untargeted.MCMCglmm_sex')

globout <- read.table('imputed.untargeted.MCMCglmm_sex') 
head(globout)

table(globout$MCMC_convergence)

globout$fdr.sex <- p.adjust(globout$sex_P, method='BH')
table(globout$fdr.sex < 0.1)

convergence.result2 <- numeric()
pvals2 <- numeric()
post.mode.fixed2 <- numeric()
post.mode.rand2 <- matrix(nrow=length(features), ncol=2)

# limit the 5x10^6 iterations below to those features with and initial FDR of <0.1:
iterations <- 5000000
head(globout)
second_round_features <- rownames(globout[globout$fdr.sex < 0.1, ])

for(i in 1:length(second_round_features)){
  glob.dat[ ,colnames(glob.dat) == second_round_features[i]] -> y
  o <- summary(m <- MCMCglmm(y ~ sex, random = ~ strain + species, family="gaussian", ginverse=list(species= inv.phylo$Ainv), prior=species_strain_priors, data= glob.dat, nitt=iterations, burnin=10000, thin=500))
  o$solutions[2,5] -> pvals2[i]
  posterior.mode(m$Sol)[2] -> post.mode.fixed2[i]
  posterior.mode(m$VCV)[1:2] -> post.mode.rand2[i, ]
  
  heidel.diag(m$Sol) -> conv.test
  sum(conv.test[ ,1], conv.test[ ,4]) -> convergence.result2[i] # 1=pass, 0=fail for stationary test (first row) and for for halfwidth test (2nd row), all pass: sum=4, sum<4 means at least one of the estimates failed
}

second_round <- cbind(post.mode.fixed2, pvals2, convergence.result2)
rownames(second_round) <- second_round_features

globout2 <- as.data.frame(second_round)
head(globout2)
colnames(globout2) <- c('sex_effect', 'sex_P', 'MCMC_convergence')

table(globout2$MCMC_convergence) # all converged

head(globout)
head(globout2)

sexP <- globout$sex_P
names(sexP) <- rownames(globout)
sexP2 <- globout2$sex_P
names(sexP2) <- rownames(globout2)

plot(-log(sexP[names(sexP2)], 10), -log(sexP2, 10)) #  
abline(0,1)
min(sexP)
min(sexP2) # 10x iterations in the second round allowed Pvalue that are 10x smaller, and presumably more precise

sexP[names(sexP2)] <- sexP2 # replace 1st round Ps with 2nd round Ps
plot(-log(globout$sex_P, 10), -log(sexP, 10), pch=16)
qqman::qq(sexP) # still 'maxed' out

globout$sex_P <- sexP
globout$fdr.sex <- p.adjust(sexP, method='BH')
table(globout$fdr.sex < 0.05) # this second round seems to have better resolved (reduced) the features at metabolome-wide FDR<0.1

colSums(is.na(globout))

# remember to save the output!!
# write.table(globout, 'imputed.untargeted.MCMCglmm_sex_up.to.5.million.iterations')


### read the data and continue
globout <- read.table('imputed.untargeted.MCMCglmm_sex_up.to.5.million.iterations') 
head(globout)

# write an input file for mummichog:
globout[grepl('pos', rownames(globout)), c(1:4,6)] -> positive.features
globout[grepl('neg', rownames(globout)), c(1:4,6)] -> negative.features


# merge retention times:
read.table('Negative ion mode metabolite info.csv', sep=',' , stringsAsFactors = F, header=T) -> neg.retention.times
head(neg.retention.times)
paste0('neg_', neg.retention.times$Mass) -> neg.retention.times$feature

rownames(negative.features) -> negative.features$feature
head(negative.features)
head(neg.retention.times)
table(negative.features$feature %in% neg.retention.times$feature)
merge(negative.features, neg.retention.times, all.x=T, by='feature') -> negative.features
head(negative.features)

# columns required for mummichog: c('mz', 'rtime', 'p-value', 't-score')
negative.features[ ,c('Mass', 'Retention.Time', 'sex_P', 'sex_effect')] -> neg.sex.input 
head(neg.sex.input)
colnames(neg.sex.input) <- c('mz', 'rtime', 'p-value', 't-score')

colSums(is.na(neg.sex.input))
rowSums(is.na(neg.sex.input))
neg.sex.input <- neg.sex.input[rowSums(is.na(neg.sex.input)) == 0, ]

write.table(neg.sex.input, file='mummichog-1.0.9/phylo_mz/neg.sex.input.txt', sep='\t', row.names = F)

# merge retention times:
read.table('Positive ion mode metabolite info.csv', sep=',' , stringsAsFactors = F, header=T) -> pos.retention.times
head(pos.retention.times) 
paste0('pos_', pos.retention.times$Mass) -> pos.retention.times$feature

rownames(positive.features) -> positive.features$feature
head(positive.features)
head(pos.retention.times)
table(positive.features$feature %in% pos.retention.times$feature)
merge(positive.features, pos.retention.times, all.x=T, by='feature') -> positive.features
head(positive.features)

# columns required: c('mz', 'rtime', 'p-value', 't-score')
positive.features[ ,c('Mass', 'Retention.Time', 'sex_P', 'sex_effect')] -> pos.sex.input 
head(pos.sex.input)
colnames(pos.sex.input) <- c('mz', 'rtime', 'p-value', 't-score')
head(pos.sex.input)

write.table(pos.sex.input, file='mummichog-1.0.9/phylo_mz/pos.sex.input.txt', sep='\t', row.names = F)

# find a Pvalue threshold to pass to mummichog that reflects the metabolome-wide signifcance level:
table(globout$fdr.sex <= 0.05) # FDR <= 5%
max(globout$sex_P[globout$fdr.sex <= 0.05]) # maximum Pvalue among the metabolites with a sex effect at FDR <= 5%

