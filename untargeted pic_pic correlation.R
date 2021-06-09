# Examine covariance among the metabolite PICs; i.e. what is the covariance of metabolites across the 50 million year old phylogeny?

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.11")
BiocManager::install("limma")
BiocManager::install('WGCNA')
BiocManager::install('org.Mm.eg.db')
BiocManager::install('HEMDAG')

library(plyr) # data wrangling
library(splitstackshape)
library(dplyr)
library(tidyr)
library(reshape2)
library(readr)
library(RVenn)
library(purrr)
library(DMwR)
library(impute)
library(AssocTests)

library(sjstats) # modeling tools
library(lme4)
library(lmerTest)
library(nlme) 
library(MuMIn) 
library(MASS)
library(smatr)
library(lmodel2)

library(WGCNA) # cluster ID and pathway tools
library(org.Mm.eg.db)
library(FELLA)
library(HEMDAG)

library(graphics) # plotting
library(ggplot2) 
library(gplots) 
library(corrplot)
library(RColorBrewer)
library(limma)
library(igraph)
library(pheatmap)

library(ape) # phylogenetic analysis
library(phytools) 
library(geiger) 
library(caper) 


mycol7 <- brewer.pal(7, 'RdYlBu')
mycol4 = c("#EEC900", "#CD6600", "#87CEFA", "#104E8B")

setwd("~/Google Drive/My Drive/Documents/Hoffman_Metabolite_Phylogeny_Paper") 

flytree <- read.tree('fly.species.list.nwk')
plot(flytree)
axisPhylo(side = 1, root.time = NULL, backward = TRUE)


dat <- read.table('combined_modes_untargeted_inputed_590.features', header=T, sep=' ')
dat[1:6,1:10] 

features <- colnames(dat)[4:ncol(dat)]
feats <- as.matrix(dat[ ,4:ncol(dat)])

# a look at the data
par(mfrow=c(2,2))
boxplot(feats[order(apply(feats, 1, median)), ], use.cols = F, main='log, mean-centered', xlab='sample', las=1)
boxplot(feats[ ,order(colMeans(feats))], use.cols = T, main='log, mean-centered', xlab='metabolite', las=1)

######################################################################################################
# assess the covariance of metabolites, after correcting for phylogenetic correlation:
## PIC including intraspecies covariance

sex <- dat$sex

for(k in unique(sex)) {
  tmp.dat <- dat[sex == k, ]
  no.species <- length(unique(tmp.dat$species))
  
  tmpics <- matrix(ncol=no.species-1, nrow=length(features))
  
  for(i in 1:length(features)) {
    feature <- features[i]
    tmp <- tmp.dat[ ,feature]
    names(tmp) <- tmp.dat$species
    tmp <- tmp[!is.na(tmp)]
    tmp <- aggregate(tmp, by=list(names(tmp)), list)
    names(tmp$x) <- tmp$Group.1
    tmp  <- tmp$x
    tmp.tree <- keep.tip(flytree, names(tmp))
    tmp  <- tmp[tmp.tree$tip.label]
    tmpics[i, ] <- pic.ortho(tmp, tmp.tree, intra = T, var.contrasts = F) }
  # method from Felsenstein (2008) Am. Nat.
  rownames(tmpics) <- features
  write.table(tmpics, file=paste('~/Google Drive/My Drive/Documents/Hoffman_Metabolite_Phylogeny_Paper/correlation among the metabolome/pics_untargeted', k))
}

picM <- read.table('correlation among the metabolome/pics_untargeted males')
picF <- read.table('correlation among the metabolome/pics_untargeted females')

features <- rownames(picM)

# add lifespan data:
picM5 <- read.table('correlation among the metabolome/pics M 5')
picF5 <- read.table('correlation among the metabolome/pics F 5')

picM <- rbind(picM5[2, ], picM)
picF <- rbind(picF5[2, ], picF)

##############################################################################

nSets = 2
PICs <- list(picF, picM)
setLabels = c("female", "male")
shortLabels = c('F', 'M')

# plot correlation matrices of metabolites for each sex x age:
F_pal <- colorRampPalette(c("grey", "white", "blue"))(n = 50)
M_pal <- colorRampPalette(c("grey", "white", "brown"))(n = 50)

setPals <- list(F_pal, M_pal)

for(set in 1:nSets) {
  tmpic <- apply(PICs[[set]], 2, unlist)
  pic.cor <- cor(t(tmpic[ ,5:ncol(tmpic)]), method='pearson')
  heatmap.2(pic.cor, trace='n', col=unlist(setPals[set]), key.title = 'cor PIC-PIC', margins = c(10,12), key.xlab = 'pearson (r)', cexRow=0.5, cexCol=0.5, keysize=1, dendrogram='none', main=paste('PIC correlation', setLabels[set])) }
##############################################################################

# PERMUTATION test: compare the distribution of correlation among feature PICs: test distribution of correlations between PICs by permuting.  In this case, try conservative permutation in which only the species labels are randomizes, after species means are calculated and, most conservatively, while conserving the metabolite-metabolite relationships within each species

traits <- names(dat)[4:ncol(dat)] # traits = LCMS features 

tmp.dat <- dat[dat$sex == 'females', ]
no.species <- length(table(tmp.dat$species))

tmpics <- matrix(ncol=no.species-1, nrow=length(traits))

par(mfrow=c(1,1))
set <- c(1:2)[setLabels == 'female']
tmpic <- apply(PICs[[set]], 2, unlist)
pic.cor <- cor(t(tmpic[ ,5:ncol(tmpic)]), method='pearson')
d <- density(pic.cor[lower.tri(pic.cor, diag=F)])
plot(d, ylim=c(0, 1), las=1, xlab='metabolite PIC, pariwise correlation', col=mycol4[2], lwd=2, ylab='kernel density', main='')

for(i in 1:20) {
  randomised.species <- sample(unique(tmp.dat$species)) # conservative permutation, where species names are shuffled prior to PIC calculation, but after species means are estimated.  THis is also keeping all metabolites within a species together [ this does not disrupt the mz-mz relationship, just the effect of phylogeny]
  
  for(i in 1:length(traits)) {
    trait <- traits[i]
    tmp <- tmp.dat[ ,trait]
    names(tmp) <- tmp.dat$species
    tmp <- tmp[!is.na(tmp)]
    tmp <- aggregate(tmp, by=list(names(tmp)), list)
    names(tmp$x) <- randomised.species
    tmp  <- tmp$x
    tmp.tree <- keep.tip(flytree, names(tmp))
    tmp  <- tmp[tmp.tree$tip.label]
    tmpics[i, ] <- pic.ortho(tmp, tmp.tree, intra = T, var.contrasts = F) }
  # method from Felsenstein (2008) Am. Nat.
  rownames(tmpics) <- traits

  corperm <- cor(t(tmpics), method='pearson')
  lines(density(corperm[lower.tri(corperm, diag=F)]), col=rgb(0,0,0,0.5))  }


# and the males:
tmp.dat <- dat[dat$sex == 'males', ]
no.species <- length(table(tmp.dat$species))

tmpics <- matrix(ncol=no.species-1, nrow=length(traits))

par(mfrow=c(1,1))
set <- c(1:2)[setLabels == 'male']
tmpic <- apply(PICs[[set]], 2, unlist)
pic.cor <- cor(t(tmpic[ ,5:ncol(tmpic)]), method='pearson')
d <- density(pic.cor[lower.tri(pic.cor, diag=F)])
plot(d, ylim=c(0, 1), las=1, xlab='metabolite PIC, pariwise correlation', col=mycol4[4], lwd=2, ylab='density', main='')

for(i in 1:20) {
  randomised.species <- sample(unique(tmp.dat$species)) # conservative permutation, where species names are shuffled prior to PIC calculation, but after species means are estimated.  THis is also keeping all metabolites within a species together [ this does not disrupt the mz-mz relationship, just the effect of phylogeny]
  
  for(i in 1:length(traits)) {
    trait <- traits[i]
    tmp <- tmp.dat[ ,trait]
    names(tmp) <- tmp.dat$species
    tmp <- tmp[!is.na(tmp)]
    tmp <- aggregate(tmp, by=list(names(tmp)), list)
    names(tmp$x) <- randomised.species
    tmp  <- tmp$x
    tmp.tree <- keep.tip(flytree, names(tmp))
    tmp  <- tmp[tmp.tree$tip.label]
    tmpics[i, ] <- pic.ortho(tmp, tmp.tree, intra = T, var.contrasts = F) }
  # method from Felsenstein (2008) Am. Nat.
  rownames(tmpics) <- traits
  
  tmpics[1:4,1:4]
  
  corperm <- cor(t(tmpics), method='pearson')
  lines(density(corperm[lower.tri(corperm, diag=F)]), col=rgb(0,0,0,0.5))  }

rm(corperm)
rm(pic.cor)
rm(tmp.pic)
rm(tmp.dat)
rm(randomised.species)
##############################################################################


##############################################################################
# plot a few of the top ls ~ mz PIC correlations
par(mar=c(5,4,4,2))
par(mfrow=c(4,3))

for(set in 1:nSets){

  tmpic <- apply(PICs[[set]], 2, unlist)
  pic.cor <- cor(t(tmpic), method='pearson')
  
  trait.of.interest <- 'mean.ls'
  
  tops <- names(head(pic.cor[ ,trait.of.interest][rev(order(pic.cor[ ,trait.of.interest]))]))
  ifelse(length(-grep('weight', tops)) > 0, tops <- tops[-grep('weight', tops)], tops <- tops)
  ifelse(length(-grep('ls', tops)) > 0, tops <- tops[-grep('ls', tops)], tops <- tops)
  bots <- names(tail(pic.cor[ ,trait.of.interest][rev(order(pic.cor[ ,trait.of.interest]))]))
  hits <- c(tops, bots)
  hits
  
  # plot the 'top' 3 'hits'; the features with the largest (or smallest) correlation coefficients [not the same as the lowest P]  for each set:
  for(i in 1:3){
    mz.of.interest <- hits[i]
    
    plot(tmpic[trait.of.interest, ] ~  tmpic[mz.of.interest, ], pch=19, ylab=paste('PIC', trait.of.interest), xlab=paste('PIC', mz.of.interest), las=1, main=setLabels[set], col=mycol4[set])
    major.axis.regression <- sma(tmpic[trait.of.interest, ] ~  tmpic[mz.of.interest, ] -1)
    pvalue <- paste0('P=', round(unlist(major.axis.regression$pval), 3))
    fdr <- paste0('FDR=', round(p.adjust(c(major.axis.regression$pval, runif(length(features)-1, min=0.1, max=1)), method='BH')[1], 3))
    legend(ifelse(unlist(major.axis.regression$coef)[2] >0, 'bottomright', 'topright'), legend=c(pvalue, fdr), bty='n')
    abline(sma(tmpic[trait.of.interest, ] ~  tmpic[mz.of.interest, ]-1))}
}
# seeing a few metabolites whose PICs correlate with lifespan trait PICs




#############################################################
# try ls ~ mz univariate major axis regression models 
#############################################################

ls.by.mz <- function(x) {
  m <- sma(ls ~ x -1)
  return(c(m$r2, m$pval)) }

ls_mz <- list()

for(set in 1:nSets){
  tmpic <- apply(PICs[[set]], 2, unlist)  
  ls <- t(tmpic)[ ,'mean.ls']
  mz.pics <- t(tmpic)[ ,2:nrow(tmpic)]
  
  out <- matrix(nrow=nrow(PICs[[set]])-1, ncol=2)
  for(i in 1:length(features)) {
    mz <- mz.pics[ ,i]
    out[i, ] <- unlist(ls.by.mz(mz)) }
  
  out <- as.data.frame(out)
  rownames(out) <- features
  colnames(out) <- c('r.sq','pval')
  out$FDR <- p.adjust(out$pval, method='BH') 
  ls_mz[[set]] <- out
  
  write.table(out, file=paste0('~/Google Drive/My Drive/Documents/Hoffman_Metabolite_Phylogeny_Paper/correlation among the metabolome/lifespanPIC_by_mzPIC_major.axis.regression_untargeted_', setLabels[set], '.csv'), sep=',') }


ls_mz[[1]]

ls_by_mz.table <- cbind(ls_mz[[1]], ls_mz[[2]])
colnames(ls_by_mz.table) <- paste(colnames(ls_by_mz.table), setLabels[c(1,1,1,2,2,2)])

ls_by_mz.table$mz <- rownames(ls_by_mz.table)
ls_by_mz.table[1:4, ]

ls_by_mz.table <- ls_by_mz.table[ ,c(7, 1:6)]

table(ls_by_mz.table$`FDR female` <0.5 | ls_by_mz.table$`FDR male` < 0.5)
table(ls_by_mz.table$`FDR female` <0.25 | ls_by_mz.table$`FDR male` < 0.25)
table(ls_by_mz.table$`FDR female` <0.1 | ls_by_mz.table$`FDR male` < 0.1)
table(ls_by_mz.table$`FDR female` <0.05 | ls_by_mz.table$`FDR male` < 0.05)

write.table(ls_by_mz.table, file='~/Google Drive/My Drive/Documents/Hoffman_Metabolite_Phylogeny_Paper/correlation among the metabolome/lifespanPIC_by_untargeted_mzPIC_major.axis.regression_Table.csv', sep=',', row.names=F) 


rsq.mat <- as.matrix(cbind(ls_mz[[1]][1], ls_mz[[2]][1]))
fdr.mat <- as.matrix(cbind(ls_mz[[1]][3], ls_mz[[2]][3]))
colnames(rsq.mat) <- setLabels
colnames(fdr.mat) <- setLabels

colfunc <- colorRampPalette(c('white', 'blue'))

heatmap.2(rsq.mat[rowSums(fdr.mat < 0.1) > 0, ], margin=c(8,8), trace='n', col=colfunc(20), Colv=F, key.xlab = 'r-squared', density.info = 'n', keysize = 0.8, key.title = '', cexRow=1, cexCol=1)

rsq.mat[rowSums(fdr.mat < 0.1) > 0, ]

### plot the 'fdr curve':  
fdrs <- unlist(c(ls_mz[[1]][3], ls_mz[[2]][3]))
table(fdrs <0.5)
par(mar=c(5,4,4,2))
par(mfrow=c(1,1))

plot(unlist(ls_mz[[1]][3])[order(unlist(ls_mz[[1]][3]))], pch=16, las=1, main='ls ~ LC-MS feature FDR (n=598 comparisons)', xlab='LC/MS feature rank', ylab='FDR', col=2)
points(unlist(ls_mz[[2]][3])[order(unlist(ls_mz[[2]][3]))], col=4)
legend('topleft', pch=16, col=c(2,4), legend=c('female', 'male'), bty='n')


### plot the 'top hit' from each set:
par(mfrow=c(1,2))
par(xpd=FALSE)

for(set in 1:nSets){
  minFDR <- min(ls_mz[[set]][3])
  hit <- rownames(ls_mz[[set]][ls_mz[[set]][3] == minFDR, ]) # find rowname (lc-ms feature) with the minimum FDR value
  plot(unlist(PICs[[set]]['mean.ls', ]) ~  unlist(PICs[[set]][hit, ]), pch=19, ylab=paste('PIC of mean lifespan'), xlab=paste('PIC of', hit), las=1, main=paste("top 'hit' from", setLabels[set], 'untargeted metabolome')) 
  legend('top', legend=paste('FDR=', round(minFDR, 4)), bty='n')
  abline(sma(unlist(PICs[[set]]['mean.ls', ]) ~  unlist(PICs[[set]][hit, ]) -1) )}


#############################################################

########################################################################################################################
# send the the associated features to die inside the mummichog!!!
########################################################################################################################

# to get LC-MS feature information: retention time and mass/charge ratio:
neg.retention.times <- read.table('~/Google Drive/My Drive/Documents/Hoffman_Metabolite_Phylogeny_Paper/Negative ion mode metabolite info.csv', sep=',' , stringsAsFactors = F, header=T) 
head(neg.retention.times)
neg.retention.times$feature <- paste0('neg_', neg.retention.times$Mass)
head(neg.retention.times)

pos.retention.times <- read.table('~/Google Drive/My Drive/Documents/Hoffman_Metabolite_Phylogeny_Paper/Positive ion mode metabolite info.csv', sep=',' , stringsAsFactors = F, header=T) 
head(pos.retention.times) 
pos.retention.times$feature <- paste0('pos_', pos.retention.times$Mass)
head(pos.retention.times) 

## look to see how many of the LC-MS features that have been analysed above are included in the feature meta data:
table(features %in% neg.retention.times$feature) 
table(features %in% pos.retention.times$feature)
features[!features %in% pos.retention.times$feature & !features %in% neg.retention.times$feature] # one of the 589 LC-MS features is not present in the feature meta data 

pos.retention.times <- pos.retention.times[pos.retention.times$feature %in% features, ] # cut down meta data to only those analyzed here (the ones measured in all strains)
neg.retention.times <- neg.retention.times[neg.retention.times$feature %in% features, ] # cut down meta data to only those analyzed here (the ones measured in all strains)

meta.data <- rbind(pos.retention.times, neg.retention.times)

head(meta.data)
head(ls_by_mz.table)

lstab <- merge(meta.data, ls_by_mz.table, by.x='feature', by.y='mz')
head(lstab)

lstab$mode <- as.factor(ifelse(1:length(lstab$feature) %in% grep('pos', lstab$feature), 'positive', 'negative'))
head(lstab)

# columns required for mummichog input: c('mz', 'rtime', 'p-value', 't-score')

# make input file for the female data [NOTE there were few associations in the male data, so just do the femme data]:
female_neg_mummichog_input <- lstab[lstab$mode == 'negative' ,c('Mass', 'Retention.Time', 'pval female', 'r.sq female')]
head(female_neg_mummichog_input)
colnames(female_neg_mummichog_input) <- c('mz', 'rtime', 'p-value', 't-score')

write.table(female_neg_mummichog_input , file='~/Google Drive/My Drive/Applications/mummichog-1.0.9/phylo_mz/fem_neg_lifespan_mummichog_input.txt', sep='\t', row.names = F)

female_pos_mummichog_input <- lstab[lstab$mode == 'positive' ,c('Mass', 'Retention.Time', 'pval female', 'r.sq female')]
colnames(female_pos_mummichog_input) <- c('mz', 'rtime', 'p-value', 't-score')
head(female_pos_mummichog_input)

write.table(female_pos_mummichog_input, file='~/Google Drive/My Drive/Applications/mummichog-1.0.9/phylo_mz/fem_pos_lifespan_mummichog_input.txt', sep='\t', row.names = F)




########################################################################################################################
########################################################################################################################


#############################################################
# multivariate association.  What about ls ~ metabolite modules?
#############################################################

## use WGCNA to pick modules
setwd("~/Google Drive/My Drive/Documents/Hoffman_Metabolite_Phylogeny_Paper")
options(stringsAsFactors = FALSE)

# Form a multi-set structure that will hold the metabolite data and the organismal phenotypes (i.e. lifespan and weight).
multiPIC = vector(mode = "list", length = nSets) # empty entity

multiPIC[[1]] = list(data = as.data.frame(t(picF[-c(1), ]))) # a multi-set list holding the lc-ms feature PICs for each group (sex in this case)
names(multiPIC[[1]]$data) = rownames(picF)[2:nrow(picF)]

multiPIC[[2]] = list(data = as.data.frame(t(picM[-c(1), ]))) # a multi-set list holding the lc-ms feature PICs for each group (sex in this case)
names(multiPIC[[2]]$data) = rownames(picM)[2:nrow(picM)]

Traits = vector(mode="list", length = nSets)
Traits[[1]] <- t(picF[1, ])
Traits[[2]] <- t(picM[1, ])

collectGarbage()

exprSize = checkSets(multiPIC) # Check that the data has the correct format for many functions operating on multiple sets (I think this will give an error if there is a problem)

sampleTrees = list()
for (set in 1:nSets)
{sampleTrees[[set]] = hclust(dist(multiPIC[[set]]$data), method = "average")}

par(mfrow=c(2,2))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets)
  plot(sampleTrees[[set]], main = paste("Sample clustering in", setLabels[set]), xlab="", sub="", cex = 0.7, las=1)

collectGarbage()

# Check the size of the leftover data
exprSize = checkSets(multiPIC)

# Define data set dimensions
nGenes = exprSize$nGenes
nSamples = exprSize$nSamples

save(multiPIC, Traits, nGenes, nSamples, setLabels, shortLabels, exprSize, file = "correlation among the metabolome/untargeted-metabolome-WGCNA-dataInput.RData") # saves the data build so far

load("correlation among the metabolome/untargeted-metabolome-WGCNA-dataInput.RData") # load the data

nSets = checkSets(multiPIC)$nSets

# Evaluate a set of soft-thresholding powers
#################################################################
# Constructing a weighted gene network entails the choice of the soft thresholding power β to which co-expression similarity is raised to calculate adjacency [3]. The authors of [3] have proposed to choose the soft thresholding power based on the criterion of approximate scale-free topology. We refer the reader to that work for more details; here we illustrate the use of the function pickSoftThreshold that performs the analysis of network topology and aids the user in choosing a proper soft-thresholding power. The user chooses a set of candidate powers (the function provides suitable default values), and the function returns a set of network indices that should be inspected.

# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets)

powers = c(seq(4,10,by=1), seq(12,20, by=2)) # make a set of power levels to calculate over

# Call the network topology analysis function for each set in turn
for (set in 1:nSets)
  powerTables[[set]] = list(data = pickSoftThreshold(multiPIC[[set]]$data, powerVector=powers, verbose = 2)[[2]])

collectGarbage()

# Plot the results:
colors = c(1:nSets)
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity")
# Get the minima and maxima of the plotted points

ylim = matrix(NA, nrow = 2, ncol = 2)

for (set in 1:nSets){
  for (col in 1:length(plotCols)) {
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE)
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE) }}

# Plot model fit and mean connectivity over soft power: 
par(mfrow = c(1, 2))
par(mar = c(4.2, 4.2 , 2.2, 0.5))

for (col in 1:length(plotCols)) for (set in 1:nSets)
{
  if (set==1)
  { plot(powerTables[[set]]$data[ ,1], -sign(powerTables[[set]]$data[ ,3]) * powerTables[[set]]$data[ ,2],
         xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
         main = colNames[col]) }
  if (col==1)
  { text(powerTables[[set]]$data[ ,1], -sign(powerTables[[set]]$data[ ,3])*powerTables[[set]]$data[ ,2], labels=powers, col=colors[set])
  } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[ ,plotCols[col]], labels=powers, col=colors[set])
  if (col==1)
  { legend("bottomright", legend = setLabels, col = colors, pch = 20) 
  } else
    legend("topright", legend = setLabels, col = colors, pch = 20) 
}


par(mfrow = c(1,1))

for (set in 1:nSets) {
  if (set==1) {
    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],  xlab="Soft Threshold (cor^x)", type="b", ylim=c(-0.3,1), col=set, pch=16, las=1, xaxp  = c(0, 20, 20))
  } else 
    lines(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],type="b", col=set, pch=16) 
}
legend('topleft', legend=setLabels, col=1:nSets, pch=16)


#####################################
softPower = 6 # chosen power 
#####################################

# how to check the scale-freeness of the network?  Scale-free = linear relationship between: frequency of log(conectivity) ~ log(conectivity), a.k.a a 'power law' relationship.  Do this check once the network is made:

# Network construction starts by calculating the adjacencies in the individual sets, using the soft thresholding power
# build an empty array to hold the adjacencies
adjacencies = array(0, dim = c(nSets, nGenes, nGenes));

# Calculate adjacencies in each individual data set
for (set in 1:nSets)
  adjacencies[set, , ] = abs(cor(multiPIC[[set]]$data, use = "p"))^softPower # calculates the adjacency matrix, in this case the pearson correlation raised to the 'softpower'.  I don't know the purpose of the absolute value [abs()] part of this line, since x ^ power is positive.

# [optional]: analyze network properties, e.g. degree distribution:
# Initialize an appropriate array to hold the topological overlap measures (TOM)s
TOM = array(0, dim = c(nSets, nGenes, nGenes))

# Calculate TOMs in each individual data set
for (set in 1:nSets)
  TOM[set, , ] = TOMsimilarity(adjacencies[set, , ]) # Calculation of the topological overlap matrix, and the corresponding dissimilarity, from a given adjacency matrix.

###############################################################################################################################################################################
# [an OPTIONAL QUEST]: Scaling of Topological Overlap Matrices to make them comparable across sets. Topological Overlap Matrices of different data sets may have different statistical properties. For example, the TOM in the male data may be systematically lower than the TOM in female data. Since consensus is defined as the component-wise minimum of the two TOMs, a bias may result. Here we illustrate a simple scaling that mitigates the effect of different statistical properties to some degree. We scale the remaining TOMs such that their 95th percentile equals the 95th percentile of the first TOM:

scaleP = 0.95 # Define the reference percentile
scaleQuant = rep(1, nSets) # These are TOM values at reference percentile
scalePowers = rep(1, nSets) # Scaling powers to equalize reference TOM values
scaledTOM = array(0, dim = c(nSets, nGenes, nGenes))
forScaling = list()

# Loop over sets
for (set in 1:nSets)
{
  forScaling[[set]] = as.dist(TOM[set, , ]) # Select the sampled TOM entries
  scaleQuant[[set]] = quantile(forScaling[[set]], probs = scaleP, type = 8) # Calculate the 95th percentile
  
  if(set==1) {scaledTOM[set, , ] = TOM[set, , ]}
  if(set>1) { scalePowers[[set]] = log(scaleQuant[1])/log(scaleQuant[set])
  scaledTOM[set, ,] = TOM[set, ,]^scalePowers[set]
  }
}

# The array scaledTOM now contains the scaled TOMs. To see what the scaling achieved, we form a quantile-quantile plot of the male and female topological overlaps before and after scaling:

# plot of the scaled and unscaled samples (set 1 values, scaled and unscaled, should fall on the isometric line):
par(mfrow=c(2,2))
for(set in 1:nSets) {
  qqUnscaled = qqplot(TOM[1, , ], TOM[set, , ], plot.it = T, col = "black", cex = 0.6, pch = 20, ylim=c(0, 0.6), xlim=c(0, 0.6))
  qqScaled = qqplot(scaledTOM[1, , ], scaledTOM[set, , ], plot.it = F)
  points(qqScaled$x, qqScaled$y, col = "red", cex = 0.6, pch = 20);
  abline(a=0, b=1, col = "blue")
  legend("topleft", legend = c("Unscaled TOM", "Scaled TOM"), pch = 20, col = c("black", "red")) }
# The result of scaling the sets is shown. Question to ask yourself: Did scaling bring the data closer to the isometric line (shown in blue), or otherwise adjust for some systematic deviation?
###############################################################################################################################################################################

rm(forScaling)
rm(TOM) # remove either TOM or scaledTOM (rm the one you don't use) to save space


### analyze TOM network properties
network_colors <-  c("#EEC900", "#CD6600", "#87CEFA", "#104E8B")

par(mfrow=c(4,5))

cutoffs <- c(0.05, 0.1, 0.15, 0.2, 0.4)
par(mar=c(4,4,4,4))
for(set in 1:nSets) {
  for(cut in cutoffs) {
    adj <-  scaledTOM[set, , ] # I chose to use the scaled TOM based on the qq-plots above.
    dimnames(adj) = list(colnames(multiPIC[[set]]$data), colnames(multiPIC[[set]]$data))
    diag(adj) <- NA
    adj[adj > cut] = 1  
    adj[adj != 1] = 0
    
    network <- graph.adjacency(adj)
    network <- igraph::simplify(network)  # removes self-loops
    network <- delete.vertices(network, degree(network)==0) # removes nodes without an edge
    
    plot(network, edge.arrow.size = 0, vertex.size=6, vertex.label=NA, vertex.color=network_colors[set], main=paste0(setLabels[set], ', cor^7 > ', cut))
    
  }}

for(set in 1:nSets) {
  for(cut in cutoffs) {
    adj <-  scaledTOM[set, , ]
    dimnames(adj) = list(colnames(multiPIC[[set]]$data), colnames(multiPIC[[set]]$data))
    diag(adj) <- NA
    adj[adj > cut] = 1  
    adj[adj != 1] = 0
    
    network <- graph.adjacency(adj)
    network <- igraph::simplify(network)  # removes self-loops
    network <- delete.vertices(network, degree(network)==0) # removes nodes without an edge
    par(mar=c(4,6,4,5))
    plot(degree.distribution(network, mode='in'), pch=19, log='xy', las=1, main=paste0(setLabels[set], ', cor^7 > ', cut), ylab='degree distribution') # should show a negative relationship that is linear on this log-log plot
  }}



par(mfrow=c(nSets,2))

cutoff = 0.15 # cutoff chosen based on the analysis above
# WGCNA will use its own criteria for module identification, but this cutoff can be used later for network visualization (ie. in cytoscape or igraph)

for(set in 1:nSets) {
  adj <-  scaledTOM[set, , ]
  dimnames(adj) = list(colnames(multiPIC[[set]]$data), colnames(multiPIC[[set]]$data))
  diag(adj) <- NA
  adj[adj > cutoff] = 1  
  adj[adj != 1] = 0
  
  network <- graph.adjacency(adj)
  network <- igraph::simplify(network)  # removes self-loops
  network <- delete.vertices(network, degree(network)==0) # removes nodes without an edge
  par(mar=c(3,4,4,2))
  plot(network, edge.arrow.size = 0, vertex.size=6, vertex.label=NA, vertex.color=network_colors[set], main=paste0(setLabels[set], ', cor^7 > ', cutoff))
  par(mar=c(6,6,6,6))
  plot(degree.distribution(network, mode='in'), pch=19, log='xy', las=1, ylab='degree distribution') # should show a negative relationship that is linear on this log-log plot
}

rm(adj)
rm(adjacencies)


#############################
# Identify Modules:
#############################
# loop though the sets to identify modules:
par(mfrow=c(2,2))
par(mar=c(6,6,6,6))

minModuleSize = 8 # set module size; I used the output of the code below to help choose the minimum mod size.  I used the heatmap with the modules shown on the margins to reveal how the modules related to the covariance among features.  I favored a number that was small enough to give clusters that were mosly non 'broken up', but a number that was not so small that what looked like separate clusters in the heatmap came out as modules in WGCNA

dynamicMods <- list()

for(set in 1:nSets) {
  mzTree <- hclust(as.dist(1-scaledTOM[set, , ]), method = "average")
  plot(mzTree, xlab="", sub="", main = setLabels[set], labels = FALSE, las=1, hang = 0.04) 
  
  # Module identification and assignment using dynamic tree cut:
  dynamicMods[[set]] = cutreeDynamic(dendro = mzTree, distM = 1-scaledTOM[set, , ], deepSplit = 4, pamRespectsDendro = T, minClusterSize = minModuleSize) }

dynamicMods[[set]] # module assignments for the metabolites (shown here as an example)

lapply(dynamicMods, table) # # module assignments for metabolites, class 0 are metabolites that weren't assigned to a module

dynamicColors = lapply(dynamicMods, labels2colors) # Convert numeric labels into colors
lapply(dynamicColors, table)

plotDendroAndColors(mzTree, dynamicColors[[set]], "Dynamic Modules", hang = 0.03, las=1, main=setLabels[set])


## find the total number of modules
tmp <- lapply(dynamicColors, table)
nMods <- list()
for(i in 1:nSets) {nMods[[i]] <- length(tmp[[i]])}
nMods <- max(unlist(nMods))


#WGCNA assigns module colors in order of module size 'turquoise' being the largest.  [NOTE]: 'grey' is reserved for all features that don't make it into a module, and grey may actually be larger that turquoise.

# replace module colors with letters, A,B,C.... depending order of module size (but keep 'grey' as 'grey')
module_renaming_vector <- c('grey', LETTERS[1:(nMods-1)]) # out of letters..
module_renaming_vector
added_letters <- LETTERS[1:table(is.na(module_renaming_vector))[2]]
module_renaming_vector[is.na(module_renaming_vector)] <- paste0(added_letters, added_letters)
names(module_renaming_vector) <- 0:(nMods-1)
module_renaming_vector

renamedMods <- dynamicMods
for(set in 1:nSets) {
  renamedMods[[set]] <- module_renaming_vector[dynamicMods[[set]]+1] }

renamedMods


mzTree.order <- list()

# use trees and/or heatmaps to evaluate the modules:
for(set in 1:nSets) {
  tmptree <- hclust(as.dist(1-scaledTOM[set, , ]), method='average')
  mzTree.order[[set]] <- tmptree$order
  plotDendroAndColors(tmptree, dynamicColors[[set]], "Dynamic Modules", hang = 0.03, las=1, main=setLabels[set]) 
}

for(set in 1:nSets) {
  colfunc <- colorRampPalette(c("grey88", mycol4[set]))  
  
  x <- abs(cor(multiPIC[[set]]$data[ ,mzTree.order[[set]]]))
  x[1:4,1:4]
  diag(x) = NA
  
  legend_names <- unique(unlist(renamedMods[set]))
  legend_names[legend_names=='grey'] <- NA
  legend_names <- sort(legend_names, na.last=T)
  legend_names[is.na(legend_names)] <- 'non-modular'
  
  legend_colors <- unique(unlist(dynamicColors[[set]]))
  legend_colors[legend_colors=='grey'] <- NA
  legend_colors <- sort(legend_colors, na.last=T)
  legend_colors[is.na(legend_colors)] <- 'grey'
  
  jpeg(file=paste(setLabels[set], "heatmap plus modules.jpg"), width = 800, height = 700, quality=100, res=100)
  heatmap.2(x, trace='none', dendrogram = 'none', main=setLabels[set], RowSideColors= dynamicColors[[set]][mzTree.order[[set]]], Rowv=F, Colv=F, col=colorRampPalette(c('white', mycol4[c(2,4)][set])), key.xlab='abs(cor)', margin=c(12,12), keysize=1, key.title = '', labRow = NULL, labCol = NULL) 
  legend("left", legend = legend_names, col = legend_colors, lwd = 7, bty='n', title='WGCNA modules', inset=c(-0.1,-0.01), xpd=T, cex=0.8)
  
  dev.off()}


## Look at the intersection between modules between groups:
l <- lapply(renamedMods, as.factor)
df <- data.frame(matrix(unlist(l), nrow=length(l[[1]]), byrow=F))
head(df)
colnames(df) <- setLabels
rownames(df) <- features
head(df)

for(k in 1:nSets){
  mods <- list()
  for(i in 1:length(unique(df[ ,setLabels[k]]))) {
    levs <- unique(df[ ,setLabels[k]])
    lev <- levs[i]
    mods[[i]] <- rownames(df)[df[ ,setLabels[k]] == lev] }
  names(mods) <- paste(setLabels[k], unique(df[ ,setLabels[k]]))
  ifelse(k==1, out <- mods, out <- c(out, mods)) }

out
setmap(Venn(out), element_clustering = T, set_clustering = T, method='ward.D') # the column clustering here looks like a rationale for renaming modules.

out2 <- out[!grepl('grey', names(out))] # remove 'grey' modules, it is composed of features that aren't assigned o a module and so will not be given another name.
setmap(Venn(out2), element_clustering = T, set_clustering = T, method='ward.D') # the column clustering here looks like a rationale for renaming modules.

# a network visulization of the overlap between modules in each group (age,sex):
a <- print(overlap_pairs(Venn(out2)))
b <- strsplit(names(a), split='\\...')
c <- strsplit(sapply(b, "[[", 1), split=' ')
c <- lapply(c, function(x) paste(x, collapse="_"))
d <- strsplit(sapply(b, "[[", 2), split=' ')
d <- lapply(d, function(x) paste(x, collapse="_"))

tmp <- data.frame(set1 = unlist(c), set2 = unlist(d), joint = unlist(lapply(a, length)))
head(tmp)
rownames(tmp) = NULL
tmp

write.table(tmp, 'correlation among the metabolome/untargeted_wgcna_module_overlap.txt', row.names=F, quote=F, col.names=F)

require(HEMDAG)
W <- weighted.adjacency.matrix(file='correlation among the metabolome/untargeted_wgcna_module_overlap.txt')

W[1:7,1:7]
table(W)
g <- graph_from_adjacency_matrix(W, mode='undirected', weighted = T)
g <- igraph::simplify(g)
plot(g, vertex.size=2, vertex.label=NA, edge.width=E(g)$weight)

E(g)$weight

node_sex <- 1:length(V(g))
vertex_attr(g)$name <- paste0('_', vertex_attr(g)$name)
node_sex[grepl('_female', vertex_attr(g)$name)] = 'female'
node_sex[grepl('_male', vertex_attr(g)$name)] = 'male'
node_sex <- as.factor(node_sex)
node_sex


plot(g, vertex.color = c(mycol4[2], mycol4[4])[node_sex], edge.width=E(g)$weight, vertex.label=NA, vertex.size=5)

g <- set_vertex_attr(g, 'sex', index = V(g), node_sex)

rm(a)
rm(b)
rm(c)
rm(d)

# follow steps to export to Cytoscape
require(RCy3)
# 1. Download the latest Cytoscape from http://www.cytoscape.org/download.php
# 2. Complete installation wizard
# 3. Launch Cytoscape
cytoscapePing() # confirm that this line gives: "You are connected to Cytoscape!"

createNetworkFromIgraph(g, collection = "My Igraph Network Collection") # this should open in Cytoscape automatically

# NOTE: there is extensive and complex overlap between the modules of different groups.  Rather than reassign group names based on similarlity (which will be misleading).

###############################

dev.off()
#############################################
# Calculate eigengenes (PC1 of each module)
eigenfeatures <- list()
modColors <- list()
trees <- list()

for(set in 1:nSets) {
  MEList = moduleEigengenes(as.data.frame(multiPIC[[set]]), colors = renamedMods[[set]])
  MEs = MEList$eigengenes
  MEDiss = 1-cor(MEs) # Calculate dissimilarity of module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average") # Cluster module eigengenes
  
  # clustering of eigenmetabolites (first PC of each module)
  plot(METree, main = "Clustering of module eigenmetabolites", xlab = "", sub = "")
  
  # [OPTIONAL] merging step: modules that are similar enough here can be merged
  MEDissThres = 0.20 # threshold for merging
  abline(h=MEDissThres, col = "red")
  
  # Call an automatic merging function
  merge = mergeCloseModules(as.data.frame(multiPIC[[set]]), dynamicColors[[set]], cutHeight = MEDissThres, verbose = 3) # using conventional WGCNA mod names (colors)
  merge.mods = mergeCloseModules(as.data.frame(multiPIC[[set]]), renamedMods[[set]], cutHeight = MEDissThres, verbose = 3) # using custom module names
  
  mergedColors = merge$colors # The merged module colors
  
  # visualize the effect of merging
  mzTree <- hclust(as.dist(1-scaledTOM[set, , ]), method = "average")
  plotDendroAndColors(mzTree, cbind(dynamicColors[[set]], mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), hang = 0.03, las=1)
  
  merge.mods$colors # The merged module names
  
  merged.modMEs = merge.mods$newMEs # Eigengenes of the new merged modules [NOTE here you are making the choice to use the merged modules or the previous unmerged modules]
  moduleColors =  merge.mods$colors # make new 'moduleColors' based on the merged modules
  
  merge.mods$newMEs 
  merged.modMEs
  
 # save network entities
  eigenfeatures[[set]] <- merged.modMEs
  modColors[[set]] <- moduleColors
  trees[[set]] <- mzTree
  
  # and on the disk:
  save(merged.modMEs, moduleColors, mzTree, file = paste0('correlation among the metabolome/', shortLabels[set], '_networkConstruction.RData'))
}


########################################################################
# Relate modules to traits and identifying important metabolites:
########################################################################

# load eigenmodules:
eigenfeatures <- list()
for(set in 1:nSets) {
  shortLabels[set]
  load(paste0('correlation among the metabolome/', shortLabels[set], '_networkConstruction.RData'))
  eigenfeatures[[set]] <- merged.modMEs
}

eigenfeatures


ls_mod <- list()

ls.by.mz <- function(x) {
  m <- sma(ls ~ x -1)
  return(c(m$r2, m$pval)) }

for(set in 1:nSets){
  tmpic <- apply(PICs[[set]], 2, unlist)  
  ls <- t(tmpic)[ ,'mean.ls']
  out <- matrix(nrow=ncol(eigenfeatures[[set]]), ncol=2)
  for(i in 1:ncol(eigenfeatures[[set]])) {
    eigenmod <- unlist(eigenfeatures[[set]][i])
    out[i, ] <- unlist(ls.by.mz(eigenmod)) }
  
  out <- as.data.frame(out)
  rownames(out) <- colnames(eigenfeatures[[set]])
  colnames(out) <- c('r.sq','pval')
  out$FDR <- p.adjust(out$pval, method='BH') # conservative, becuase we allow the grey module to be included
  out$group <- setLabels[set]
  out <- out[ ,c(4,1:3)]
  ls_mod[[set]] <- out
  
  write.table(out, file=paste('~/Google Drive/My Drive/Documents/Hoffman_Metabolite_Phylogeny_Paper/correlation among the metabolome/lifespanPIC_by_eigenmodules_major.axis.regression_untargeted', setLabels[set])) }


ls_eigenmod_table <- rbind(ls_mod[[1]], ls_mod[[2]])
ls_eigenmod_table$module <- rownames(ls_eigenmod_table)
head(ls_eigenmod_table)
ls_eigenmod_table <- ls_eigenmod_table [ ,c(1, 5, 2:4)]

write.table(ls_eigenmod_table, file='~/Google Drive/My Drive/Documents/Hoffman_Metabolite_Phylogeny_Paper/correlation among the metabolome/lifespanPIC_by_eigenmodules_major.axis.regression_untargeted_Table.csv', sep=',', row.names = F)

ls_eigenmod_table

##########################################################################
# the following choices are manually coded based on a module and trait of interest:
dev.off()
par(xpd=FALSE)

set_of_interest <- 'female'
trait_of_interest <- 'mean-ls'
module_of_interest <- 'K'
##########################################################################
setLabels == set_of_interest

set <- seq(1:nSets) [setLabels == set_of_interest]
y <- Traits[[set]]
x <- eigenfeatures[[set]][ ,paste0('ME', module_of_interest)]

plot(y ~ x, pch=16, las=1, main=set_of_interest,  xlab=paste('eigenmodule', module_of_interest), ylab=trait_of_interest)
abline(sma(y ~ x -1))
summary(s <- sma(y ~ x -1))
legend('bottomright', bty='n', legend=c('type II regression', paste('P=', round(unlist(s$pval), 4)), paste('r^2=', round(unlist(s$r2), 3))))


################################################################################################################
# use mummichog to search for pathways whose activity could explain the features in the modules
################################################################################################################

# to get LC-MS feature information: retention time and mass/charge ratio:

neg.retention.times <- read.table('~/Google Drive/My Drive/Documents/Hoffman_Metabolite_Phylogeny_Paper/Negative ion mode metabolite info.csv', sep=',' , stringsAsFactors = F, header=T) 
head(neg.retention.times)
neg.retention.times$feature <- paste0('neg_', neg.retention.times$Mass)
head(neg.retention.times)

pos.retention.times <- read.table('~/Google Drive/My Drive/Documents/Hoffman_Metabolite_Phylogeny_Paper/Positive ion mode metabolite info.csv', sep=',' , stringsAsFactors = F, header=T) 
head(pos.retention.times) 
pos.retention.times$feature <- paste0('pos_', pos.retention.times$Mass)
head(pos.retention.times) 

## look to see how many of the LC-MS features that have been analysed above are included in the feature meta data:
table(features %in% neg.retention.times$feature) 
table(features %in% pos.retention.times$feature)
features[!features %in% pos.retention.times$feature & !features %in% neg.retention.times$feature] # one of the 589 LC-MS features is not present in the feature meta data 

pos.retention.times <- pos.retention.times[pos.retention.times$feature %in% features, ] # cut down meta data to only those analyzed here (the ones measured in all strains)
neg.retention.times <- neg.retention.times[neg.retention.times$feature %in% features, ] # cut down meta data to only those analyzed here (the ones measured in all strains)

meta.data <- rbind(pos.retention.times, neg.retention.times)

names(modColors[[2]]) <- NULL
modColors

module_assignment <- data.frame(modColors[[1]], modColors[[2]])
colnames(module_assignment) <- paste0(setLabels, '_modules')
module_assignment$feature <- features
head(module_assignment)

head(meta.data)
x <- merge(module_assignment, meta.data, by='feature')
x$mode <- ifelse(1:length(x$feature) %in% grep('pos', x$feature), 'positive', 'negative')
head(x)
meta.data <- x
head(meta.data)
rm(x)


## add module-level association stats to the meta.data
ls_mod[[1]][1:4, ]
ls_mod[[1]][ ,5] <- rownames(ls_mod[[1]])
ls_mod[[1]][1:4, ]
ls_mod[[2]][ ,5] <- rownames(ls_mod[[2]])
ls_mod[[2]][1:4, ]
ls_module_assoc <- rbind(ls_mod[[1]], ls_mod[[2]])
ls_module_assoc <- cSplit(ls_module_assoc, 'V5', sep='ME')
ls_module_assoc  <- ls_module_assoc [ ,c(1:4, 6)]
colnames(ls_module_assoc)[5] <- 'module'
head(ls_module_assoc)

out_F <- merge(ls_module_assoc[ls_module_assoc$group == 'female', ], meta.data[ ,c(1,2,4,5,7)], by.x='module', by.y='female_modules')
head(out_F)
table(out_F$module, out_F$mode) 


out_M <- merge(ls_module_assoc[ls_module_assoc$group == 'male', ], meta.data[ ,c(1,3,4,5,7)], by.x='module', by.y='male_modules')
head(out_M)
table(out_M$module, out_M$mode) 

# columns required for mummichog input: c('mz', 'rtime', 'p-value', 't-score')

# make input file for the male data:
fem_neg_mummichog_input <- out_F[out_F$mode == 'negative' ,c('Mass', 'Retention.Time', 'pval', 'r.sq')]
colnames(fem_neg_mummichog_input) <- c('mz', 'rtime', 'p-value', 't-score')

write.table(fem_neg_mummichog_input, file='~/Google Drive/My Drive/Applications/mummichog-1.0.9/phylo_mz/fem_neg_mummichog_input.txt', sep='\t', row.names = F)

fem_pos_mummichog_input <- out_F[out_F$mode == 'positive' ,c('Mass', 'Retention.Time', 'pval', 'r.sq')]
colnames(fem_pos_mummichog_input) <- c('mz', 'rtime', 'p-value', 't-score')
head(fem_pos_mummichog_input)

write.table(fem_pos_mummichog_input, file='~/Google Drive/My Drive/Applications/mummichog-1.0.9/phylo_mz/fem_pos_mummichog_input.txt', sep='\t', row.names = F)


# make input file for the male data:
male_neg_mummichog_input <- out_M[out_M$mode == 'negative' ,c('Mass', 'Retention.Time', 'pval', 'r.sq')]
colnames(male_neg_mummichog_input) <- c('mz', 'rtime', 'p-value', 't-score')

write.table(male_neg_mummichog_input, file='~/Google Drive/My Drive/Applications/mummichog-1.0.9/phylo_mz/male_neg_mummichog_input.txt', sep='\t', row.names = F)

male_pos_mummichog_input <- out_M[out_M$mode == 'positive' ,c('Mass', 'Retention.Time', 'pval', 'r.sq')]
colnames(male_pos_mummichog_input) <- c('mz', 'rtime', 'p-value', 't-score')
head(male_pos_mummichog_input)

write.table(male_pos_mummichog_input, file='~/Google Drive/My Drive/Applications/mummichog-1.0.9/phylo_mz/male_pos_mummichog_input.txt', sep='\t', row.names = F)


# mummichog uses the 'p-value' comlumn to select the features of interest from among the list of features.  We may need to alter the p-value column to force mummichog to interpret the modules we want it to.  

#[Example]: I want to analyze enrichment among female module K
table(out_F$module, out_F$mode) 

head(out_F) 
z <- unique(out_F$pval[out_F$module == 'B'])
table(fem_neg_mummichog_input$`p-value` == z)

fem_negB <- fem_neg_mummichog_input
fem_negB$`p-value` <- ifelse(fem_negB$`p-value` == z, 0.0001, 1)
write.table(fem_negB, file='~/Google Drive/My Drive/Applications/mummichog-1.0.9/phylo_mz/fem_negB.txt', sep='\t', row.names = F)
