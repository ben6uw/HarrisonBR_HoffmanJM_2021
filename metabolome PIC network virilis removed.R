# Examine covariance among the metabolite PICs; i.e. what is the covariance of metabolties across the 50 million year old phylogeny?

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.11")
BiocManager::install("limma")
BiocManager::install('WGCNA')
BiocManager::install('org.Mm.eg.db')

library(WGCNA)
library(org.Mm.eg.db)
library(igraph)
library(pheatmap)
library(RColorBrewer)
library(gplots)
library(smatr)
library(FELLA)
library(org.Dm.eg.db)

library(plyr) # data wrangling
library(splitstackshape)
library(dplyr)
library(tidyr)
library(reshape2)
library(readr)

library(DMwR) # imputation
library(impute)

library(sjstats) # modeling tools
library(lme4)
library(lmerTest)
library(nlme) 
library(MuMIn) 
library(MASS)
library(smatr)
library(lmodel2)

library(graphics) # plotting
library(ggplot2) 
library(gplots) 
library(corrplot)
library(RColorBrewer)
library(limma)

library(ape) # phylogenetic analysis
library(phytools) 
library(geiger) 
library(caper) 
library(AssocTests)
library(igraph)

mycol4 = c("#EEC900", "#CD6600", "#87CEFA", "#104E8B")

setwd("~/Google Drive File Stream/My Drive/Documents/Hoffman_Metabolite_Phylogeny_Paper") 

flytree <- read.tree('fly.species.list.nwk')
plot(flytree)
axisPhylo(side = 1, root.time = NULL, backward = TRUE)

dat <- read.table('comparative.mz.targeted.data.log_only', header=T, sep='\t')
head(dat) # a leading X is added to the metabolite names that begin with a number

names(dat)[1:18]
dat[, 18:ncol(dat)] -> mets
mets[1:4,1:4]
mets <- as.matrix(mets)
mets <- mets[ ,colSums(is.na(mets)) == 0] # rid matrix of mzs with missing data
# remove metabolites identified by Danijel as artifacts
ghost_peaks <- c('4-Trimethylammoniobutanoate', 'X4.Trimethylammoniobutanoate', 'Sorbitol', 'Xanthosine', 'Oxypurinol', 'Uridine', 'Trimethylamine', 'Valine', 'Adenine', 'Pipecolic acid', 'Pipecolic.acid') # these were one member of a few pair of highly-correlated metabolties (r>0.98 in another analysis),Danijel suspected that these metabolties are actualy LC-MS artifacts
mets <- mets[ ,!colnames(mets) %in% ghost_peaks] 
centered_mets <- t(scale(t(mets), scale=F)) # center, but don't scale metabolites by sample

# a look at the data
par(mfrow=c(2,2))
boxplot(mets[order(apply(mets, 1, median)), ], use.cols = F, main='log only', xlab='sample', las=1)
boxplot(centered_mets[order(apply(mets, 1, median)), ], use.cols = F, main='log, mean-centered', xlab='sample', las=1)
boxplot(mets[ ,order(colMeans(mets))], use.cols = T, main='log only', xlab='metabolite', las=1)
boxplot(centered_mets[ ,order(colMeans(mets))], use.cols = T, main='log, mean-centered', xlab='metabolite', las=1)

centered_mets[1:4, 1:2]
# replace leading X in metabolite names:
colnames(centered_mets)[grepl('X', colnames(centered_mets))]
colnames(centered_mets) <- gsub('X1', '1', colnames(centered_mets))
colnames(centered_mets) <- gsub('X2', '2', colnames(centered_mets))
colnames(centered_mets) <- gsub('X3', '3', colnames(centered_mets))
colnames(centered_mets) <- gsub('X4', '4', colnames(centered_mets))
colnames(centered_mets) <- gsub('X5', '5', colnames(centered_mets))
colnames(centered_mets)[grepl('X', colnames(centered_mets))]

centered_mets[1:4, 1:2]
dat <- cbind(dat[ ,1:17], centered_mets)
dat[1:5, 1:20]


########################################################################################
# should we use residual lifespan from LS ~ fly mass ?
par(mfrow=c(2,2))
par(mar=c(5,4,4,2))
plot(median.ls ~ weight.mg, dat, pch=16, col=as.factor(dat$sex), las=1, cex.lab=1.5)
summary(sma(median.ls ~ weight.mg * sex, dat)) # slopes do not differ between the sexes
abline(sma(median.ls ~ weight.mg, dat[dat$sex=='F', ]), col=1)
abline(sma(median.ls ~ weight.mg, dat[dat$sex=='M', ]), col=2)
abline(sma(median.ls ~ weight.mg, dat[dat$sex=='F' & dat$species != 'virilis', ]), col=1, lty=2)
abline(sma(median.ls ~ weight.mg, dat[dat$sex=='M' & dat$species != 'virilis', ]), col=2, lty=2)
legend('topleft', legend=c('female', 'male', 'dash = w/o virilis'), pch=16, col=c(1:2, 0), bty='n')

plot(mean.ls ~ weight.mg, dat, pch=16, col=as.factor(dat$sex), las=1, cex.lab=1.5)
summary(sma(mean.ls ~ weight.mg * sex, dat)) # slopes do not differ between the sexes
abline(sma(mean.ls ~ weight.mg, dat[dat$sex=='F', ]), col=1)
abline(sma(mean.ls ~ weight.mg, dat[dat$sex=='M', ]), col=2)
abline(sma(mean.ls ~ weight.mg, dat[dat$sex=='F' & dat$species != 'virilis', ]), col=1, lty=2)
abline(sma(mean.ls ~ weight.mg, dat[dat$sex=='M' & dat$species != 'virilis', ]), col=2, lty=2)
legend('topleft', legend=c('female', 'male', 'dash = w/o virilis'), pch=16, col=c(1:2, 0), bty='n')

plot(max.ls ~ weight.mg, dat, pch=16, col=as.factor(dat$sex), las=1, cex.lab=1.5)
summary(sma(max.ls ~ weight.mg * sex, dat)) # slopes do not differ between the sexes
abline(sma(max.ls ~ weight.mg, dat[dat$sex=='F', ]), col=1)
abline(sma(max.ls ~ weight.mg, dat[dat$sex=='M', ]), col=2)
abline(sma(max.ls ~ weight.mg, dat[dat$sex=='F' & dat$species != 'virilis', ]), col=1, lty=2)
abline(sma(max.ls ~ weight.mg, dat[dat$sex=='M' & dat$species != 'virilis', ]), col=2, lty=2)
legend('topleft', legend=c('female', 'male', 'dash = w/o virilis'), pch=16, col=c(1:2, 0), bty='n')

# a positive LS ~ mass relationship is driven mostly by D.virilis, and removing virilis actually gives a negative relationship.  Seems there is no good way to 'remove the effect' of LS ~ mass. 

# perhaps a phylogenetic model will reduce this effect, either way, a later analysis with and without virilis seems appropriate, along with a mention of the confounding (and complex) LS ~ mass relationship. 
###################################################

###################################################
# repeat this analysis without virilis.  Analyze the effect of virilis removal on the ls ~ mz relationships:
###################################################

no.vir <- dat[dat$species != 'virilis', ] # remove virilis from data

## PIC including intraspecies covariance
traits <- names(no.vir)[c(10:13, 18:ncol(no.vir))] 
groups <- paste(no.vir$sex, no.vir$age.days)
group <- names(table(groups))

for(k in group) {
  tmp.no.vir <- no.vir[groups == k, ]
  no.species <- length(table(tmp.no.vir$species))
  
  tmpics <- matrix(ncol=no.species-1, nrow=length(traits))
  
  for(i in 1:length(traits)) {
    trait <- traits[i]
    tmp <- tmp.no.vir[ ,trait]
    names(tmp) <- tmp.no.vir$species
    tmp <- tmp[!is.na(tmp)]
    tmp <- aggregate(tmp, by=list(names(tmp)), list)
    names(tmp$x) <- tmp$Group.1
    tmp  <- tmp$x
    tmp.tree <- keep.tip(flytree, names(tmp))
    tmp  <- tmp[tmp.tree$tip.label]
    tmpics[i, ] <- pic.ortho(tmp, tmp.tree, intra = T, var.contrasts = F) }
  # method from Felsenstein (2008) Am. Nat.
  rownames(tmpics) <- traits
  write.table(tmpics, file=paste('~/Google Drive File Stream/My Drive/Documents/Hoffman_Metabolite_Phylogeny_Paper/correlation among the metabolome/pics', k, 'virilis removed'))
}

nv.picM5 <- read.table('correlation among the metabolome/pics M 5 virilis removed')
nv.picM31 <- read.table('correlation among the metabolome/pics M 31 virilis removed')
nv.picF5 <- read.table('correlation among the metabolome/pics F 5 virilis removed')
nv.picF31 <- read.table('correlation among the metabolome/pics F 31 virilis removed')

##############################################################################

nSets = 4
nv.PICs <- list(nv.picF5, nv.picF31, nv.picM5, nv.picM31)
nv.setLabels = c("young female - virilis removed", "old female - virilis removed", 'young male - virilis removed', 'old male - virilis removed')
nv.shortLabels = c('nv.F5', 'nv.F31' ,'nv.M5', 'nv.M31')

# plot correlation matricies of metabolites for each sex x age:
F5_pal <- colorRampPalette(c("grey", "white", "blue"))(n = 50)
F31_pal <- colorRampPalette(c("grey", "white", "purple"))(n = 50)
M5_pal <- colorRampPalette(c("grey", "white", "brown"))(n = 50)
M31_pal <- colorRampPalette(c("grey", "white", "darkgoldenrod"))(n = 50)
setPals <- list(F5_pal, F31_pal, M5_pal, M31_pal)

for(set in 1:nSets) {
  tmpic <- apply(nv.PICs[[set]], 2, unlist)
  pic.cor <- cor(t(tmpic), method='pearson')
  heatmap.2(pic.cor, trace='n', col=unlist(setPals[set]), key.title = 'cor PIC-PIC', margins = c(10,12), key.xlab = 'spearman rho', cexRow=0.5, cexCol=0.5, keysize=1, dendrogram='none', main=paste('PIC correlation', nv.setLabels[set])) }

##############################################################################
# compare ls PIC ~ metabolite PIC from PICs with and without virilis

pic.cors <- list()
nv.pic.cors <- list()

for(set in 1:nSets){
  tmpic <- apply(nv.PICs[[set]], 2, unlist)
  nv.pic.cors[[set]] <- cor(t(tmpic), method='pearson')
  
  tmpic <- apply(PICs[[set]], 2, unlist)
  pic.cors[[set]] <- cor(t(tmpic), method='pearson') }

traits.of.interest <- list('mean.ls', 'median.ls', 'max.ls')

par(mfrow=c(4,3))

for(set in 1:nSets){
  for(i in 1:3){
plot(pic.cors[[set]][ ,traits.of.interest[[i]]], nv.pic.cors[[set]][ ,traits.of.interest[[i]]], main=paste(setLabels[set], '-', traits.of.interest[[i]]), col=mycol4[set], pch=16, ylab=paste('virilis removed', '(cor', traits.of.interest[[i]], '~ mzPIC)'), xlab=paste('all species', '(cor', traits.of.interest[[i]], '~ mzPIC)')) 
    legend('bottomright', legend=paste('cor=', round(cor(pic.cors[[set]][ ,traits.of.interest[[i]]], nv.pic.cors[[set]][ ,traits.of.interest[[i]]]), 3)), bty='n')} }

out <- matrix(nrow=4, ncol=3)
for(set in 1:nSets){
  for(i in 1:3){
    out[set, i] <- cor(pic.cors[[set]][ ,traits.of.interest[[i]]], nv.pic.cors[[set]][ ,traits.of.interest[[i]]]) }}

out
colMeans(out)

# of the three lifespan traits, mean lifespan seems to have the smallest effect of virilis on the relationship between PIC(ls) ~ PIC(mz).  I say this because, the correlation betweeen cor(ls ~ mz) for the with and without virilis PICs are themselves corrleated with at least cor=0.486 (for the young female metabolome) with an average cor across groups of 0.66 [also, see the plots generated in the loops above].

# I therefore choose to analyze mean lifespan and to note that the ls ~ mz relationships are note broadly effected by removing virilis.
#######################################################################################################################################################################################

#############################################################
## use WGCNA to pick modules
#############################################################

# add KEGG and other database names to metabolites [these will be useful later for pathway analysis with FELLA]: 
setwd("~/Google Drive File Stream/My Drive/Documents/Hoffman_Metabolite_Phylogeny_Paper/correlation among the metabolome")
options(stringsAsFactors = FALSE)
mz.info <- read.table('mz.info.txt', header=T, sep=',')
head(mz.info)

rownames(nv.picM5)[!rownames(nv.picM5) %in% mz.info$metabolite]

rownames(nv.picM5) <- gsub('.ac', ' ac', rownames(nv.picM5), fixed=T)
rownames(nv.picM5) <- gsub('.Acid', ' acid', rownames(nv.picM5), fixed=T)
rownames(nv.picM5) <- gsub('.', '-', rownames(nv.picM5), fixed=T)
rownames(nv.picM5) <- gsub('D-Glyceraldehyde-3-phosphate', 'D-Glyceraldehyde 3-phosphate', rownames(nv.picM5), fixed=T)
rownames(nv.picM5) <- gsub('gamma-Aminobutyric acid', 'gamma-Aminobutyric Acid', rownames(nv.picM5), fixed=T)
rownames(nv.picM5) <- gsub('Glucose-6-phosphate', 'Glucose 6-phosphate', rownames(nv.picM5), fixed=T)
rownames(nv.picM5) <- gsub('Glycerol-3-phosphate', 'Glycerol 3-phosphate', rownames(nv.picM5), fixed=T)
rownames(nv.picM5) <- gsub('Oxidized-glutathione', 'Oxidized glutathione', rownames(nv.picM5), fixed=T)

rownames(nv.picM5)[!rownames(nv.picM5) %in% mz.info$metabolite]
mz.info$metabolite[!mz.info$metabolite %in% rownames(nv.picM5)]

mzs <- rownames(nv.picM5)

rownames(nv.picF5) <- mzs 
rownames(nv.picF31) <- mzs 
rownames(nv.picM31) <- mzs 


# Form a multi-set structure that will hold the metabolite data and the organismal phenotypes (i.e. lifespan and weight).
nv.multiPIC = vector(mode = "list", length = nSets) # empty entity

nv.multiPIC[[1]] = list(data = as.data.frame(t(nv.picF5[-c(1:4), ]))) # remove the trait data (first 4 rows) and transpose
names(nv.multiPIC[[1]]$data) = rownames(nv.picF5)[5:nrow(nv.picF5)]
nv.multiPIC[[2]] = list(data = as.data.frame(t(nv.picF31[-c(1:4), ]))) # remove the trait data (first 4 rows) and transpose
names(nv.multiPIC[[2]]$data) = rownames(nv.picF31)[5:nrow(nv.picF31)]
nv.multiPIC[[3]] = list(data = as.data.frame(t(nv.picM5[-c(1:4), ]))) # remove the trait data (first 4 rows) and transpose
names(nv.multiPIC[[3]]$data) = rownames(nv.picM5)[5:nrow(nv.picM5)]
nv.multiPIC[[4]] = list(data = as.data.frame(t(nv.picM31[-c(1:4), ]))) # remove the trait data (first 4 rows) and transpose
names(nv.multiPIC[[4]]$data) = rownames(nv.picM31)[5:nrow(nv.picM31)]

nv.Traits = vector(mode="list", length = nSets)

nv.Traits[[1]] <- t(nv.picF5[1:4, ])
nv.Traits[[2]] <- t(nv.picF31[1:4, ])
nv.Traits[[3]] <- t(nv.picM5[1:4, ])
nv.Traits[[4]] <- t(nv.picM31[1:4, ])

collectGarbage()

exprSize = checkSets(nv.multiPIC) # Check that the data has the correct format for many functions operating on multiple sets (I think this will give an error if there is a problem)

sampleTrees = list()
for (set in 1:nSets)
{sampleTrees[[set]] = hclust(dist(nv.multiPIC[[set]]$data), method = "average")}

par(mfrow=c(2,2))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets)
  plot(sampleTrees[[set]], main = paste("Sample clustering in", nv.setLabels[set]), xlab="", sub="", cex = 0.7, las=1)

collectGarbage()

# Check the size of the leftover data
exprSize = checkSets(nv.multiPIC)

# Define data set dimensions
nGenes = exprSize$nGenes
nSamples = exprSize$nSamples

save(nv.multiPIC, nv.Traits, nGenes, nSamples, nv.setLabels, nv.shortLabels, exprSize, file = "virilis_removed metabolomeConsensus-dataInput.RData") # saves the data in a particular format

nSets = checkSets(nv.multiPIC)$nSets

# Evaluate a set of soft-thresholding powers
#################################################################
# Constructing a weighted gene network entails the choice of the soft thresholding power β to which co-expression similarity is raised to calculate adjacency [3]. The authors of [3] have proposed to choose the soft thresholding power based on the criterion of approximate scale-free topology. We refer the reader to that work for more details; here we illustrate the use of the function pickSoftThreshold that performs the analysis of network topology and aids the user in choosing a proper soft-thresholding power. The user chooses a set of candidate powers (the function provides suitable default values), and the function returns a set of network indices that should be inspected.

# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets)

powers = c(seq(4,10,by=1), seq(12,20, by=2)) # make a set of power levels to calculate over

# Call the network topology analysis function for each set in turn
for (set in 1:nSets)
  powerTables[[set]] = list(data = pickSoftThreshold(nv.multiPIC[[set]]$data, powerVector=powers, verbose = 2)[[2]])

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
  { legend("bottomright", legend = nv.setLabels, col = colors, pch = 20) 
  } else
    legend("topright", legend = nv.setLabels, col = colors, pch = 20) 
}


par(mfrow = c(1,1))
for (set in 1:nSets) {
  if (set==1) {
    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],  xlab="Soft Threshold (cor^x)", type="b", ylim=c(-0.3,1), col=set, pch=16, las=1, xaxp  = c(0, 20, 20))
  } else 
    lines(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],type="b", col=set, pch=16) 
}
legend('topleft', legend=nv.setLabels, col=1:nSets, pch=16)


#####################################
softPower = 10 # chosen power 
#####################################

# how to check the scale-freeness of the network?  Scale-free = linear relationship between: frequency of log(conectivity) ~ log(conectivity), a.k.a a 'power law' relationship.  Do this check once the network is made:

# Network construction starts by calculating the adjacencies in the individual sets, using the soft thresholding power
# build an empty array to hold the adjacencies
adjacencies = array(0, dim = c(nSets, nGenes, nGenes));

# Calculate adjacencies in each individual data set
for (set in 1:nSets)
  adjacencies[set, , ] = abs(cor(nv.multiPIC[[set]]$data, use = "p"))^softPower # calculates the adjacency matrix, in this case the pearson correlation raised to the 'softpower'.  I don't know the purpose of the absolute value [abs()] part of this line, since x ^ power is positive.

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

### analyze TOM network properties
network_colors <-  c("#EEC900", "#CD6600", "#87CEFA", "#104E8B")

par(mfrow=c(4,5))

cutoffs <- c(0.05, 0.1, 0.15, 0.2, 0.4)
par(mar=c(4,4,4,4))
for(set in 1:nSets) {
  for(cut in cutoffs) {
    adj <-  scaledTOM[set, , ]
    dimnames(adj) = list(colnames(nv.multiPIC[[set]]$data), colnames(nv.multiPIC[[set]]$data))
    diag(adj) <- NA
    adj[adj > cut] = 1  
    adj[adj != 1] = 0
    
    network <- graph.adjacency(adj)
    network <- igraph::simplify(network)  # removes self-loops
    network <- delete.vertices(network, degree(network)==0) # removes nodes without an edge
    
    plot(network, edge.arrow.size = 0, vertex.size=6, vertex.label=NA, vertex.color=network_colors[set], main=paste0(nv.setLabels[set], ', cor^7 > ', cut))
    
  }}

for(set in 1:nSets) {
  for(cut in cutoffs) {
    adj <-  scaledTOM[set, , ]
    dimnames(adj) = list(colnames(nv.multiPIC[[set]]$data), colnames(nv.multiPIC[[set]]$data))
    diag(adj) <- NA
    adj[adj > cut] = 1  
    adj[adj != 1] = 0
    
    network <- graph.adjacency(adj)
    network <- igraph::simplify(network)  # removes self-loops
    network <- delete.vertices(network, degree(network)==0) # removes nodes without an edge
    par(mar=c(4,6,4,5))
    plot(degree.distribution(network, mode='in'), pch=19, log='xy', las=1, main=paste0(nv.setLabels[set], ', cor^7 > ', cut), ylab='degree distribution') # should show a negative relationship that is linear on this log-log plot
  }}



par(mfrow=c(4,2))

cutoff = 0.1 # cutoff chosen based on the analysis above
# WGCNA will use its own criteria for module identification, but this cutoff can be used later for network visualization (ie. in cytoscape or igraph)

for(set in 1:nSets) {
  adj <-  scaledTOM[set, , ]
  dimnames(adj) = list(colnames(nv.multiPIC[[set]]$data), colnames(nv.multiPIC[[set]]$data))
  diag(adj) <- NA
  adj[adj > cutoff] = 1  
  adj[adj != 1] = 0
  
  network <- graph.adjacency(adj)
  network <- igraph::simplify(network)  # removes self-loops
  network <- delete.vertices(network, degree(network)==0) # removes nodes without an edge
  par(mar=c(3,4,4,2))
  plot(network, edge.arrow.size = 0, vertex.size=6, vertex.label=NA, vertex.color=network_colors[set], main=paste0(nv.setLabels[set], ', cor^7 > ', cutoff))
  par(mar=c(6,6,6,6))
  plot(degree.distribution(network, mode='in'), pch=19, log='xy', las=1, ylab='degree distribution') # should show a negative relationship that is linear on this log-log plot
}


#############################
# Identify Modules:
#############################
# loop though the sets to identify modules:
dynamicMods <- list()
par(mfrow=c(2,2))
par(mar=c(6,6,6,6))

minModuleSize = 8 # set module size; I used the output of the code below to help choose the minimum mod size.  I used the heatmap with the modules shown on the margins to reveal how the modules related to the covariance among mzs.  I favored a number that was small enough to give clusters that were mosly non 'broken up', but a number that was not so small that what looked like separate clusters in the heatmap came out as modules in WGCNA

for(set in 1:nSets) {
  mzTree <- hclust(as.dist(1-scaledTOM[set, , ]), method = "average")
  plot(mzTree, xlab="", sub="", main = nv.setLabels[set], labels = FALSE, las=1, hang = 0.04) 
  
  # Module identification and assignment using dynamic tree cut:
  dynamicMods[[set]] = cutreeDynamic(dendro = mzTree, 
                                     distM = 1-scaledTOM[set, , ],
                                     deepSplit = 2, 
                                     pamRespectsDendro = T,
                                     minClusterSize = minModuleSize) }

dynamicMods[[set]] # module assignments for the metabolites (shown here as an example)

lapply(dynamicMods, table) # # module assignments for metabolites, class 0 are metabolites that weren't assigned to a module
dynamicColors = lapply(dynamicMods, labels2colors) # Convert numeric labels into colors
lapply(dynamicColors, table)

mzTree.order <- list()

# use trees and/or heatmaps to evaluate the modules:
for(set in 1:nSets) {
  tmptree <- hclust(as.dist(1-scaledTOM[set, , ]), method='average')
  mzTree.order[[set]] <- tmptree$order
  plotDendroAndColors(tmptree, dynamicColors[[set]], "Dynamic Modules", hang = 0.03, las=1, main=nv.setLabels[set]) 
}

for(set in 1:nSets) {
  colfunc <- colorRampPalette(c("grey88", mycol4[set]))  
  
  x <- abs(cor(nv.multiPIC[[set]]$data[ ,mzTree.order[[set]]]))
  x[1:4,1:4]
  diag(x) = NA
  
  legend_names <- unique(unlist(dynamicColors[set]))
  legend_names[legend_names=='grey'] <- NA
  legend_names <- sort(legend_names, na.last=T)
  legend_names[is.na(legend_names)] <- 'non-modular'
  
  legend_colors <- unique(unlist(dynamicColors[set]))
  legend_colors[legend_colors=='grey'] <- NA
  legend_colors <- sort(legend_colors, na.last=T)
  legend_colors[is.na(legend_colors)] <- 'grey'
  
  jpeg(file=paste(nv.setLabels[set], "heatmap plus modules virilis removed.jpg"), width = 800, height = 700, quality=100, res=100)
  heatmap.2(x, trace='none', dendrogram = 'none', main=nv.setLabels[set], RowSideColors= dynamicColors[[set]][mzTree.order[[set]]], Rowv=F, Colv=F, col=colfunc, key.xlab='abs(cor)', margin=c(12,12), keysize=1, key.title = '') 
  legend("left", legend = legend_names, col = legend_colors, lwd = 7, bty='n', title='WGCNA modules', inset=c(-0.1,-0.01), xpd=T, cex=0.8)
  
  dev.off()}

###############################

###############################
# Calculate eigengenes (PC1 of each module)
eigenMZs <- list()
modColors <- list()
modLabels <- list()
trees <- list()

for(set in 1:nSets) {
  MEList = moduleEigengenes(as.data.frame(nv.multiPIC[[set]]), colors = dynamicColors[[set]])
  MEs = MEList$eigengenes
  MEDiss = 1-cor(MEs) # Calculate dissimilarity of module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average") # Cluster module eigengenes
  
  # clustering of eigenmetabolites (first PC of each module)
  plot(METree, main = "Clustering of module eigenmetabolites", xlab = "", sub = "")
  
  # modules that are similar enough here can be merged
  MEDissThres = 0.25 # threshold for merging
  # Plot the cut line into the dendrogram
  abline(h=MEDissThres, col = "red")
  # Call an automatic merging function
  merge = mergeCloseModules(as.data.frame(nv.multiPIC[[set]]), dynamicColors[[set]], cutHeight = MEDissThres, verbose = 3)
  
  mergedColors = merge$colors # The merged module colors
  mergedMEs = merge$newMEs # Eigengenes of the new merged modules:
  
  mzTree <- hclust(as.dist(1-scaledTOM[set, , ]), method = "average")
  plotDendroAndColors(mzTree, cbind(dynamicColors[[set]], mergedColors),
                      c("Dynamic Tree Cut", "Merged dynamic"), hang = 0.03, las=1)
  
  moduleColors = mergedColors # Rename to moduleColors
  colorOrder = c("grey", standardColors(50)) # Construct numerical labels corresponding to the colors
  moduleLabels = match(moduleColors, colorOrder)-1
  MEs = mergedMEs
  
  # save network entities
  eigenMZs[[set]] <- MEs
  modColors[[set]] <- moduleColors
  modLabels[[set]] <- moduleLabels
  trees[[set]] <- mzTree
  
  # and on the disk:
  save(MEs, moduleLabels, moduleColors, mzTree, file = paste0(nv.setLabels[set], 'virilis removed_networkConstruction.RData'))
}

########################################################################
# Relate modules to traits and identifying important metabolites:
########################################################################

# Define numbers of genes and samples
nMZ = ncol(nv.multiPIC[[1]]$data)
nPIC.nodes <- numeric() # will define later

# Define/examine the organismal traits
nv.Traits

# Recalculate MEs with color labels
eigenMZs[[set]] 
modColors[[set]] 
modLabels[[set]]
trees[[set]] 

MEs <- list()

modmembcols <- colorRampPalette(c("blue", "white", "red"))(n = 50)

for(set in 1:nSets) {
  nPIC.nodes[set] <- nrow(nv.multiPIC[[set]]$data)
  MEs0 = moduleEigengenes(as.data.frame(nv.multiPIC[[set]]), modColors[[set]])$eigengenes
  MEs[[set]] = orderMEs(MEs0)
  
  moduleTraitCor = cor(MEs[[set]], nv.Traits[[set]], method='pearson')
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nPIC.nodes[set])
  
  # Will display correlations and their p-values
  textMatrix =  paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(7, 12, 3, 3))
  par(mfrow=c(1,1))
  # Display the correlation values within a heatmap plot
  
  jpeg(file=paste(nv.setLabels[set], " virilis removed mod-trait labeled heatmap.jpg"), width = 700, height = 600, quality=100, res=100)
  
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(nv.Traits[[set]]),
                 yLabels = names(MEs[[set]]),
                 ySymbols = names(MEs[[set]]),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = T,
                 cex.text = 1,
                 cex.lab=1,
                 zlim = c(-1,1),
                 main = paste(nv.setLabels[set], "\ncorrelation [eigenmetabolite - trait] (P value)"))
  
  dev.off() # the plot won't save unless you execute this line 
}

par(mar = c(5,4,4,2)) # reset normal-ish margins 


##########################################################################
# the following choices are manually coded
par(xpd=FALSE)

set_of_interest <- 'young male - virilis removed'
trait_of_interest <- 'mean-ls'
module_of_interest <- 'blue'

##########################################################################

set <- seq(1:4) [nv.setLabels == set_of_interest]
y <- nv.Traits[[set]][ ,trait_of_interest]
x <- MEs[[set]][ ,paste0('ME', module_of_interest)]

plot(y ~ x, pch=16, las=1, main=set_of_interest,  xlab=paste('eigenmodule', module_of_interest), ylab=trait_of_interest, col=module_of_interest)
abline(sma(y ~ x -1))
summary(s <- sma(y ~ x -1))
legend('bottomright', bty='n', legend=c('type II regression', paste('P=', round(unlist(s$pval), 4)), paste('r^2=', round(unlist(s$r2), 3))))

########################################################################
# network visualization using cytoscape
########################################################################

for(set in 1:nSets) {
  
  cyt.all = exportNetworkToCytoscape(scaledTOM[set, , ], 
                                     edgeFile = paste('Cytoscape', nv.setLabels[set], "edges.txt"), 
                                     nodeFile = paste('Cytoscape', nv.setLabels[set], "nodes.txt"),
                                     weighted = TRUE,
                                     threshold = 0.15, # determines extent of connectivity, this can be chosen based on the plots of degree distribution (above)
                                     nodeNames = mzs[5:101],
                                     altNodeNames = mzs[5:101],
                                     nodeAttr = modColors[[set]])
}

cyt.all



######################################################################################################################
# use FELLA, tutorial at: http://127.0.0.1:11793/library/FELLA/doc/musmusculus.pdf (modify code for Dme)
######################################################################################################################

set.seed(1)

graph <- buildGraphFromKEGGREST(organism = "dme", filter.path = c("01100", "01200", "01210", "01212", "01230"))  # Build network and filter out the large overview pathways:

tmpdir <- "~/Google Drive File Stream/My Drive/Documents/Hoffman_Metabolite_Phylogeny_Paper/FELLA_database"

unlink(tmpdir, recursive = TRUE) # this line allows you to rewrite the: buildataFromGraph()

buildDataFromGraph(graph, databaseDir = tmpdir, internalDir = FALSE, matrices = 'diffusion', normality = "diffusion", niter = 1000) # using matrices=listMethods() here leads to three matrices being made: 'diffusion', ' hypergemetric' and 'pagerank'.  It takes longer than building with a single method.  Note that 'normality' can either be 'diffusion' or 'pagerank', and this will determine which of those two analyses can run based on these data.  I ran normality as 'diffusion' and this caused an error saying that Z-scores would not be available for the 'pagerank' method.  It did however make a 'pagerank' matrix, so I suspect that I can still do non-parametric analysis by pagerank, rather than a parametric analysis that might depend on the outcome of this line. [NOTE]: I ultimately ran simulations (permutations) rather than paramertic analysis, so I don' think the 'normality' constrution was later used.

###[add back the three lines from: "# I'm not sure if/when these three lines are employed in the analysis" (above) if this doens't work]

fella.data <- loadKEGGdata(databaseDir = tmpdir, internalDir = FALSE, loadMatrix = 'diffusion') # using listMethods() here loads all 3 matrices, which in turn allows all 3 analysis methods

getInfo(fella.data) # prints the version KEGG used (Sergio says it always uses the most recent release)
fella.data

id.cpd <- getCom(fella.data, level = 5, format = "id") %>% names # 3968 KEGG compounds for Dmel
id.rx <- getCom(fella.data, level = 4, format = "id") %>% names # 748 KEGG enzymes for Dmel
id.ec <- getCom(fella.data, level = 3, format = "id") %>% names # 5417 KEGG reactions for Dmel

# get compounds from WGCNA modules:
set
nv.setLabels
set = 3

mzs = colnames(nv.multiPIC[[set]]$data) 

MOI = "turquoise" # Select module
ModforFella = is.finite(match(modColors[[set]], MOI))
ModforFella
ModforFella<- mzs[ModforFella]

# get KEGG compound names for metabolites:
mz.info <- read.table('~/Google Drive File Stream/My Drive/Documents/Hoffman_Metabolite_Phylogeny_Paper/correlation among the metabolome/mz.info.txt', header=T, sep=',')
head(mz.info)
rownames(mz.info) <- mz.info$metabolite

module_compounds <- mz.info[ModforFella, ]$KEGG
background_compounds <-  mz.info[mzs, ]$KEGG

c <- defineCompounds(module_compounds, background_compounds, data=fella.data) # prior to running many simulations (permutations), use this line to check that the compounds are on the KEGG graph (not all are usually on the graph, but this line along with the next will tell you if you might have incorrect names)
getExcluded(c) # show compounds that aren't ont he KEGG graph

module_analysis <- enrich(compounds = module_compounds, compoundsBackground = background_compounds, data = fella.data, method = 'diffusion', approx = "simulation", niter=10000)
getExcluded(module_analysis) # assoc. compounds were not analyzed because they were not on the KEGG graph, sould be the same as the 

## Results Table and Results Plot,
module_table <- generateResultsTable(object=module_analysis, data=fella.data)

path <- paste0(nv.setLabels[set], '_', MOI, '_module_table.txt')

write.table(module_table, file=path, sep='\t', row.names = F, quote=F)

table(module_compounds %in% module_table$KEGG.id) # are any of the input compound on the diffusion network nodes?  This is not required, but is a curiosity
table(background_compounds[!background_compounds%in% module_compounds] %in% module_table$KEGG.id) # are any of the background metabolites on the diffusion network nodes?

par(mfrow=c(1,1))
par(mar=c(4,4,4,4))

networkName <- paste(nv.setLabels[set], '-',MOI, 'module')

plot(module_analysis, method = "diffusion", main=networkName, threshold=0.15, data = fella.data, plotLegend = T)

############################################################
# export igraph to Cytoscape
module_graph <- generateResultsGraph(object = module_analysis, method = "diffusion", threshold = 0.15, data = fella.data) # make igraph object
plot(module_graph, edge.arrow.size = 0, vertex.size=6, vertex.label=NA, vertex.color=3)

# follow steps to export to Cytoscape
require(RCy3)
# 1. Download the latest Cytoscape from http://www.cytoscape.org/download.php
# 2. Complete installation wizard
# 3. Launch Cytoscape
cytoscapePing() # confirm that this line gives: "You are connected to Cytoscape!"

createNetworkFromIgraph(module_graph, title = networkName, collection = "My Igraph Network Collection") # this should open in Cytoscape automatically




