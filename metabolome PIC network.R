# Examine covariance among the metabolite PICs; i.e. what is the covariance of metabolites across the 50 million year old phylogeny?

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.11")
BiocManager::install("limma")
BiocManager::install('WGCNA')
BiocManager::install('HEMDAG')
BiocManager::install("diffcoexp")
BiocManager::install('ComplexHeatmap')

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
library(org.Dm.eg.db)
library(FELLA)
library(HEMDAG)
library(diffcoexp)
library(clusterRepro)

library(graphics) # plotting
library(ggplot2) 
library(gplots) 
library(corrplot)
library(RColorBrewer)
library(limma)
library(igraph)
library(pheatmap)
library(ComplexHeatmap)
library(gridExtra)

library(ape) # phylogenetic analysis
library(phytools) 
library(geiger) 
library(caper) 

dev.off()
gc()
set.seed(1)

# build colors
mycol7 <- brewer.pal(7, 'RdYlBu')
mycol4 = c("#EEC900", "#CD6600", "#87CEFA", "#104E8B")
reorderedmycol4 <- mycol4[c(2,1,4,3)]
mytranscol4 = c(rgb(238, 201, 0, 100, maxColorValue = 255), rgb(205, 102, 0, 100, maxColorValue = 255), rgb(135, 206, 250, 100, maxColorValue = 255), rgb(16, 78, 139, 100, maxColorValue = 255))


setwd("~/Google Drive/My Drive/Documents/Hoffman_Metabolite_Phylogeny_Paper/Harrison_Hoffman_GitHub") 

flytree <- read.tree('fly.species.list.nwk')
plot(flytree)
axisPhylo(side = 1, root.time = NULL, backward = TRUE)

dat <- read.table('comparative.mz.targeted.data.log_only', header=T, sep='\t')
head(dat) # a leading X is added to the metabolite names that begin with a number

names(dat)[1:18]
dat[, 17:ncol(dat)] -> mets
mets[1:4,1:4]
mets <- as.matrix(mets)
mets <- mets[ ,colSums(is.na(mets)) == 0] # rid matrix of mzs with missing data
# remove metabolites identified by Danijel as artifacts
ghost_peaks <- c('4-Trimethylammoniobutanoate', 'X4.Trimethylammoniobutanoate', 'Sorbitol', 'Xanthosine', 'Oxypurinol', 'Uridine', 'Trimethylamine', 'Valine', 'Adenine', 'Pipecolic acid', 'Pipecolic.acid') # these were one member of a few pair of highly-correlated metabolties (r>0.98 in another analysis),Danijel suspected that these metabolties are actualy LC-MS artifacts
mets <- mets[ ,!colnames(mets) %in% ghost_peaks] 
centered_mets <- t(scale(t(mets), scale=F)) # center, but don't scale metabolites by sample
rm(ghost_peaks)

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
colnames(dat)[1:18]
dat <- cbind(dat[ ,1:16], centered_mets)
rm(mets)
dat[1:5, 1:20]
rm(centered_mets)

######################################################################################################
# assess the covariance of metabolites, after correcting for phylogenetic correlation:

## PIC including intraspecies covariance
colnames(dat)[1:18]
traits <- names(dat)[c(9:12, 17:ncol(dat))] 

groups <- paste(dat$sex, dat$age.days)
group <- names(table(groups))

for(k in group) {
tmp.dat <- dat[groups == k, ]
no.species <- length(table(tmp.dat$species))

tmpics <- matrix(ncol=no.species-1, nrow=length(traits))

for(i in 1:length(traits)) {
  trait <- traits[i]
  tmp <- tmp.dat[ ,trait]
  names(tmp) <- tmp.dat$species
  tmp <- tmp[!is.na(tmp)]
  tmp <- aggregate(tmp, by=list(names(tmp)), list)
  names(tmp$x) <- tmp$Group.1
  tmp  <- tmp$x
  tmp.tree <- keep.tip(flytree, names(tmp))
  tmp  <- tmp[tmp.tree$tip.label]
  tmpics[i, ] <- pic.ortho(tmp, tmp.tree, intra = T, var.contrasts = F) }
   # method from Felsenstein (2008) Am. Nat.
rownames(tmpics) <- traits
write.table(tmpics, file=paste('correlation among the metabolome/pics', k))
}

rm(tmp)
rm(tmp.tree)
rm(tmp.dat)

picM5 <- read.table('correlation among the metabolome/pics M 5')
picM31 <- read.table('correlation among the metabolome/pics M 31')
picF5 <- read.table('correlation among the metabolome/pics F 5')
picF31 <- read.table('correlation among the metabolome/pics F 31')

mzs <- rownames(picF5)[5:101]

##############################################################################
nSets = 4
PICs <- list(picF5, picF31, picM5, picM31)
setLabels = c("young female", "old female", 'young male', 'old male')
shortLabels = c('F 5', 'F 31' ,'M 5', 'M 31')

rm(list=c('picF31', 'picF5', 'picM31', 'picM5'))

# plot correlation matrices of metabolites for each sex x age:
F5_pal <- colorRampPalette(c("grey", "white", "blue"))(n = 50)
F31_pal <- colorRampPalette(c("grey", "white", "purple"))(n = 50)
M5_pal <- colorRampPalette(c("grey", "white", "brown"))(n = 50)
M31_pal <- colorRampPalette(c("grey", "white", "darkgoldenrod"))(n = 50)
setPals <- list(F5_pal, F31_pal, M5_pal, M31_pal)

for(set in 1:nSets) {
  tmpic <- apply(PICs[[set]], 2, unlist)
  pic.cor <- cor(t(tmpic[ ,5:ncol(tmpic)]), method='pearson') # excludes traits 1:4, which are weight and lifespan data
  heatmap.2(pic.cor, trace='n', col=unlist(setPals[set]), key.title = 'cor PIC-PIC', margins = c(10,12), key.xlab = 'pearson (r)', cexRow=0.5, cexCol=0.5, keysize=1, dendrogram='none', main=paste('PIC correlation', setLabels[set])) }
##############################################################################

# PERMUTATION test: compare the distribution of correlation among feature PICs:

setLabels
heatmapcol4 <- c("#CD6600", "#a10202", "#2a44ad", "#AE02E6") # darker colors used later for heatmaps:

for(set in 1:nSets) {
  tmpic <- apply(PICs[[set]], 2, unlist)
  pic.cor <- cor(t(tmpic[ ,5:ncol(tmpic)]), method='pearson')

# plot density of correlations among real data in comparison to correlations among 100 permutations

  d <- density(pic.cor[lower.tri(pic.cor, diag=F)])
if(set==1) { plot(d, ylim=c(0, 1), las=1, xlab='metabolite PIC, pariwise correlation', col=heatmapcol4[set], lwd=2, ylab='density', main='')}
if(set>1) { lines(d, col=heatmapcol4[set], lwd=2) } 
}
legend('topright', legend = c(setLabels, 'permutation'), pch=16, col=c(heatmapcol4, 8), bty='n', cex=1)

# test distribution of correlations between PICs by permuting.  In this case, try conservative permutation in which only the species labels are randomizes, after species means are calculated and, most conservatively, while conserving the metabolite-metabolite relatiopnships within each species

k <- shortLabels[4] # I chose to permute the 'old male' metabolome becuase it showed the largest excess of correalted PICs
tmp.dat <- dat[groups == k, ]
no.species <- length(table(tmp.dat$species))

tmpics <- matrix(ncol=no.species-1, nrow=length(traits))

for(i in 1:20) {
  randomised.species <- sample(unique(tmp.dat$species)) # conservative permutation, where species names are shuffled prior to PIC calculation, but after species means are estimated.  THis is also keeping all metabolites within a speices together [ this does not disrupt the mz-mz relationship, just the effect of phylogeny]
  
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
  lines(density(corperm[lower.tri(corperm, diag=F)]), col=rgb(0,0,0,0.1), lwd=2)  }

# re-plot the data so that its 'on top' of the permuted data on the plot:
for(set in 1:nSets) {
  tmpic <- apply(PICs[[set]], 2, unlist)
  pic.cor <- cor(t(tmpic[ ,5:ncol(tmpic)]), method='pearson')
  d <- density(pic.cor[lower.tri(pic.cor, diag=F)])
  lines(d, col=reorderedmycol4[set], lwd=2) }


var.test(pic.cor[lower.tri(pic.cor, diag=F)], corperm[lower.tri(corperm, diag=F)]) # an F-test to compare the variances in both populations

rm(corperm)
rm(pic.cor)
rm(randomised.species)
rm(mytranscol4)

#############################################################
# try ls ~ mz univariate major axis regression models 
#############################################################
ls_mz <- list()

trait.of.interest <- 'mean.ls'

ls.by.mz <- function(x) {
  m <- sma(ls ~ x -1)
  return(c(m$r2, m$pval)) }

for(set in 1:nSets){
  tmpic <- apply(PICs[[set]], 2, unlist)  
  ls <- t(tmpic)[ ,trait.of.interest]
  mz.pics <- t(tmpic)[ ,5:101]
  
  out <- matrix(nrow=97, ncol=2)
  for(i in 1:97) {
    mz <- mz.pics[ ,i]
    out[i, ] <- unlist(ls.by.mz(mz)) }
  
  out <- as.data.frame(out)
  rownames(out) <- mzs
  colnames(out) <- c('r.sq','pval')
  out$FDR <- p.adjust(out$pval, method='BH', n=194) # NOTE there are 194 comparisons to each set of mean lifespan PICs, this is becuase each sex (for which there is only one mean lifespan) is compared to 97 metabolties from 2 ages.  2*97 = 194
  ls_mz[[set]] <- out

write.table(out, file=paste0(trait.of.interest, '_lifespanPIC_by_mzPIC_major.axis.regression_', setLabels[set])) }


fdrs <- sapply(ls_mz, getElement, name="FDR")
table(fdrs < 0.1)
apply(fdrs, 2, function(x) {x[order(x)][1:5] }) 

pees <- sapply(ls_mz, getElement, name="pval")
apply(pees, 2, function(x) {x[order(x)][1:5] }) # the lowest 5 P values per set
apply(pees, 2, function(x) {mzs[order(x)[1:3]] }) # the metabolites with the lowest 3 P values per set


# plot the 'top' 3 'hits'; the mzs with the largest (or smallest) correlation coefficients [not the same as the lowest P]  for each set:
par(mfrow=c(4,3))

top.mzs <- apply(pees, 2, function(x) {mzs[order(x)[1:3]] }) # the metabolites with the lowest 3 pvalues per set

for(set in 1:nSets){
  tmpic <- apply(PICs[[set]], 2, unlist)

for(i in 1:3){
  mz.of.interest <- top.mzs[i, set]
  
  plot(tmpic['mean.ls', ] ~  tmpic[mz.of.interest, ], pch=19, ylab=ifelse(i==1, paste(trait.of.interest, '(PIC)'), ''), xlab=paste('metabolite (PIC)'), las=1, main=mz.of.interest, col=reorderedmycol4[set])
  major.axis.regression <- sma(tmpic['mean.ls', ] ~  tmpic[mz.of.interest, ] -1)
  pvalue <- paste0('P=', round(unlist(major.axis.regression$pval), 4))
  fdr <- paste0('FDR=', round(p.adjust(c(major.axis.regression$pval, runif(96, min=0.1, max=1)), method='BH')[1], 3))
  legend(ifelse(unlist(major.axis.regression$coef)[2] >0, 'bottomright', 'topright'), legend=c(pvalue, fdr), bty='n', cex=0.8)
  abline(sma(tmpic['mean.ls', ] ~  tmpic[mz.of.interest, ]-1), col=8)}
}


ls_by_mz.table <- cbind(ls_mz[[1]], ls_mz[[2]], ls_mz[[3]], ls_mz[[4]])
colnames(ls_by_mz.table) <- paste(colnames(ls_by_mz.table), setLabels[c(1,1,1,2,2,2,3,3,3,4,4,4)])

ls_by_mz.table$mz <- rownames(ls_by_mz.table)
ls_by_mz.table <- ls_by_mz.table[ ,c(13,1:12)]

write.table(ls_by_mz.table, file=paste0(trait.of.interest, 'lifespanPIC_by_mzPIC_major.axis.regression_Table.csv'), sep=',', row.names=F) 
rsq.mat <- as.matrix(cbind(ls_mz[[1]][1], ls_mz[[2]][1], ls_mz[[3]][1], ls_mz[[4]][1]))
fdr.mat <- as.matrix(cbind(ls_mz[[1]][3], ls_mz[[2]][3], ls_mz[[3]][3], ls_mz[[4]][3]))
colnames(rsq.mat) <- setLabels
colnames(fdr.mat) <- setLabels

colfunc <- colorRampPalette(c('white', 'blue'))

heatmap.2(rsq.mat, margin=c(12,10), trace='n', col=colfunc(20), Colv=F, main='type II regression\nPIC ls ~ PIC mz', key.xlab = 'r-squared')

colSums(fdr.mat < 0.5)
colSums(fdr.mat < 0.25)
rowSums(fdr.mat < 0.5)[rev(order(rowSums(fdr.mat < 0.5)))][1:12]

heatmap.2(rsq.mat[rowSums(fdr.mat < 0.5) > 0, ], margin=c(12,14), trace='n', col=colfunc(20), Colv=F, key.xlab = 'r-squared', density.info = 'n', keysize = 1, key.title = '', cexRow=1, cexCol=1)

par(mfrow=c(2,1))
fdrs <- unlist(c(ls_mz[[1]][3], ls_mz[[2]][3], ls_mz[[3]][3], ls_mz[[4]][3]))
fdrs[order(fdrs)]
par(mar=c(5,4,4,2))
plot(fdrs[order(fdrs)], pch=19, las=1, main='FDR for 97 over 4 groups metabolites (n=194 comparisons)')

tmp <- rbind(colSums(fdr.mat < 0.5), 97-colSums(fdr.mat < 0.5))
rownames(tmp) <- c('FDR<0.5', 'FDR>0.5')
tmp
barplot(tmp, col=c(2, 'lightgrey'), border="white", space=0.04, font.axis=2, las=1, ylab='metabolites', main='metabolites FDR<0.5')

rm(fdr)
rm(fdrs)
rm(ls_mz)
rm(ls_by_mz.table)
rm(rsq.mat)
rm(fdr.mat)
rm(tmpic)
rm(tmpics)
rm(top.mzs)
rm(out)
#############################################################


#############################################################
# multivariate association.  What about ls ~ metabolite modules?
#############################################################

## use WGCNA to pick modules
nSets=4
picM5 <- read.table('correlation among the metabolome/pics M 5')
picM31 <- read.table('correlation among the metabolome/pics M 31')
picF5 <- read.table('correlation among the metabolome/pics F 5')
picF31 <- read.table('correlation among the metabolome/pics F 31')

mzs <- rownames(picF5)[5:101]

# add KEGG and other database names to metabolites:
mz.info <- read.table('correlation among the metabolome/mz.info.txt', header=T, sep=',')
head(mz.info)

rownames(picM5)[!rownames(picM5) %in% mz.info$metabolite]

rownames(picM5) <- gsub('.ac', ' ac', rownames(picM5), fixed=T)
rownames(picM5) <- gsub('.Acid', ' acid', rownames(picM5), fixed=T)
rownames(picM5) <- gsub('.', '-', rownames(picM5), fixed=T)
rownames(picM5) <- gsub('D-Glyceraldehyde-3-phosphate', 'D-Glyceraldehyde 3-phosphate', rownames(picM5), fixed=T)
rownames(picM5) <- gsub('gamma-Aminobutyric acid', 'gamma-Aminobutyric Acid', rownames(picM5), fixed=T)
rownames(picM5) <- gsub('Glucose-6-phosphate', 'Glucose 6-phosphate', rownames(picM5), fixed=T)
rownames(picM5) <- gsub('Glycerol-3-phosphate', 'Glycerol 3-phosphate', rownames(picM5), fixed=T)
rownames(picM5) <- gsub('Oxidized-glutathione', 'Oxidized glutathione', rownames(picM5), fixed=T)

rownames(picM5)[!rownames(picM5) %in% mz.info$metabolite]
mz.info$metabolite[!mz.info$metabolite %in% rownames(picM5)]

mzs <- rownames(picM5)

rownames(picF5) <- mzs 
rownames(picF31) <- mzs 
rownames(picM31) <- mzs 


# Form a multi-set structure that will hold the metabolite data and the organismal phenotypes (i.e. lifespan and weight).
multiPIC = vector(mode = "list", length = nSets) # empty entity

multiPIC[[1]] = list(data = as.data.frame(t(picF5[-c(1:4), ]))) # remove the trait data (first 4 rows) and transpose
names(multiPIC[[1]]$data) = rownames(picF5)[5:nrow(picF5)]
multiPIC[[2]] = list(data = as.data.frame(t(picF31[-c(1:4), ]))) # remove the trait data (first 4 rows) and transpose
names(multiPIC[[2]]$data) = rownames(picF31)[5:nrow(picF31)]
multiPIC[[3]] = list(data = as.data.frame(t(picM5[-c(1:4), ]))) # remove the trait data (first 4 rows) and transpose
names(multiPIC[[3]]$data) = rownames(picM5)[5:nrow(picM5)]
multiPIC[[4]] = list(data = as.data.frame(t(picM31[-c(1:4), ]))) # remove the trait data (first 4 rows) and transpose
names(multiPIC[[4]]$data) = rownames(picM31)[5:nrow(picM31)]

Traits = vector(mode="list", length = nSets)

Traits[[1]] <- t(picF5[1:4, ])
Traits[[2]] <- t(picF31[1:4, ])
Traits[[3]] <- t(picM5[1:4, ])
Traits[[4]] <- t(picM31[1:4, ])

gc()

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

save(multiPIC, Traits, nGenes, nSamples, setLabels, shortLabels, exprSize, file = "correlation among the metabolome/metabolomeConsensus-dataInput.RData") # saves the data in a particular format

rm(multiPIC)
rm(Traits)
rm(nGenes)
rm(nSamples)
rm(setLabels)
rm(shortLabels)
rm(exprSize)


load("correlation among the metabolome/metabolomeConsensus-dataInput.RData")

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
softPower = 7 # chosen power 
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

rm(adjacencies)


# a comparison of module membership differences between groups 
##############################################################################

### analyze TOM network properties
network_colors <-  c("#EEC900", "#CD6600", "#87CEFA", "#104E8B")

par(mfrow=c(4,5))
cutoffs <- c(0.05, 0.1, 0.15, 0.2, 0.4)
par(mar=c(4,4,4,4))

for(set in 1:nSets) {
  for(cut in cutoffs) {
    adj <-  TOM[set, , ] 
    dimnames(adj) = list(colnames(multiPIC[[set]]$data), colnames(multiPIC[[set]]$data))
    diag(adj) <- NA
    adj[adj > cut] = 1  
    adj[adj != 1] = 0
    
    network <- graph.adjacency(adj)
    network <- igraph::simplify(network)  # removes self-loops
    
    plot(network, edge.arrow.size = 0, vertex.size=6, vertex.label=NA, vertex.color=network_colors[set], main=paste0(setLabels[set], ', cor^7 > ', cut))
    
  }}

for(set in 1:nSets) {
  for(cut in cutoffs) {
    adj <-  TOM[set, , ]
    dimnames(adj) = list(colnames(multiPIC[[set]]$data), colnames(multiPIC[[set]]$data))
    diag(adj) <- NA
    adj[adj > cut] = 1  
    adj[adj != 1] = 0
    
    network <- graph.adjacency(adj)
    network <- igraph::simplify(network)  # removes self-loops
    
    par(mar=c(4,6,4,5))
    plot(degree.distribution(network, mode='in'), pch=19, log='xy', las=1, main=paste0(setLabels[set], ', cor^7 > ', cut), ylab='degree distribution') # should show a negative relationship that is linear on this log-log plot
  }}



par(mfrow=c(4,2))

cutoff = 0.15 # cutoff chosen based on the analysis above
# WGCNA will use its own criteria for module identification, but this cutoff can be used later for network visualization (ie. in cytoscape or igraph)

for(set in 1:nSets) {
  adj <-  TOM[set, , ]
  dimnames(adj) = list(colnames(multiPIC[[set]]$data), colnames(multiPIC[[set]]$data))
  diag(adj) <- NA
  adj[adj > cutoff] = 1  
  adj[adj != 1] = 0
  
  network <- graph.adjacency(adj)
  network <- igraph::simplify(network)  # removes self-loops
  
  par(mar=c(3,4,4,2))
  plot(network, edge.arrow.size = 0, vertex.size=6, vertex.label=NA, vertex.color=network_colors[set], main=paste0(setLabels[set], ', cor^7 > ', cutoff))
  par(mar=c(6,6,6,6))
  plot(degree.distribution(network, mode='in'), pch=19, log='xy', las=1, ylab='fraction of nodes', xlab='degree') # should show a negative relationship that is linear on this log-log plot
}


#############################
# Identify Modules:
#############################
# loop though the sets to identify modules:
dynamicMods <- list()
par(mfrow=c(2,2))
par(mar=c(6,6,6,6))

mzs <- colnames(multiPIC[[4]]$data)

minModuleSize = 8 # set module size; I used the output of the code below to help choose the minimum mod size.  I used the heatmap with the modules shown on the margins to reveal how the modules related to the covariance among mzs.  I favored a number that was small enough to give clusters that were mosly non 'broken up', but a number that was not so small that what looked like separate clusters in the heatmap came out as modules in WGCNA

for(set in 1:nSets) {
  mzTree <- hclust(as.dist(1-TOM[set, , ]), method = "average")
  plot(mzTree, xlab="", sub="", main = setLabels[set], labels = F, las=1, hang = 0.04) 

    # Module identification and assignment using dynamic tree cut:
  dynamicMods[[set]] = cutreeDynamic(dendro = mzTree, distM = 1-TOM[set, , ], deepSplit = 3, pamRespectsDendro = T, minClusterSize = minModuleSize) }

lapply(dynamicMods, table) # # module assignments for metabolites, class 0 are metabolites that weren't assigned to a module

dynamicColors = lapply(dynamicMods, labels2colors) # Convert numeric labels into colors
lapply(dynamicColors, table)

plotDendroAndColors(mzTree, dynamicColors[[set]], "Dynamic Modules", hang = 0.03, las=1, main=setLabels[set], dendroLabels = mzs)

#WGCNA assigns module colors in order of module size 'turquoise' being the largest.  [NOTE]: 'grey' is reserved for all features that don't make it into a module, and grey may actually be larger that turquoise.

# replace module colors with letters, A,B,C.... depending order of module size (but keep 'grey' as 'grey')
module_renaming_vector <- c('grey', LETTERS[1:6])
names(module_renaming_vector) <- 0:6
module_renaming_vector

renamedMods <- dynamicMods
for(set in 1:nSets) {
  renamedMods[[set]] <- module_renaming_vector[dynamicMods[[set]]+1] }

renamedMods


mzTree.order <- list()

# use trees and heatmaps to evaluate the modules:
for(set in 1:nSets) {
  tmptree <- hclust(as.dist(1-TOM[set, , ]), method='average')
  mzTree.order[[set]] <- tmptree$order
  plotDendroAndColors(tmptree, dynamicColors[[set]], "Modules", hang = 0.03, las=1, main=setLabels[set], dendroLabels = mzs) 
}


dev.off() 

# heatmap for figure (large enough to see the metabolite names along with partitions to distinguish modules)
for(set in 1:nSets) {
  x <- abs(cor(multiPIC[[set]]$data[ ,mzTree.order[[set]]]))
  x[1:4,1:4]
  diag(x) = NA
  
  # cumbersome code to make partitions between modules on the ol' heatmap:
  z <- as.numeric(as.factor(dynamicColors[[set]][mzTree.order[[set]]]))
  zz <- z-c(0, z)!=0 
  zzz <- as.logical(c(zz[2:length(z)], 'FALSE') )
  
  png(file=paste(setLabels[set], "heatmap plus modules.png"), width = 2000, height = 1700, units = "px", res=200, bg = "transparent")
  
  pheatmap(x, color = colorRampPalette(brewer.pal(n = 7, name ="PuRd"))(100), fontsize=6, cluster_cols = F, cluster_rows = F, gaps_row = which(zzz), gaps_col = which(zzz), border_color=NA, show_colnames = F, main=setLabels[set])
  
  dev.off()}

renamedMods[[set]][mzTree.order[[set]]]# list of the module membership for the metabolites in 'set'


## replot smaller and without metabolite labels to use as main figure (which will be too small to show metabolite names)
for(set in 1:nSets) {
  colfunc <- colorRampPalette(c('white', heatmapcol4[set]))  
  
  x <- abs(cor(multiPIC[[set]]$data[ ,mzTree.order[[set]]]))
  x[1:4,1:4]
  diag(x) = NA
  x[lower.tri(x)] <- x[lower.tri(x)]^7 # this line makes the lower triangle show the correlations raise to the 'softPower' used in the module identification
  
  legend_names <- unique(unlist(renamedMods[set]))
  legend_names[legend_names=='grey'] <- NA
  legend_names <- sort(legend_names, na.last=T)
  legend_names[is.na(legend_names)] <- 'non-modular'
  
  legend_colors <- unique(unlist(dynamicColors[[set]]))
  legend_colors[legend_colors=='grey'] <- NA
  legend_colors <- sort(legend_colors, na.last=T)
  legend_colors[is.na(legend_colors)] <- 'grey'
  
  jpeg(file=paste(setLabels[set], "heatmap plus modules_small for main figure.jpg"), width = 1200, height = 1200, units = "px", res=100,
       quality = 1000)

heatmap.2(x, trace='none', dendrogram = 'none', main=setLabels[set], RowSideColors= dynamicColors[[set]][mzTree.order[[set]]], Rowv=F, Colv=F, col=colfunc, key.xlab='abs(cor)', margin=c(5,5), keysize=1, key.title = '', density.info = 'n', labRow='', labCol='', xlab=expression(paste("r"^"7")), ylab='r') 

dev.off()}



save(TOM, dynamicColors, renamedMods, file = 'correlation among the metabolome/WGCNA_module_identification.RData')
rm(TOM)
rm(dynamicColors)
rm(renamedMods)


#############################################################################################
#############################################################################################
load("correlation among the metabolome/metabolomeConsensus-dataInput.RData") # load data
load('correlation among the metabolome/WGCNA_module_identification.RData') # load modules
lapply(dynamicColors, table)

picM5 <- read.table('correlation among the metabolome/pics M 5')
picM31 <- read.table('correlation among the metabolome/pics M 31')
picF5 <- read.table('correlation among the metabolome/pics F 5')
picF31 <- read.table('correlation among the metabolome/pics F 31')
mzs <- rownames(picF5)[5:101]

PICs <- list(picF5, picF31, picM5, picM31)
rm(list=c('picF31', 'picF5', 'picM31', 'picM5'))

nSets = 4
setLabels = c("young female", "old female", 'young male', 'old male')
shortLabels = c('F 5', 'F 31' ,'M 5', 'M 31')


########################################################
# analyze the similarity/disimiarity between networks:
########################################################
# several approaches
# approach 1: use the WGCNA 'overlapTable()' to get the overlap between modules in each group
combinations <- print(combn(1:4,2)) # matrix of unique pairwise combinations

global_colors = circlize::colorRamp2(c(0, 1), c("white", "red")) # set a global color scale for all heatmaps that covers 0 to 1

# plot heatmaps of the number of intersecting metabolites:
htlist=NULL

for(i in 1:3){
tmp <- overlapTable(renamedMods[[combinations[1, i]]], renamedMods[[combinations[2, i]]], ignore = 'grey')
ct <- tmp$countTable
column_ha = HeatmapAnnotation('n' = anno_barplot(colSums(ct), border=F))
row_ha = rowAnnotation('n' = anno_barplot(rowSums(ct), border=F))

# to calculate the proportion of intersecting metabolites, you'll need to know for each pair, which module has the LEAST metabolites.  This number defines the MAXIMUM intersection:
r <- matrix(rowSums(ct), byrow=F, nr=nrow(ct), nc=ncol(ct))
y <- matrix(colSums(ct), byrow=T, nr=nrow(ct), nc=ncol(ct))
y[which(r <= y, arr.ind = T)] <- r[which(r <= y, arr.ind = T)]

htlist = htlist + Heatmap(ct/y, name = "intersection", top_annotation = column_ha, right_annotation = row_ha, cluster_rows = F, cluster_columns = F, col=global_colors, row_title=setLabels[combinations[1, i]], column_title = setLabels[combinations[2, i]], row_names_side = 'left', column_title_side = 'bot', width = unit(4, "cm"), height = unit(4, "cm")) 
}

draw(htlist, ht_gap = unit(2, "cm"))


htlist2=NULL

for(i in 4:5){
  tmp <- overlapTable(renamedMods[[combinations[1, i]]], renamedMods[[combinations[2, i]]], ignore = 'grey')
  ct <- tmp$countTable
  column_ha = HeatmapAnnotation('n' = anno_barplot(colSums(ct), border=F))
  row_ha = rowAnnotation('n' = anno_barplot(rowSums(ct), border=F))
  
  # to calculate the proportion of intersecting metabolites, you'll need to know for each pair, which module has the LEAST metabolites.  This number defines the MAXIMUM intersection:
  x <- matrix(rowSums(ct), byrow=F, nr=nrow(ct), nc=ncol(ct))
  y <- matrix(colSums(ct), byrow=T, nr=nrow(ct), nc=ncol(ct))
  y[which(x <= y, arr.ind = T)] <- x[which(x <= y, arr.ind = T)]
  y
  
  htlist2 = htlist2 + Heatmap(ct/y, name = "intersection", top_annotation = column_ha, right_annotation = row_ha, cluster_rows = F, cluster_columns = F,  col=global_colors, row_title=setLabels[combinations[1, i]], column_title = setLabels[combinations[2, i]], row_names_side = 'left', column_title_side = 'bot', width = unit(4, "cm"), height = unit(4, "cm")) 
}

draw(htlist2, ht_gap = unit(2, "cm"))


htlist3=NULL

for(i in 6){
  tmp <- overlapTable(renamedMods[[combinations[1, i]]], renamedMods[[combinations[2, i]]], ignore = 'grey')
  ct <- tmp$countTable
  column_ha = HeatmapAnnotation('n' = anno_barplot(colSums(ct), border=F))
  row_ha = rowAnnotation('n' = anno_barplot(rowSums(ct), border=F))
   
  # to calculate the proportion of intersecting metabolites, you'll need to know for each pair, which module has the LEAST metabolites.  This number defines the MAXIMUM intersection:
  x <- matrix(rowSums(ct), byrow=F, nr=nrow(ct), nc=ncol(ct))
  y <- matrix(colSums(ct), byrow=T, nr=nrow(ct), nc=ncol(ct))
  y[which(x <= y, arr.ind = T)] <- x[which(x <= y, arr.ind = T)]
  y
  
  htlist3 = htlist3 + Heatmap(ct/y, name = "intersection", top_annotation = column_ha, right_annotation = row_ha, cluster_rows = F, cluster_columns = F,  col=global_colors, row_title=setLabels[combinations[1, i]], column_title = setLabels[combinations[2, i]], row_names_side = 'left', column_title_side = 'bot', width = unit(4, "cm"), height = unit(4, "cm")) 
}

draw(htlist3, ht_gap = unit(2, "cm"))



tmp <- overlapTable(renamedMods[[combinations[1, i]]], renamedMods[[combinations[2, i]]], ignore = 'grey')
ct <- tmp$countTable
pt <- tmp$pTable <= 0.05

ct

# to calculate the proportion of intersecting metabolites, you'll need to know for each pair, which module has the LEAST metabolites.  This number defines the MAXIMUM intersection:
x <- matrix(rowSums(ct), byrow=F, nr=nrow(ct), nc=ncol(ct))
y <- matrix(colSums(ct), byrow=T, nr=nrow(ct), nc=ncol(ct))
y[which(x <= y, arr.ind = T)] <- x[which(x <= y, arr.ind = T)]
y

Heatmap(ct/y, name = "intersection", top_annotation = column_ha, right_annotation = row_ha, cluster_rows = F, cluster_columns = F, col=colorRampPalette(c("white", "red"))(10), row_title=setLabels[combinations[1, i]], column_title = setLabels[combinations[2, i]], row_names_side = 'left', column_title_side = 'bot', width = unit(4, "cm"), height = unit(4, "cm"), 
  cell_fun = function(j, i, x, y, w, h, fill) {
  if(pt[i, j] == T) {
    grid.text("*", x, y)
  }
})


#############################################################################################
# Calculate eigengenes (PC1 of each module)
#############################################################################################
eigenMZs <- list()
modColors <- list()
modLabels <- list()
trees <- list()

for(set in 1:nSets) {
  MEList = moduleEigengenes(as.data.frame(multiPIC[[set]]), colors = renamedMods[[set]])
  MEs = MEList$eigengenes
  MEDiss = 1-cor(MEs) # Calculate dissimilarity of module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average") # Cluster module eigengenes
  
  # clustering of eigenmetabolites (first PC of each module)
  plot(METree, main = "Clustering of module eigenmetabolites", xlab = "", sub = "")
  
  # optional merging step: modules that are similar enough here can be merged
  MEDissThres = 0.25 # threshold for merging
  # Plot the cut line into the dendrogram
  abline(h=MEDissThres, col = "red")
  # Call an automatic merging function
  merge = mergeCloseModules(as.data.frame(multiPIC[[set]]), dynamicColors[[set]], cutHeight = MEDissThres, verbose = 3)
  merge.mods = mergeCloseModules(as.data.frame(multiPIC[[set]]), renamedMods[[set]], cutHeight = MEDissThres, verbose = 3)
  
  mergedColors = merge$colors # The merged module colors
  merge.mods$colors # The merged module names
  merged.modMEs = merge.mods$newMEs # Eigengenes of the new merged modules:
  
  mzTree <- hclust(as.dist(1-TOM[set, , ]), method = "average")
  plotDendroAndColors(mzTree, cbind(dynamicColors[[set]], mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), hang = 0.03, las=1)
  
  
  merge.mods$colors
  moduleColors =  merge.mods$colors # Rename to moduleColors
  colorOrder = c('grey', LETTERS[1:20]) # Construct numerical labels corresponding to the colors
  moduleLabels = match(moduleColors, colorOrder)-1
  MEs =  merged.modMEs
  
  # save network entities
  eigenMZs[[set]] <- MEs
  modColors[[set]] <- moduleColors
  modLabels[[set]] <- moduleLabels
  trees[[set]] <- mzTree
  
  # and on the disk:
  save(MEs, moduleLabels, moduleColors, mzTree, file = paste0('correlation among the metabolome/', shortLabels[set], '_networkConstruction.RData')) }
#############################################################################################

########################################################################
# Relate modules to traits by major axis regression and identifying important metabolites:
########################################################################

# load eigen modules:
eigenMZs <- list()
for(set in 1:nSets) {
  shortLabels[set]
  load(paste0('correlation among the metabolome/', shortLabels[set], '_networkConstruction.RData'))
  eigenMZs[[set]] <- MEs
}
eigenMZs


ls_mod <- list()

ls.by.mz <- function(x) {
  m <- sma(ls ~ x -1)
  return(c(m$r2, m$pval)) }

nModulesFemale <- ncol(eigenMZs[[1]]) + ncol(eigenMZs[[2]])
nModulesMale <- ncol(eigenMZs[[3]]) + ncol(eigenMZs[[4]])

for(set in 1:nSets){
  tmpic <- apply(PICs[[set]], 2, unlist)  
  ls <- t(tmpic)[ ,'mean.ls']
  
  out <- matrix(nrow=ncol(eigenMZs[[set]]), ncol=2)
  for(i in 1:ncol(eigenMZs[[set]])) {
    eigenmod <- unlist(eigenMZs[[set]][i])
    out[i, ] <- unlist(ls.by.mz(eigenmod)) }
  
  out <- as.data.frame(out)
  rownames(out) <- colnames(eigenMZs[[set]])
  colnames(out) <- c('r.sq','pval')
  out$FDR <- p.adjust(out$pval, method='BH', n=ifelse(set %in% c(1:2), nModulesFemale, nModulesMale)-2) # I subtract 2 because we will not test the grey module from each sex/age
  out$group <- setLabels[set]
  out <- out[ ,c(4,1:3)]
  ls_mod[[set]] <- out
  
  write.table(out, file=paste0('correlation among the metabolome/lifespanPIC_by_eigenmodules_major.axis.regression', setLabels[set])) }


ls_eigenmod_table <- rbind(ls_mod[[1]], ls_mod[[2]], ls_mod[[3]], ls_mod[[4]])
ls_eigenmod_table$module <- rownames(ls_eigenmod_table)
ls_eigenmod_table <- ls_eigenmod_table [ ,c(1, 5, 2:4)]

write.table(ls_eigenmod_table, file='correlation among the metabolome/lifespanPIC_by_eigenmodules_major.axis.regression_Table.csv', sep=',', row.names = F)


##########################################################################
# plot to examine the relationship between trait (i.e. mean-ls) and a module eigenmetabolite

# the following choices are manually coded based on a module and trait of interest:
par(xpd=FALSE)

set_of_interest <- 'old male'
trait_of_interest <- 'mean-ls'
module_of_interest <- 'D'
##########################################################################
set <- seq(1:4) [setLabels == set_of_interest]
y <- Traits[[set]][ ,trait_of_interest]
x <- eigenMZs[[set]][ ,paste0('ME', module_of_interest)]

par(mfrow=c(1,1))
plot(y ~ x, pch=16, las=1, main=set_of_interest,  xlab=paste('eigenmodule', module_of_interest), ylab=trait_of_interest)
abline(sma(y ~ x -1))
summary(s <- sma(y ~ x -1))
legend('bottomright', bty='n', legend=c('type II regression', paste('P=', round(unlist(s$pval), 4)), paste('r^2=', round(unlist(s$r2), 3))))

########################################################################
# network visualization using cytoscape
########################################################################
for(set in 1:nSets) {
  
  cyt.all = exportNetworkToCytoscape(TOM[set, , ], 
                                     edgeFile = paste('Cytoscape', setLabels[set], "edges.txt"), 
                                     nodeFile = paste('Cytoscape', setLabels[set], "nodes.txt"),
                                     weighted = TRUE,
                                     threshold = 0.15, # determines extent of connectivity, this can be chosen based on the plots of degree distribution (above)
                                     nodeNames = mzs[5:101],
                                     altNodeNames = mzs[5:101],
                                     nodeAttr = modColors[[set]])
}

cyt.all

# three ideas for labeling the modules, 1) keep them as-is (i.e. modules A, B ...D, etc), 2) label the modules based on the most connected (hub) metabolite (inspired by Guimera and Amaral Nature 2005), or 3) use the top pathway from enricment analysis.

# identify the hubs:

cutoff = 0.15 # cutoff chosen based on degree distribution analysis (way above this line)

metabolite<- mzs[5:101]

for(set in 1:nSets) {
adj <-  TOM[set, , ]
dimnames(adj) = list(colnames(multiPIC[[set]]$data), colnames(multiPIC[[set]]$data))
diag(adj) <- NA
adj[adj > cutoff] = 1  
adj[adj != 1] = 0
network <- graph.adjacency(adj)
network <- igraph::simplify(network)  # removes self-loops

tmp <- data.frame(metabolite, module=as.factor(modColors[[set]]))
deg <- degree(network, v = V(network), mode = "all", loops = F, normalized = F)
tmp$degree <- deg[tmp$metabolite]
tmp <- merge(aggregate(degree ~ module, tmp, max), tmp, by=c("module","degree"))
tmp$group <- setLabels[set] 
ifelse(set==1, hubs <- tmp, hubs <- rbind(hubs, tmp))}

table(hubs$module, hubs$group)
hubs[hubs$group == 'old female', ]

  
cyt.all = exportNetworkToCytoscape(TOM[set, , ], 
                                     edgeFile = paste('Cytoscape', setLabels[set], "withhubs_edges.txt"), 
                                     nodeFile = paste('Cytoscape', setLabels[set], "withhubs_nodes.txt"),
                                     weighted = TRUE,
                                     threshold = 0.15, 
                                     nodeNames = mzs[5:101],
                                     altNodeNames = mzs[5:101],
                                     nodeAttr = list(modColors[[set]], tmp$hub))

cyt.all


######################################################################################################################
# use FELLA, tutorial at: http://127.0.0.1:11793/library/FELLA/doc/musmusculus.pdf (modify code for Dme)
######################################################################################################################

set.seed(1)

graph <- buildGraphFromKEGGREST(organism = "dme", filter.path = c("01100", "01200", "01210", "01212", "01230"))  # Build network and filter out the large overview pathways:

tmpdir <- "~/Google Drive File Stream/My Drive/Documents/Hoffman_Metabolite_Phylogeny_Paper/FELLA_database"

unlink(tmpdir, recursive = TRUE) # this line allows you to rewrite the: buildataFromGraph()

buildDataFromGraph(graph, databaseDir = tmpdir, internalDir = FALSE, matrices = 'diffusion', normality = "diffusion", niter = 1000) # using matrices=listMethods() here leads to three matrices being made: 'diffusion', ' hypergemetric' and 'pagerank'.  It takes longer than building with a single method.  Note that 'normality' can either be 'diffusion' or 'pagerank', and this will determine which of those two analyses can run based on these data.  I ran normality as 'diffusion' and this caused an error saying that Z-scores would not be available for the 'pagerank' method.  It did however make a 'pagerank' matrix, so I suspect that I can still do non-parametric analysis by pagerank, rather than a parametric analysis that might depend on the outcome of this line. [NOTE]: I ultimately ran simulations (permutations) rather than paramertic analysis, so I don' think the 'normality' construction was later used.

fella.data <- loadKEGGdata(databaseDir = tmpdir, internalDir = FALSE, loadMatrix = 'diffusion') # using listMethods() here loads all 3 matrices, which in turn allows all 3 analysis methods

getInfo(fella.data) # prints the version KEGG used (Sergio says it always uses the most recent release)
fella.data

# get compounds from WGCNA modules:
mzs = colnames(multiPIC[[set]]$data) 

MOI = "B" # Select module
ModforFella = is.finite(match(modColors[[set]], MOI))
ModforFella
ModforFella<- mzs[ModforFella]

# get KEGG compound names for metabolites:
mz.info <- read.table('correlation among the metabolome/mz.info.txt', header=T, sep=',')
head(mz.info)
rownames(mz.info) <- mz.info$metabolite

module_compounds <- mz.info[ModforFella, ]$KEGG
background_compounds <-  mz.info[mzs, ]$KEGG

c <- defineCompounds(module_compounds, background_compounds, data=fella.data) # prior to running many simulations (permutations), use this line to check that the compounds are on the KEGG graph (not all are usually on the graph, but this line along with the next will tell you if you might have incorrect names)
getExcluded(c) # show compounds that aren't on the KEGG graph

module_analysis <- enrich(compounds = module_compounds, compoundsBackground = background_compounds, data = fella.data, method = 'diffusion', approx = "simulation", niter=10000)
getExcluded(module_analysis) # assoc. compounds were not analyzed because they were not on the KEGG graph

## Results Table and Results Plot,
module_table <- generateResultsTable(object=module_analysis, data=fella.data)

path <- paste0(setLabels[set], '_', MOI, '_module_table.txt')

write.table(module_table, file=path, sep='\t', row.names = F, quote=F)

par(mfrow=c(1,1))
par(mar=c(4,4,4,4))

networkName <- paste(setLabels[set], '-',MOI, 'module')

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




######################################################################################################################################################
# repeat this analysis without virilis.  Analyze the effect of virilis removal on ls ~ mz and all other relationships:
######################################################################################################################################################

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
  write.table(tmpics, file=paste('correlation among the metabolome/pics', k, 'virilis removed'))
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
    plot(pic.cors[[set]][ ,traits.of.interest[[i]]], nv.pic.cors[[set]][ ,traits.of.interest[[i]]], main=paste(setLabels[set], '-', traits.of.interest[[i]]), col=reorderedmycol4[set], pch=16, ylab=paste('virilis removed', '(cor', traits.of.interest[[i]], '~ mzPIC)'), xlab=paste('all species', '(cor', traits.of.interest[[i]], '~ mzPIC)')) 
    legend('bottomright', legend=paste('cor=', round(cor(pic.cors[[set]][ ,traits.of.interest[[i]]], nv.pic.cors[[set]][ ,traits.of.interest[[i]]]), 3)), bty='n')} }

out <- matrix(nrow=4, ncol=3)
for(set in 1:nSets){
  for(i in 1:3){
    out[set, i] <- cor(pic.cors[[set]][ ,traits.of.interest[[i]]], nv.pic.cors[[set]][ ,traits.of.interest[[i]]]) }}

out
colMeans(out)

# of the three lifespan traits, mean lifespan seems to have the smallest effect of virilis on the relationship between PIC(ls) ~ PIC(mz).  I say this because, the correlation between cor(ls ~ mz) for the with and without virilis PICs are themselves correleated with at least cor=0.486 (for the young female metabolome) with an average cor across groups of 0.66 [also, see the plots generated in the loops above].

# I therefore choose to analyze mean lifespan and to note that the ls ~ mz relationships are not broadly effected by removing virilis.
#######################################################################################################################################################################################
