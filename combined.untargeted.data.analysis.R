if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")

library(plyr) # data wrangling
library(splitstackshape)
library(dplyr)
library(tidyr)
library(reshape2)
library(data.table)
library(DMwR) 

library(sjstats) # mixed modeling tools
library(lme4)
library(lmerTest)
library(nlme) 
library(MuMIn)
library(ggpubr) 

library(graphics) # plotting
library(ggplot2) 
library(gplots) 
library(corrplot)

library(ape) # phylogenetic analysis
library(phytools) 
library(geiger) 
library(caper) 
library(dendextend)
library(AssocTests)

gc()
mycol4 = c("#EEC900", "#CD6600", "#87CEFA", "#104E8B")

setwd("/Volumes/GoogleDrive/My Drive/Documents/Hoffman_Metabolite_Phylogeny_Paper/Harrison_Hoffman_GitHub") 

flytree <- read.tree('fly.Species.list.nwk')

#############################################################################
# load data from positive and negative modes:
#############################################################################
read.table('Comparative metabolomics.untargeted neg mode logged-only.csv', sep=',', stringsAsFactors = F, header=T) -> neg
read.table('Comparative metabolomics.untargeted pos mode logged-only.csv', sep=',', stringsAsFactors = F, header=T) -> pos

pos[1:4,1:10]
neg[1:4,1:10]

colnames(pos) <- gsub('X', 'pos_', colnames(pos))
colnames(neg) <- gsub('X', 'neg_', colnames(neg))

neg2 <- neg[neg$New.Vial.number.for.Metabolomics %in% pos$New.Vial.number.for.Metabolomics, ]
match(neg2$New.Vial.number.for.Metabolomics, pos$New.Vial.number.for.Metabolomics)

glob <- cbind(neg2[match(neg2$New.Vial.number.for.Metabolomics, pos$New.Vial.number.for.Metabolomics), ], pos[ ,10:ncol(pos)])

glob <- subset(glob, select=-c(Use)) # remove the 'Use' column

write.table(glob, 'combined.untargeted.data.csv', row.names=F, quote=F, sep=',')
rm(pos)
rm(neg)
rm(neg2)

glob[1:4,1:10]

## get means for each Strain & Sex (this is a means to study variation over Species while still including info from each Strain):

x <- glob$Line
x <- strsplit(x, '-')
glob$Line <- paste0(glob$Species, '_', sapply(x, "[[", 2))
colnames(glob)[colnames(glob) == 'Line'] <- 'Strain'
glob[1:4,1:10]

table(glob$Strain, glob$Sex)

glob.Strain.means <- aggregate(glob[, 9:ncol(glob)], list(glob$Sex, glob$Strain), mean, na.rm=T)
glob.Strain.means[1:5, 1:5]
names(glob.Strain.means)[1:2] <- c('Sex', 'Strain')
separate(glob.Strain.means, 'Strain', into=c('Species', 'Strain')) -> glob.Strain.means
glob.Strain.means[1:4, 1:10]
glob.Strain.means$Strain <- paste(glob.Strain.means$Species, glob.Strain.means$Strain, sep='_')

glob.Species.means <- aggregate(glob.Strain.means[ ,4:ncol(glob.Strain.means)], list(glob.Strain.means$Sex, glob.Strain.means$Species), mean, na.rm=T)
glob.Species.means[1:5, 1:10]
names(glob.Species.means)[1:2] <- c('Sex', 'Species')

table(colSums(is.na(glob.Strain.means))[4:ncol(glob.Strain.means)]) # 362 features present at least once in each Strain AND in BOTH Sexes
table(colSums(is.na(glob.Strain.means[glob.Strain.means$Sex == 'females', ]))[4:ncol(glob.Strain.means)]) # 554 features present at least once in each Strain in females
table(colSums(is.na(glob.Strain.means[glob.Strain.means$Sex == 'males', ]))[4:ncol(glob.Strain.means)]) # 609 features present at least once in each Strain in males
#############################################################################

#############################################################################
# missingness by sample, Species, Sex, etc:
#############################################################################
tmp <- glob.Strain.means
(rowSums(is.na(tmp)))/ncol(tmp) -> tmp$missing
par(mfrow=c(2,2))
plot(missing ~ as.factor(Species), tmp, las=1)
stripchart(missing ~ as.factor(Species), tmp, add=T, vertical=T, pch=16)
plot(missing ~ as.factor(Sex), tmp, las=1)
stripchart(missing ~ as.factor(Sex), tmp, add=T, vertical=T, pch=16)
plot(missing ~ as.factor(Strain), tmp, las=1)
stripchart(missing ~ as.factor(Strain), tmp, add=T, vertical=T, pch=16)
plot(tmp$missing[order(tmp$missing)], pch=16, ylab='missing', xlab='sample')

tmp[1:4,1:4]

missing <- apply(tmp[ ,4:4421], 2, function(z) 100*(sum(is.na(z))/58))
signal <- apply(tmp[ ,4:4421], 2, function(z) mean(z, na.rm=T))

par(mfrow=c(2,2))
plot(log(signal) ~ missing, cex=0.5, pch=16, las=1, ylab='log(LCMS signal)', xlab='# missing', main='major axis regression')
lin.mod2 <- ma(log(signal) ~ missing)
summary(lin.mod2)
abline(lin.mod2, col=2, lwd=2)

plot(missing ~ log(signal), cex=0.5, pch=16, las=1, xlab='log(LCMS signal)', ylab='# missing', main='major axis regression')
lin.mod1 <- ma(missing ~ log(signal))
summary(lin.mod1)
abline(lin.mod1, col=2, lwd=2)


pf(596.7, 1, 4416, lower.tail = F) # calculate the exact pvalue from the F distribution

# phylogenetic signal in missingness?
missing <- tmp$missing
names(missing) <- tmp$Species
plotTree.boxplot(flytree, missing)
phylosig(flytree, missing, method='K', test=T, nsim=10000) # there is no evidence for phylogenetic signal in missingness
###################################################################


###################################################################
# imputation
###################################################################
tmp <- glob.Strain.means[ ,4:ncol(glob.Strain.means)]
impute.these <- as.matrix(tmp[ ,colSums(is.na(tmp)) <2])
table(colSums(is.na(impute.these)))
imputed <- knnImputation(impute.these)
dim(imputed) # we get 590 features

dat <- cbind(glob.Strain.means[ ,1:3], imputed)
str(dat)
dat$Sex <- as.factor(dat$Sex)
dat$Species <- as.factor(dat$Species)
dat$Strain <- as.factor(dat$Strain)

dat[1:4,1:10]

centered_dat <- t(scale(t(dat[ ,4:ncol(dat)]), scale=F)) # center, but don't scale metabolites by sample
dat <- cbind(dat[ ,1:3], centered_dat)

write.table(dat, 'combined_modes_untargeted_inputed_590.features', row.names=F, quote=F)

dat <- read.table('combined_modes_untargeted_inputed_590.features', header=T, stringsAsFactors = T)
dat[1:4,1:5]

###################################################################
### PCA with the untargeted metabolome:
###################################################################
dev.off()
glob.PCA <- prcomp(dat[ ,4:ncol(dat)], scale=T)
tw(glob.PCA$sdev^2, 58, criticalpoint = 2.0234)$SigntEigenL #NOTE this critical point = alpha 0.01, see ?tw for other values
plot(glob.PCA)

PCs <- cbind(dat[ ,1:3], glob.PCA$x[ ,1:6])
head(PCs)
cents <- aggregate(cbind(mean.PC1=PC1, mean.PC2=PC2, mean.PC3=PC3, mean.PC4=PC4) ~ Species + Sex, PCs, mean)
head(cents)
head(PCs)
cents <- merge(PCs, cents, by=c('Species', 'Sex')) 
head(cents)

write.table(cents, 'PCA_with_centroids', row.names = F, quote=F)


# make plots of Sexual dimorphism within PC1 and PC2:
# example of Sex dimorphism plot:
ggplot(cents, aes(PC1 , PC2, color=factor(Sex))) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_point(size=1, shape = ifelse(cents$Species %in% "ananassae", 19, 1)) +
  geom_segment(aes(x=mean.PC1, mean.PC2, xend=PC1, yend=PC2), data = subset(cents, Species %in% c("ananassae"))) +
  geom_segment(aes(x=mean.PC1, mean.PC2, xend=PC1, yend=PC2), data = subset(cents, Species %in% tmp)) +
  geom_line(aes(x=mean.PC1, mean.PC2, group=1), data = subset(cents, Species %in% c("ananassae")), col='grey', linetype = "dashed") +
  ggtitle(paste('D.', "ananassae")) +
  labs(col = "Sex")



# create multi-panel figure
plot_list = list()
for(i in 1:length(unique(PCs$Species))){
  unique(PCs$Species)[i] -> tmp
  tmp_plot = ggplot(cents, aes(PC1 , PC2, color=factor(Sex))) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position = 'none') +
    geom_point(size=ifelse(cents$Species %in% tmp, 1.5, 1), shape = ifelse(cents$Species %in% tmp, 19, 1)) +
    geom_segment(aes(x=mean.PC1, mean.PC2, xend=PC1, yend=PC2), data = subset(cents, Species %in% tmp)) +
    geom_line(aes(x=mean.PC1, mean.PC2, group=1), data = subset(cents, Species %in% tmp), col='grey', linetype = "dashed") +
    ggtitle(paste('D.', tmp)) +
    labs(col = "Sex")
  
  plot_list[[i]] = tmp_plot
}

# create pdf where each page is a separate plot.
pdf("Sexual.dimophism.plots.pdf")
for (i in 1:length(unique(PCs$Species))) {
  print(plot_list[[i]])
}
dev.off()

# create multipanel plot:
ggarrange(plot_list[[1]], 
          plot_list[[2]], 
          plot_list[[3]], 
          plot_list[[4]],
          plot_list[[5]],
          plot_list[[6]],
          plot_list[[7]],
          plot_list[[8]],
          plot_list[[9]],
          plot_list[[10]],
          plot_list[[11]],
          ncol = 3, nrow = 4)


head(cents)

x <- as.data.frame(pivot_wider(cents, id_cols = c(Strain, Species), names_from = Sex, values_from= c(PC1, PC2)))
head(x)

PC1_PC2_dimorphism <- sqrt((x$PC1_females - x$PC1_males)^2 + (x$PC2_females - x$PC2_males)^2) # calculate euclidian distance between male and female samples
x$PC1_PC2_dimorphism <- PC1_PC2_dimorphism

plot(PC1_PC2_dimorphism ~ Species, x)
names(PC1_PC2_dimorphism) <- x$Species

plot(flytree)
nodelabels()
plotTree.boxplot(rotateNodes(flytree, nodes=c(15, 17)), PC1_PC2_dimorphism, args.boxplot=list(xlab="Sexual dimophism (PC1, PC2)"))

phylosig(flytree, PC1_PC2_dimorphism, test=T, nsim=1000)

###################################################################################
# Species level bar plots
###################################################################################
# phylogenetic signal
###################################################################################
head(PCs)
agg <- aggregate(PCs[ ,4:9], by=list(PCs$Species, PCs$Sex), mean, na.rm=T) 
head(agg)
names(agg)[1:2] <- c('Species', 'Sex')
ses <- aggregate(PCs[ ,4:9], by=list(PCs$Species, PCs$Sex), se, na.rm=T)
names(ses) <- c('Species', 'Sex', 'sePC1', 'sePC2', 'sePC3', 'sePC4', 'sePC5', 'sePC6')
mean.se <- apply(ses[ ,3:8], 2, mean, na.rm=T)
ses[ses$Species == 'erecta' | ses$Species == 'yakuba', 3:8] <- mean.se
z <- cbind(agg, ses[ ,3:8])
head(z)

cast <- dcast(setDT(z), Species ~ Sex, value.var = c("PC1", "sePC1", 'PC2', 'sePC2', 'PC3', 'sePC3'))
cast <- data.frame(cast)
rownames(cast) <- cast$Species
head(cast)

cast <- cast[flytree$tip.label, ]
head(cast)
plotTree.barplot(flytree, cast[ ,c(2,3)], args.barplot = list(col = mycol4[c(2, 4)], beside = T, xlab='PC1 (20.0%)', legend.text = F))
phylosig(flytree, cast[ ,2], method='K', test=T, nsim=10000, se=cast[ ,4]) 
phylosig(flytree, cast[ ,3], method='K', test=T, nsim=10000, se=cast[ ,5]) 

plotTree.barplot(flytree, cast[ ,c(6,7)], args.barplot = list(col = mycol4[c(2, 4)], beside = T, xlab='PC2 (10.4%)', legend.text = F))
phylosig(flytree, cast[ ,6], method='K', test=T, nsim=10000, se=cast[ ,4]) 
phylosig(flytree, cast[ ,7], method='K', test=T, nsim=10000, se=cast[ ,5]) 

plotTree.barplot(flytree, cast[ ,c(10, 11)], args.barplot = list(col = mycol4[c(2, 4)], beside = T, xlab='PC3 (7.4%)', legend.text = F))
phylosig(flytree, cast[ ,10], method='K', test=T, nsim=10000, se=cast[ ,4]) 
phylosig(flytree, cast[ ,11], method='K', test=T, nsim=10000, se=cast[ ,5]) 
###############################################################################################


