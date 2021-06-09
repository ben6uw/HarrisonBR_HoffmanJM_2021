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


mycol4 = c("#EEC900", "#CD6600", "#87CEFA", "#104E8B")


setwd("/Volumes/GoogleDrive/My Drive/Documents/Hoffman_Metabolite_Phylogeny_Paper/Harrison_Hoffman_GitHub") 

flytree <- read.tree('fly.species.list.nwk')

#############################################################################
# load data from positive and negative modes:
#############################################################################
read.table('Comparative metabolomics.untargeted neg mode logged-only.csv', sep=',', stringsAsFactors = F, header=T) -> neg
read.table('Comparative metabolomics.untargeted pos mode logged-only.csv', sep=',', stringsAsFactors = F, header=T) -> pos

colnames(pos)[6] <- 'sex'
colnames(neg)[6] <- 'sex'

pos[1:4,1:9]
neg[1:4,1:9]
colnames(pos) <- gsub('X', 'pos_', colnames(pos))
colnames(neg) <- gsub('X', 'neg_', colnames(neg))

table(neg$New.Vial.number.for.Metabolomics)
table(pos$New.Vial.number.for.Metabolomics)
table(neg$New.Vial.number.for.Metabolomics %in% pos$New.Vial.number.for.Metabolomics)
neg2 <- neg[neg$New.Vial.number.for.Metabolomics %in% pos$New.Vial.number.for.Metabolomics, ]
match(neg2$New.Vial.number.for.Metabolomics, pos$New.Vial.number.for.Metabolomics)

glob <- cbind(neg2[match(neg2$New.Vial.number.for.Metabolomics, pos$New.Vial.number.for.Metabolomics), ], pos[ ,9:ncol(pos)])

write.table(glob, 'combined.untargeted.data.csv', row.names=F, quote=F, sep=',')
rm(pos)
rm(neg)
rm(neg2)

glob <- read.table('combined.untargeted.data.csv', sep=',', stringsAsFactors = F, header=T)
glob[1:4,1:10]

table(glob$species, glob$sex)


## get means for each strain & sex (this is a means to study variation over species while still including info from each strain):
glob$strain <- paste0(glob$species, '_', glob$strain)
table(glob$strain, glob$sex)
glob[1:4,1:10]

glob.strain.means <- aggregate(glob[, 9:ncol(glob)], list(glob$sex, glob$strain), mean, na.rm=T)
glob.strain.means[1:5, 1:5]
names(glob.strain.means)[1:2] <- c('sex', 'strain')
glob.strain.means[1:5, 1:5]
separate(glob.strain.means, 'strain', into=c('species', 'strain')) -> glob.strain.means
glob.strain.means[1:4, 1:10]
glob.strain.means$strain <- paste(glob.strain.means$species, glob.strain.means$strain, sep='_')

glob.species.means <- aggregate(glob.strain.means[ ,4:ncol(glob.strain.means)], list(glob.strain.means$sex, glob.strain.means$species), mean, na.rm=T)
glob.species.means[1:5, 1:10]
names(glob.species.means)[1:2] <- c('sex', 'species')

table(colSums(is.na(glob.strain.means))[4:ncol(glob.strain.means)]) # 362 features present at least once in each strain AND in BOTH sexes
table(colSums(is.na(glob.strain.means[glob.strain.means$sex == 'females', ]))[4:ncol(glob.strain.means)]) # 554 features present at least once in each strain in females
table(colSums(is.na(glob.strain.means[glob.strain.means$sex == 'males', ]))[4:ncol(glob.strain.means)]) # 609 features present at least once in each strain in males
#############################################################################



#############################################################################
# missingness by sample, species, sex, etc:
#############################################################################
tmp <- glob.strain.means
(rowSums(is.na(tmp)))/ncol(tmp) -> tmp$missing
par(mfrow=c(2,2))
plot(missing ~ as.factor(species), tmp, las=1)
stripchart(missing ~ as.factor(species), tmp, add=T, vertical=T, pch=16)
plot(missing ~ as.factor(sex), tmp, las=1)
stripchart(missing ~ as.factor(sex), tmp, add=T, vertical=T, pch=16)
plot(missing ~ as.factor(strain), tmp, las=1)
stripchart(missing ~ as.factor(strain), tmp, add=T, vertical=T, pch=16)
plot(tmp$missing[order(tmp$missing)], pch=16, ylab='missing', xlab='sample')

par(mfrow=c(1,1))
plot(apply(tmp[ ,4:4421], 2, function(z) 100*(sum(is.na(z))/58)) ~ apply(tmp[ ,4:4421], 2, function(z) mean(z, na.rm=T)), las=1, pch=16, xlab='mean log-metabolite intensity', log='x', ylab='missing data (%)')

missing <- apply(tmp[ ,4:4421], 2, function(z) 100*(sum(is.na(z))/58))
signal <- apply(tmp[ ,4:4421], 2, function(z) mean(z, na.rm=T))


par(mfrow=c(1,2))
plot(signal ~ missing, cex=0.5, pch=16, las=1, ylab='log(LCMS signal)', xlab='# missing')

lin.mod <- lm(signal ~ missing)
summary(lin.mod)
abline(lin.mod, col=2)

pf(578.8, 1, 4416, lower.tail = F) # calculate the exact pvalue from the F distribution

# phylogenetic signal in missingness?
missing <- tmp$missing
names(missing) <- tmp$species
missing
plotTree.boxplot(flytree, missing)
phylosig(flytree, missing, method='K', test=T, nsim=10000) # there is no evidence for phylogenetic signal in missingness
###################################################################


###################################################################
# imputation
###################################################################
tmp <- glob.strain.means[ ,4:ncol(glob.strain.means)]
impute.these <- as.matrix(tmp[ ,colSums(is.na(tmp)) <2])
table(colSums(is.na(impute.these)))
imputed <- knnImputation(impute.these)
dim(imputed) # we get 590 features

dat <- cbind(glob.strain.means[ ,1:3], imputed)
str(dat)
dat$sex <- as.factor(dat$sex)
dat$species <- as.factor(dat$species)
dat$strain <- as.factor(dat$strain)

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
cents <- aggregate(cbind(mean.PC1=PC1, mean.PC2=PC2, mean.PC3=PC3, mean.PC4=PC4) ~ species + sex, PCs, mean)
head(cents)
head(PCs)
cents <- merge(PCs, cents, by=c('species', 'sex')) 
head(cents)

write.table(cents, 'PCA_withcentroids', row.names = F, quote=F)


# make plots of sexual dimorphism within PC1 and PC2:
# example of sex dimorphism plot:
ggplot(cents, aes(PC1 , PC2, color=factor(sex))) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_point(size=1, shape = ifelse(cents$species %in% "ananassae", 19, 1)) +
  geom_segment(aes(x=mean.PC1, mean.PC2, xend=PC1, yend=PC2), data = subset(cents, species %in% c("ananassae"))) +
  geom_segment(aes(x=mean.PC1, mean.PC2, xend=PC1, yend=PC2), data = subset(cents, species %in% tmp)) +
  geom_line(aes(x=mean.PC1, mean.PC2, group=1), data = subset(cents, species %in% c("ananassae")), col='grey', linetype = "dashed") +
  ggtitle(paste('D.', "ananassae")) +
  labs(col = "sex")



# create multi-panel figure
plot_list = list()
for(i in 1:length(unique(PCs$species))){
  unique(PCs$species)[i] -> tmp
  tmp_plot = ggplot(cents, aes(PC1 , PC2, color=factor(sex))) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position = 'none') +
    geom_point(size=ifelse(cents$species %in% tmp, 1.5, 1), shape = ifelse(cents$species %in% tmp, 19, 1)) +
    geom_segment(aes(x=mean.PC1, mean.PC2, xend=PC1, yend=PC2), data = subset(cents, species %in% tmp)) +
    geom_line(aes(x=mean.PC1, mean.PC2, group=1), data = subset(cents, species %in% tmp), col='grey', linetype = "dashed") +
    ggtitle(paste('D.', tmp)) +
    labs(col = "sex")
  
  plot_list[[i]] = tmp_plot
}

# create pdf where each page is a separate plot.
pdf("sexual.dimophism.plots.pdf")
for (i in 1:length(unique(PCs$species))) {
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

x <- as.data.frame(pivot_wider(cents, id_cols = c(strain, species), names_from = sex, values_from= c(PC1, PC2)))
head(x)

PC1_PC2_dimorphism <- sqrt((x$PC1_females - x$PC1_males)^2 + (x$PC2_females - x$PC2_males)^2) # calculate euclidian distance between male and female samples
x$PC1_PC2_dimorphism <- PC1_PC2_dimorphism

plot(PC1_PC2_dimorphism ~ species, x)
names(PC1_PC2_dimorphism) <- x$species

plot(flytree)
nodelabels()
plotTree.boxplot(rotateNodes(flytree, nodes=c(15, 17)), PC1_PC2_dimorphism, args.boxplot=list(xlab="sexual dimophism (PC1, PC2)"))

phylosig(flytree, PC1_PC2_dimorphism, test=T, nsim=1000)

###################################################################################
# species level bar plots
###################################################################################
# phylogenetic signal
###################################################################################
head(PCs)
agg <- aggregate(PCs[ ,4:9], by=list(PCs$species, PCs$sex), mean, na.rm=T) 
head(agg)
names(agg)[1:2] <- c('species', 'sex')
ses <- aggregate(PCs[ ,4:9], by=list(PCs$species, PCs$sex), se, na.rm=T)
names(ses) <- c('species', 'sex', 'sePC1', 'sePC2', 'sePC3', 'sePC4', 'sePC5', 'sePC6')
mean.se <- apply(ses[ ,3:8], 2, mean, na.rm=T)
ses[ses$species == 'erecta' | ses$species == 'yakuba', 3:8] <- mean.se
z <- cbind(agg, ses[ ,3:8])
head(z)

cast <- dcast(setDT(z), species ~ sex, value.var = c("PC1", "sePC1", 'PC2', 'sePC2', 'PC3', 'sePC3'))
cast <- data.frame(cast)
rownames(cast) <- cast$species
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


