# HarrisonBR_HoffmanJM_2021
This GitHub repository includes the data and code used in Modular Evolution of the Drosophila Metabolome (2021), by: Benjamin R. Harrison, Jessica M. Hoffman, Ariana Samuelson, Daniel Raftery, Daniel E.L. Promislow 

As of September 16th 2021, this repository contains R files and the data necessary to run the analysis.  The user will have to download the R files and data into their file system, code the correct file paths and update the required R packages, etc.  

The R files are named:

LC-MS DATA PROCESSING (post peak-calling and metabilte intensity measurement; raw LC-MS data and earlier analytical steps are available at: https://www.metabolomicsworkbench.org/)
targeted.data.processing.R # code to examine and normalize positive ion mode LC-MS/MS (targeted) data
untargeted.neg.data.processing.R # code to examine and normalize positive ion mode LC-TOF-MS (untargeted) data
untargeted.pos.data.processing.R # code to examine and normalize positive ion mode LC-TOF-MS (untargeted) data
combined.untargeted.data.analysis.R # code to combine positive and negative mode untargeted data

ANALYSIS
Hoffman.phylogeny.survival.analysis.R # survival analysis
phylogeny_by mz_targeted_and_global.R # phylogenetic analysis of the LC-MS/MS (targeted) and LC-TOF-MS (untargeted/global) data
metabolome_distance_analysis.R # analysis of metabolome divergence over evolutionary time
MCMCglm_work.R # Markov chain Monte Carlo mixed model analysis of sex x age effects
metabolome PIC network.R # phylogenetically-independent contrasts (PIC) analysis of metabolites and lifespan traits, including covariation


Please contact Ben Harrison: ben6@uw.edu if you need any assistance.
