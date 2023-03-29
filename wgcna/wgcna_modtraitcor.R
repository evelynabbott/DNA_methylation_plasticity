library(DESeq2)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(tidyverse)
library(dplyr)
library(viridis)
library(gridExtra)
library(RColorBrewer)
library(vegan)


#WGCNA trait correlations 
rm(list=ls())
# Load the WGCNA package
library(WGCNA)
library(tidyverse)
source('~/Dropbox/project/Projects/stsw/mdRAD_wgcna/wgcna_functions.R')
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

#tagseq results
load("~/Dropbox/project/Projects/stsw/mdRAD/ixn_term_investigation/tagseq_wgcna/wgcna_input.Rdata")
load("~/Dropbox/project/Projects/stsw/mdRAD/ixn_term_investigation/tagseq_wgcna/networkdata_signed.RData")


#match up the rownames in the expression 
rownames(datExpr)[1:10]     #this tells us the rownames as given to WGCNA
rownames(datTraits)[1:10]   #this is for the entire set of datTraits and includes samples that were not put into WGCNA
rownames(MEs) = rownames(datExpr)
#datTraits = datTraits[rownames(datTraits) %in% rownames(datExpr),]  #reduce datTraits to those we have WGCNA results for
#if everything is matched up this should say TRUE
sum(rownames(datTraits) == rownames(datExpr)) == length(rownames(datExpr))
head(datTraits)
dim(datTraits)



#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate module eigengenes with color labels
#This will provide the first principal component for expression
#behavior for the genes in each module. On that principal component,
#each sample will have a loading value. These values can be corelated
#with our known sample traits to get an idea what biological mechanism
#the co-reguated genes in a given module might be responding to
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs1 = MEs0[rownames(datTraits),]


########################################
############ SUBSET BY SIZE ############
########################################

#plot all modules
minModSizeToPlot = 0;outputMM=TRUE

# plot only large modules
#minModSizeToPlot = 600;outputMM=FALSE

########################################
########################################
########################################


module.sizes = table(moduleColors)
passing=names(module.sizes)[module.sizes>minModSizeToPlot]
MEs = orderMEs(MEs1[,paste('ME', passing, sep='')])

# 
########################################
########################################
#make a methylation option
#load("~/Dropbox/project/Projects/stsw/mdRAD/ixn_term_investigation/datTraits_meBM.Rdata")
#datTraits=full_join(datTraits,me_sums,by="sampleNumber")


#look at logCcount, CD ratio, minorLR
# datTraits = datTraits %>%
#   select(colony,age,treatment,samplingTimepoint,time_day,quant)

datTraits$c1 = ifelse(datTraits$treatment == "control" & datTraits$samplingTimepoint == "1",1,0)
datTraits$c2 = ifelse(datTraits$treatment == "control" & datTraits$samplingTimepoint == "2",1,0)

datTraits$h1 = ifelse(datTraits$treatment == "heat" & datTraits$samplingTimepoint == "1",1,0)
datTraits$h2 = ifelse(datTraits$treatment == "heat" & datTraits$samplingTimepoint == "2",1,0)

datTraits$s1 = ifelse(datTraits$treatment == "heatSwitched" & datTraits$samplingTimepoint == "1",1,0)
datTraits$s2 = ifelse(datTraits$treatment == "heatSwitched" & datTraits$samplingTimepoint == "2",1,0)

datTraits =datTraits %>%
  select(c1,h1,s1,c2,h2,s2)


#use the cor() function to get the correlations between the module eigengenes and the trait data
moduleTraitCor = cor(MEs, datTraits, use = "p");
#get p values as well
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


#write out module loadings
mEigs=MEs
rownames(mEigs) = rownames(datTraits)
#save(mEigs, file='moduleEigengenes_treatment.Rdata')



#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

#set up additional formatting variables
rows = rownames(moduleTraitCor)
sub.colors = substr(rownames(moduleTraitCor), 3, 50) #trim away the ME from the
module.sizes = paste(sub.colors, table(moduleColors)[sub.colors], sep = "\n")
#dev.off()

#coul <- coul <- rev(colorRampPalette(brewer.pal(8, "RdYlBu"))(25))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = c("control 1","heat 1","switch 1","control 2", "heat 2", "switch 2"),
               yLabels = module.sizes,
               ySymbols = module.sizes,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               plotLegend = FALSE,
               cex.text = 0.7,
               cex.lab.y = .8,
               cex.lab.x = .8,
               zlim = c(-1,1),
               font.lab.y = 1,
               main = paste("Module-treatment relationships"))


