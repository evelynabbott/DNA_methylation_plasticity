library(WGCNA)
library(flashClust)
#install.packages("flashClust")
library(ape)
options(stringsAsFactors=FALSE)
allowWGCNAThreads()
library(DESeq2)
library(ggplot2)
library(vegan)
library(tidyverse)
#install.packages("ggplot2",dependencies = T)
#WGCNA data
#load("~/Desktop/wgcna_input.Rdata")


powers = c(seq(from = 2, to=26, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5,networkType="signed")
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

s.th=12
adjacency = adjacency(datExpr, power = s.th,type="signed");
TOM = TOMsimilarity(adjacency,TOMType="signed");
dissTOM = 1-TOM
# Call the hierarchical clustering function

geneTree = flashClust(as.dist(dissTOM), method = "average");

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30; 
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
METree = flashClust(as.dist(MEDiss), method = "average");

#merge modules
MEDissThres = 0.4 # in the first pass, set this to 0 - no merging (we want to see the module-traits heatmap first, then decide which modules are telling us the same story and better be merged)
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")  # on 2nd pass: does this cut height meet your merging goals? If not, reset MEDissThres and replot

# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

# plotting the fabulous ridiculogram
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = FALSE, guideHang = 0.05,lwd=0.3)

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = flashClust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# how many genes in each module?
table(moduleColors)
# Save module colors and labels for use in subsequent parts
save(MEs, geneTree, moduleLabels, moduleColors, file = "~/Desktop/N4_networkdata_signed.RData")

###################
# plotting correlations with traits:
# load(file = "networkdata_signed.RData")
# load(file = "wgcnaData.RData")
# 
# # Define numbers of genes and samples
# nGenes = ncol(datt);
# nSamples = nrow(datt);
# # Recalculate MEs with color labels
# MEs0 = moduleEigengenes(datt, moduleColors)$eigengenes
# MEs = orderMEs(MEs0)
# 
# # correlations of genes with eigengenes
# moduleGeneCor=cor(MEs,datt)
# moduleGenePvalue = corPvalueStudent(moduleGeneCor, nSamples);
# 
# moduleTraitCor = cor(MEs, datTraits, use = "p");
# moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
# 
# # module-trait correlations
# textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
#                    signif(moduleTraitPvalue, 1), ")", sep = "");
# dim(textMatrix) = dim(moduleTraitCor)
# par(mar = c(6, 8.5, 3, 3));
# # Display the correlation values within a heatmap plot
# hm1<-labeledHeatmap(Matrix = moduleTraitCor,
#                     xLabels = names(datTraits),
#                     yLabels = names(MEs),
#                     ySymbols = names(MEs),
#                     colorLabels = FALSE,
#                     colors = blueWhiteRed(50),
#                     textMatrix = textMatrix,
#                     setStdMargins = FALSE,
#                     cex.text = 0.5,
#                     zlim = c(-1,1),
#                     main = paste("Module-trait relationships"))
# 
# print(data.frame(table(moduleColors))) # gives numbers of genes in each module
# 
# ######### Evelyn plots (heatmap of modules v. columns)
# 
# library(pheatmap)
# library(RColorBrewer)
# #Build heatmap of MEs - samples as rows and MEs as columns
# MEsm = as.matrix(MEs)
# #make rownames in order: primary = cluster, secondary = age
# datTraits1 <- datTraits[order(datTraits$C4,datTraits$C3,datTraits$C2,datTraits$C1),]
# MEsm = MEsm[match(rownames(datTraits), rownames(MEsm)),]
# 
# coul <- coul <- rev(colorRampPalette(brewer.pal(8, "RdYlBu"))(25))
# pheatmap(t(MEsm),scale="row",cluster_cols=FALSE, color = coul,
#          show_colnames = TRUE,show_rownames = TRUE,fontsize = 12)
# #pheatmap(cor(datExpr))
# #make rownames in order: primary = sample location, secondary = age
# datTraits2 <- datTraits[order(datTraits$D,datTraits$N,datTraits$O),]
# MEsm2 = MEsm[match(rownames(datTraits2), rownames(MEsm)),]
# 
# coul <- coul <- rev(colorRampPalette(brewer.pal(8, "RdYlBu"))(25))
# pheatmap(t(MEsm2),scale="row",cluster_cols=FALSE, color = coul,
#          show_colnames = TRUE,show_rownames = TRUE,fontsize = 12)
# 
# #############
# # scatterplots of gene significance (correlation-based) vs kME
# 
# load(file = "networkdata_signed.RData")
# load(file = "wgcnaData.RData")
# load(file = "wgcna_input_age.Rdata")
# traits=datTraits
# traits
# table(moduleColors)
# whichTrait="J4"
# 
# nGenes = ncol(datt);
# nSamples = nrow(datt);
# selTrait = as.data.frame(traits[,whichTrait]);
# names(selTrait) = whichTrait
# # names (colors) of the modules
# modNames = substring(names(MEs), 3)
# geneModuleMembership = as.data.frame(signedKME(datt, MEs));
# MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
# names(geneModuleMembership) = paste("MM", modNames, sep="");
# names(MMPvalue) = paste("p.MM", modNames, sep="");
# geneTraitSignificance = as.data.frame(cor(datt, selTrait, use = "p"));
# GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
# names(geneTraitSignificance) = paste("GS.", names(selTrait), sep="");
# names(GSPvalue) = paste("p.GS.", names(selTrait), sep="");
# #quartz()
# #X11(type = "cairo")
# dev.off()
# par(mfrow=c(3,3))
# counter=0
# 
# for(module in modNames[1:length(modNames)]){
#   counter=counter+1
#   if (counter>9) {
#     #quartz()
#     par(mfrow=c(3,3))
#     counter=1
#   }
#   column = match(module, modNames);
#   moduleGenes = moduleColors==module;
#   #trr="heat resistance"
#   verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
#                      abs(geneTraitSignificance[moduleGenes, 1]),
#                      xlab = paste(module,"module membership"),
#                      ylab = paste("GS for", whichTrait),
#                      col = "grey50",mgp=c(2.3,1,0))
# }
# 
# ################
# # eigengene-heatmap plot (sanity check - is the whole module driven by just one crazy sample?)
# # note: this part does not make much sense for unsigned modules
# 
# load(file = "networkdata_signed.RData")
# load(file = "wgcnaData.RData");
# 
# which.module="darkturquoise" 
# datME=MEs
# datExpr=datt
# #quartz()
# ME=datME[, paste("ME",which.module, sep="")]
# par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
# plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),
#         nrgcols=30,rlabels=F,rcols=which.module,
#         main=which.module, cex.main=2)
# par(mar=c(5, 4.2, 0, 0.7))
# barplot(ME, col=which.module, main="", cex.main=2,
#         ylab="eigengene expression",xlab="sample")
# 
# length(datExpr[1,moduleColors==which.module ]) # number of genes in chosen module
# 
# #################
# # saving selected modules for GO and KOG analysis (two-parts: Fisher test, MWU test within-module)
# 
# library(WGCNA)
# load(file = "networkdata_signed.RData") # moduleColors, MEs
# load(file = "wgcnaData.RData") # vsd table
# vsd=t(datt)
# #load(file = "../data4wgcna.RData") # vsd table
# 
# # calculating modul memberships for all genes for all modules
# allkME =as.data.frame(signedKME(datt, MEs)) 
# names(allkME)=gsub("kME","",names(allkME))
# 
# whichModule="black"
# table(moduleColors==whichModule) # how many genes are in it?
# vsd.wg=t(datt) #<<<< add this line!
# 
# # Saving data for Fisher-MWU combo test (GO_MWU)
# inModuleBinary=as.numeric(moduleColors==whichModule)
# combo=data.frame("gene"=row.names(vsd.wg),"Fish_kME"=allkME[,whichModule]*inModuleBinary)
# write.csv(combo,file=paste(whichModule,".csv",sep=""),row.names=F,quote=F)
# 
# ################
# # plotting heatmap for named top-kME genes
# 
# load(file = "networkdata_signed.RData")
# ll=load(file = "wgcnaData.RData");
# 
# allkME =as.data.frame(signedKME(datt, MEs))
# gg=read.table("amil_iso2gene.tab",sep="\t")
# library(pheatmap)
# 
# whichModule="black"
# top=30 # number of named top-kME genes to plot
# 
# datME=MEs
# datExpr=datt
# modcol=paste("kME",whichModule,sep="")
# sorted=vsd[order(allkME[,modcol],decreasing=T),]
# head(sorted)
# # selection top N named genes, attaching gene names
# gnames=c();counts=0;hubs=c()
# for(i in 1:length(sorted[,1])) {
#   if (row.names(sorted)[i] %in% gg$V1) { 
#     counts=counts+1
#     gn=gg[gg$V1==row.names(sorted)[i],2]
#     gn=paste(gn,row.names(sorted)[i],sep=".")
#     if (gn %in% gnames) {
#       gn=paste(gn,counts,sep=".")
#     }
#     gnames=append(gnames,gn) 
#     hubs=data.frame(rbind(hubs,sorted[i,]))
#     if (counts==top) {break}
#   }
# } 
# row.names(hubs)=gnames
# 
# contrasting = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
# contrasting2 = colorRampPalette(rev(c("chocolate1","chocolate1","#FEE090","grey10", "cyan3","cyan","cyan")))(100)
# contrasting3 = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan","cyan")))(100)
# 
# pheatmap(hubs,scale="row",col=contrasting,border_color=NA,treeheight_col=0.1,cex=0.9,cluster_rows=F)
# 
