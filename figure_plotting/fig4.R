library(DESeq2)
library(vegan)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(cowplot)
library(data.table)
library(limma)
library(viridis)
library(hexbin)
library(ggpubr)
rm(list=ls())

#load the stsw data
load("~/Dropbox/project/Projects/stsw/tagseq/tagseq_deseq_results.Rdata")
#coldata1=coldata
coldata=subset(coldata,coldata$samplingTimepoint == "2")
coldata=subset(coldata,coldata$treatment == "control"|coldata$treatment=="heat")
coldata$treatment=factor(coldata$treatment,levels=c("control","heat"))
coldata=subset(coldata,coldata$time_day == "morning")

rownames(coldata)=coldata$sampleNumber
keeps=rownames(coldata)
counts=counts[,keeps]


#find the red mod gene names
m=0.7 #change to min kME value
ll=load("~/Dropbox/even/moduleAssignment.Rdata")
table(moduleColors) # 634 genes in the original red module
# extracting red module genes (only kME>m)
reds=row.names(geneModuleMembership[moduleColors=="red" & geneModuleMembership$MMred>paste0(m),])
reds=intersect(reds,row.names(counts))
length(reds)

#computes the vsd 
dds = DESeqDataSetFromMatrix(countData=counts, colData=coldata,
                             design=~age+colony+treatment)
dds=DESeq(dds)
rld=vst(dds)
rld.df=assay(rld)
colnames(rld.df) = coldata$sampleNumber
vsd = rld.df

# pca before batch effect removal
rawpca=rda(t(vsd)~1)
coldata$rawPC1=scores(rawpca,choices=1,dis="sites",scaling="sites")
coldata$rawPC2=scores(rawpca,choices=2,dis="sites",scaling="sites")

# pca after batch effect removal
vsd=removeBatchEffect(vsd,batch=coldata$colony)
pca1=rda(t(vsd)~1)
coldata$PC1=scores(pca1,choices=1,dis="sites",scaling="sites")
coldata$PC2=scores(pca1,choices=2,dis="sites",scaling="sites")

# plotting PCA versions
ggplot(coldata,aes(PC1,PC2))+geom_point(aes(shape=treatment, color=treatment))+coord_equal()+scale_shape_manual(values=c(1,19))

#compute eigengene
# red module eigengene and PCA
red.vsd=vsd[reds,]
redpc=rda(t(red.vsd)~1)
# importance of axes:
plot(sqrt(redpc$CA$eig/sum(redpc$CA$eig)))
plot(redpc,main="RED",scaling="sites")
coldata$redpc=scores(redpc,choices=1,dis="sites",scaling="sites")
# plotting red module PCA with different coloring
coldata=cbind(coldata,scores=scores(redpc,dis="sites",scaling="sites"))
ggplot(coldata,aes(scores.PC1,scores.PC2,color=treatment))+geom_point()+coord_equal()+theme_classic()

stats=list(c("control","heat"))
ggplot(coldata,aes(treatment,redpc,fill=treatment))+
  theme_classic()+
  geom_boxplot(alpha=0.4)+
  xlab("treatments")+
  theme(axis.text = element_text(size = 12),
        axis.title=element_text(size=14))+
  ylab("redpc")+
  #labs(subtitle = paste("red membership >",m))+
  theme(legend.position = "none")+
  scale_fill_manual(values=c("steelblue2","brown2"))+
  stat_compare_means(comparisons=stats)
