#interaction term against GBM lfc and basemean

library(DESeq2)
library(tidyverse)
library(ggplot2)
library(viridis)
library(vegan)
library(ggpmisc)
library(gridExtra)
#install.packages("ggpointdensity")
#break up by timepoint
rm(list=ls())

#gene expression ixn
load("/Users/evelynabbott/Dropbox/Mac/Downloads/tagseq_deseq_results.Rdata")
coldata=subset(coldata,coldata$samplingTimepoint == "1"|coldata$samplingTimepoint == "2")
coldata$time_day <- factor(coldata$time_day, levels = c("morning","afternoon"))
coldata$samplingTimepoint=factor(coldata$samplingTimepoint,levels=c("1","2"))
coldata$treatment=factor(coldata$treatment,levels=c("control","heat","heatSwitched"))

rownames(coldata)=coldata$sampleNumber
keeps=rownames(coldata)
counts=counts[,keeps]

rownames(coldata)=coldata$sampleNumber
keeps=rownames(coldata)
counts=counts[,keeps]
dds = DESeqDataSetFromMatrix(countData=counts, colData=coldata, 
                             design=~colony+age+treatment+time_day*samplingTimepoint)
dds=DESeq(dds)

resultsNames(dds)
ixn=results(dds,name = "time_dayafternoon.samplingTimepoint2")
ixn = data.frame("ixn"=ixn$log2FoldChange,"gene"=rownames(ixn))

#gbm change --------
load("~/Dropbox/project/Projects/stsw/mdRAD/mdRAD_deseq_results.Rdata")
coldata=subset(coldata,coldata$samplingTimepoint == "1"|coldata$samplingTimepoint == "2")
coldata$time_day <- factor(coldata$time_day, levels = c("morning","afternoon"))
coldata$samplingTimepoint=factor(coldata$samplingTimepoint,levels=c("1","2"))
coldata$treatment=factor(coldata$treatment,levels=c("control","heat","heatSwitched"))

rownames(coldata)=coldata$sampleNumber
coldata=as.data.frame(subset(coldata,coldata$treatment==a))
rownames(coldata)=coldata$sampleNumber
keeps=rownames(coldata)
counts=counts[,keeps]
dds = DESeqDataSetFromMatrix(countData=counts, colData=coldata,
                             design=~colony+age+treatment+time_day+samplingTimepoint)
dds=DESeq(dds)
m=results(dds,name = "samplingTimepoint_2_vs_1")
m = data.frame("m"=m$log2FoldChange,"gene"=rownames(m))

all=full_join(ixn,m,by="gene")


#plot

ggplot(all, aes(x=m,y=ixn))+
  geom_hex()+
  theme_classic()+
  #xlim(-4,4)+
  #ggtitle("control")+
  scale_fill_viridis(trans="log",breaks=c(1,7,54,403))+
  ylab("GE interaction")+
  xlab("GBM change")+
  theme(legend.position = "none")+
  geom_hline(yintercept = 0, linetype = "dashed", color="gray70")+
  geom_vline(xintercept = 0, linetype = "dashed", color="gray70")+
  geom_vline(xintercept=c(0,0), linetype="dotted",color="grey70")+
  geom_hline(yintercept=c(0,0), linetype="dotted",color="grey70")+
  geom_abline(intercept=0,slope = 1,linetype="dotted",color="grey70")
