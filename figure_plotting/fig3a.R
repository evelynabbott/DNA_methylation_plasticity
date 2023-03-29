#fix the heatmap

library(DESeq2)
library(tidyverse)
library(ggplot2)
library(viridis)
library(vegan)
library(ggpmisc)
library(gridExtra)
rm(list=ls())

load("/Users/evelynabbott/Dropbox/Mac/Downloads/tagseq_deseq_results.Rdata")
coldata=subset(coldata,coldata$samplingTimepoint == "2")
coldata=subset(coldata,coldata$treatment == "control"|coldata$treatment=="heat")
coldata=subset(coldata,coldata$time_day == "morning")

coldata=coldata %>% arrange(treatment)
n=coldata$sampleNumber
coldata$sampleNumber=factor(coldata$sampleNumber,levels=n)

rownames(coldata)=coldata$sampleNumber
keeps=rownames(coldata)
counts=counts[,keeps]
dds = DESeqDataSetFromMatrix(countData=counts, colData=coldata, 
                             design=~colony+age+treatment)

dds = DESeq(dds)
res = results(dds)
rld = vst(dds)
rld.df=assay(rld)
vsd=rld.df
colnames(rld.df) = colnames(counts)

library(limma)
vsd=removeBatchEffect(vsd,batch=coldata$colony)

library(pheatmap)

pheatmap(cor(vsd),fontsize = 10, 
         labels_row=coldata$treatment,
         labels_col = coldata$treatment,
         cluster_rows = F,
         cluster_cols=F)
               
                 



