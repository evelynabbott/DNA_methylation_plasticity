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

coldata$lab1=ifelse(grepl("morning",coldata$time_day)&grepl("1",coldata$samplingTimepoint)&grepl("control",coldata$treatment),"cm1",
                    ifelse(grepl("afternoon",coldata$time_day)&grepl("1",coldata$samplingTimepoint)&grepl("control",coldata$treatment),"ca1",
                           ifelse(grepl("morning",coldata$time_day)&grepl("2",coldata$samplingTimepoint)&grepl("control",coldata$treatment),"cm2",
                                  ifelse(grepl("afternoon",coldata$time_day)&grepl("2",coldata$samplingTimepoint)&grepl("control",coldata$treatment),"ca2",
                                         ifelse(grepl("morning",coldata$time_day)&grepl("1",coldata$samplingTimepoint)&grepl("heat",coldata$treatment),"hm1",
                                                ifelse(grepl("afternoon",coldata$time_day)&grepl("1",coldata$samplingTimepoint)&grepl("heat",coldata$treatment),"ha1",
                                                       ifelse(grepl("morning",coldata$time_day)&grepl("2",coldata$samplingTimepoint)&grepl("heat",coldata$treatment),"hm2","ha2")))))))

#coldata=subset(coldata,coldata$samplingTimepoint == "1")
coldata=subset(coldata,coldata$samplingTimepoint == "2")
#coldata=subset(coldata,coldata$samplingTimepoint == "2"|coldata$samplingTimepoint =="1")
coldata=subset(coldata,coldata$treatment == "control"|coldata$treatment=="heat")
#coldata=subset(coldata,coldata$time_day == "afternoon")


coldata=coldata %>% arrange(treatment)
n=coldata$sampleNumber
coldata$sampleNumber=factor(coldata$sampleNumber,levels=n)

rownames(coldata)=coldata$sampleNumber
keeps=rownames(coldata)
counts=counts[,keeps]
dds = DESeqDataSetFromMatrix(countData=counts, colData=coldata, 
                             design=~colony+age+treatment+time_day)

dds = DESeq(dds)
res = results(dds)
rld = vst(dds)
rld.df=assay(rld)
vsd=rld.df
colnames(rld.df) = colnames(counts)

library(limma)
vsd=removeBatchEffect(vsd,batch=coldata$colony)

library(pheatmap)
p=pheatmap(cor(vsd),fontsize = 10, 
           labels_row=coldata$time_day,
           labels_col = coldata$treatment,
           cluster_rows = T,
           cluster_cols=T,
           treeheight_col = F,
           treeheight_row = F)
