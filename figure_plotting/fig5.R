rm(list=ls())
library(DESeq2)
library(ggplot2)
library(ggpmisc)
library(tidyverse)
library(dplyr)
library(viridis)
library(gridExtra)
library(RColorBrewer)
library(vegan)
library(grid)

load("~/Dropbox/project/Projects/stsw/mdRAD/mdRAD_geneFPKMs.Rdata")
rownames(geneFpkm)=sub("\\_[^.]*$", "", rownames(geneFpkm))

#get GE lfc -----
load("/Users/evelynabbott/Dropbox/Mac/Downloads/tagseq_deseq_results.Rdata")
coldata=subset(coldata,coldata$samplingTimepoint == "1"|coldata$samplingTimepoint == "2")
coldata$time_day <- factor(coldata$time_day, levels = c("morning","afternoon"))
coldata$samplingTimepoint=factor(coldata$samplingTimepoint,levels=c("1","2"))
coldata$treatment=factor(coldata$treatment,levels=c("control","heat","heatSwitched"))
rownames(coldata)=coldata$sampleNumber
keeps=rownames(coldata)
counts=counts[,keeps]
dds = DESeqDataSetFromMatrix(countData=counts, colData=coldata, 
                             design=~colony+age+treatment+time_day+samplingTimepoint)
dds=DESeq(dds)
resultsNames(dds)
g=results(dds,name="samplingTimepoint_2_vs_1")
g=data.frame("g"=g$log2FoldChange,"gene"=rownames(g))


#mdrad mean--------
load("~/Dropbox/project/Projects/stsw/mdRAD/mdRAD_geneFPKMs.Rdata")
rownames(geneFpkm)=sub("\\_[^.]*$", "", rownames(geneFpkm))
load("/Users/evelynabbott/Dropbox/Mac/Downloads/tagseq_deseq_results.Rdata")

#bimodal distribution
df=data.frame("sums"=rowMeans(geneFpkm))
hist(log(df$sums,2),breaks=100)
a=ggplot(df,aes(x=log(sums,2)))+
  geom_histogram(bins=100)+
  theme_classic()+
  ylab("gene count")+
  xlab("log2(FPKM)")

#bimodal dist with raw expression level 
load("/Users/evelynabbott/Dropbox/Mac/Downloads/tagseq_deseq_results.Rdata")
df1=data.frame("tagseqsum"=rowMeans(counts))
df1$gene=rownames(df1)
df$gene=rownames(df)
df2=full_join(df,df1,by="gene")
df2=full_join(df2,g,by="gene")

df2=na.omit(df2)

grob = grobTree(textGrob(paste("r =", 
                                round(cor(log(df2$sums,2), log(df2$tagseqsum)), 2) ), 
                          x = .05, y = 0.97, hjust = 0, 
                          gp = gpar(col = "black", fontsize = 11)))

b=ggplot(df2,aes(x=log(sums,2),y=log(tagseqsum)))+
  geom_hex(bins=50)+
  scale_fill_viridis(trans="log")+
  geom_smooth(method="lm",formula = y~x,color="blue")+
  theme_classic()+
  xlim(-10,10)+
  xlab("log2(FPKM)")+
  ylab("expression level")+
  theme(legend.position = "none")+
  annotation_custom(grob)
  # stat_fit_glance(method = "lm",
  #                 label.y = "bottom",
  #                 method.args = list(formula = y ~ x),
  #                 mapping = aes(label = sprintf('r^2~"="~%.3f~~italic(p)~"="~%.3g',
  #                                               stat(r.squared), stat(p.value))),
  #                 parse = TRUE)


grob1 = grobTree(textGrob(paste("r =", 
                                round(cor(log(df2$sums,2), abs(df2$g)), 2) ), 
                          x = .05, y = 0.97, hjust = 0, 
                          gp = gpar(col = "black", fontsize = 11)))


c=ggplot(df2,aes(x=log(sums,2),y=abs(g)))+
  #geom_point(alpha=0.2)+
  geom_hex(bins=50)+
  scale_fill_viridis(trans="log")+
  geom_smooth(method="lm",formula = y~x,color="blue")+
  theme_classic()+
  xlim(-10,10)+
  ylim(0,1.5)+
  xlab("log2(FPKM)")+
  ylab("abs. expression difference")+
  theme(legend.position = "none")+
  annotation_custom(grob1)
  # stat_fit_glance(method = "lm",
  #                 label.y = "bottom",
  #                 method.args = list(formula = y ~ x),
  #                 mapping = aes(label = sprintf('r^2~"="~%.3f~~italic(p)~"="~%.3g',
  #                                               stat(r.squared), stat(p.value))),
  #                 parse = TRUE)

grid.arrange(a,b,c,nrow=1)


  


