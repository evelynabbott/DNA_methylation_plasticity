#module GBM
rm(list=ls())
library(DESeq2)
library(ggplot2)
library(ggpmisc)
library(tidyverse)
library(dplyr)
library(viridis)
library(gridExtra)
library(RColorBrewer)
library(patchwork)
library(magrittr)
library(ggprism)


load("~/Dropbox/project/Projects/stsw/mdRAD/ixn_term_investigation/tagseq_wgcna/tagseq_genecolors.Rdata")
genecolors$color <- factor(genecolors$color, levels = c("magenta","black","pink","yellow","brown","greenyellow","blue","salmon","grey60"))

load("~/Dropbox/project/Projects/stsw/mdRAD/mdRAD_deseq_results.Rdata")
coldata=subset(coldata,coldata$samplingTimepoint == "1"|coldata$samplingTimepoint == "2")
coldata$time_day <- factor(coldata$time_day, levels = c("morning","afternoon"))
coldata$samplingTimepoint=factor(coldata$samplingTimepoint,levels=c("1","2"))
coldata$treatment=factor(coldata$treatment,levels=c("heat","control","heatSwitched"))

rownames(coldata)=coldata$sampleNumber
keeps=rownames(coldata)
counts=counts[,keeps]
dds = DESeqDataSetFromMatrix(countData=counts, colData=coldata,
                             design=~colony+age+time_day+treatment+samplingTimepoint)
dds=DESeq(dds)
resultsNames(dds)
m=results(dds,name = "samplingTimepoint_2_vs_1")
m = data.frame("m"=m$log2FoldChange,"gene"=rownames(m))
bm=results(dds,name = "Intercept")
bm = data.frame("bm"=bm$baseMean,"gene"=rownames(bm))

all=full_join(m,genecolors,by="gene")
all=full_join(bm,all)

#make dataframe comparing magenta to all -------
all=na.omit(all) #15404 
df=all
df$color= ifelse(grepl("magenta",df$color),"magenta","all")

#what are the means?
dfa=subset(df,df$color=="all")
dfm=subset(df,df$color=="magenta")
mean(dfa$m)
mean(dfm$m)

#MAGENTA ~ ALL --------
#stats for m
result <- t.test(m ~ color, data = df)$p.value
result <- signif(result, digits = 3)
result
df_p_val <- data.frame(
  group1 = "all",
  group2 = "magenta",
  label = result,
  y.position = 6
)

a=ggplot(df, aes(x=color, y=m)) +
  geom_violin(aes(fill=color),alpha=0.2)+
  #scale_fill_manual(values=c("magenta","black","pink","yellow","brown","greenyellow","blue","salmon","grey60")) +
  scale_fill_manual(values=c("black","magenta"))+
  ggtitle("")+
  theme_classic()+
  theme(legend.position = "none")+
  geom_hline(yintercept = 0,linetype="dotted")+
  #ylim(-.3,.3)+
  ylab("GBM change")+
  xlab("")+
  #theme(axis.text.x = element_blank())+
  theme(axis.text.x=element_text(size=15))+
  theme(axis.title.y=element_text(size=15))+
  stat_summary(fun="mean",geom="point", size=2, color="black",alpha=0.6)+
  geom_text(data=am, aes(label=bm,y=20))+
  add_pvalue(df_p_val,
             label = "p = {label}",
             remove.bracket = TRUE,
             x = 1.5,
             y=1.3,
             label.size = 4.5)


#stats for bm
result <- t.test(bm ~ color, data = df)$p.value
result <- signif(result, digits = 3)
result
df_p_val1 <- data.frame(
  group1 = "all",
  group2 = "magenta",
  label = result,
  y.position = 6
)

meanbm=aggregate(bm~color,df,mean)
meanbm = meanbm %>% 
  mutate_if(is.numeric,
            round,
            digits=2)
ggplot(df, aes(x=color, y=bm)) +
  geom_violin(aes(fill=color),alpha=0.2)+
  #scale_fill_manual(values=c("magenta","black","pink","yellow","brown","greenyellow","blue","salmon","grey60")) +
  scale_fill_manual(values=c("black","magenta"))+
  ggtitle("")+
  theme_classic()+
  theme(legend.position = "none")+
  geom_hline(yintercept = 0,linetype="dotted")+
  ylim(0,200)+
  ylab("GBM basemean")+
  xlab("")+
  theme(axis.text.x=element_text(size=15))+
  theme(axis.title.y=element_text(size=15))+
  stat_summary(fun="mean",geom="point", size=2, color="black",alpha=0.6)+
  geom_text(data=meanbm, aes(label=bm,y=20))+
  add_pvalue(df_p_val1,
             label = "p = {label}",
             remove.bracket = TRUE,
             x = 1.5,
             y=170,
             label.size = 4.5)

meanbm=aggregate(m~color,df,mean)
meanbm = meanbm %>% 
  mutate_if(is.numeric,
            round,
            digits=2)
ggplot(df, aes(x=color, y=m)) +
  geom_violin(aes(fill=color),alpha=0.2)+
  #scale_fill_manual(values=c("magenta","black","pink","yellow","brown","greenyellow","blue","salmon","grey60")) +
  scale_fill_manual(values=c("black","magenta"))+
  ggtitle("")+
  theme_classic()+
  theme(legend.position = "none")+
  geom_hline(yintercept = 0,linetype="dotted")+
  ylim(-1,1)+
  ylab("GBM change")+
  xlab("")+
  theme(axis.text.x=element_text(size=15))+
  theme(axis.title.y=element_text(size=15))+
  stat_summary(fun="mean",geom="point", size=2, color="black",alpha=0.6)+
  geom_text(data=meanbm, aes(label=bm,y=20))+
  add_pvalue(df_p_val1,
             label = "p = {label}",
             remove.bracket = TRUE,
             x = 1.5,
             y=170,
             label.size = 4.5)
