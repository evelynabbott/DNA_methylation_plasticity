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


load("/Users/evelynabbott/Dropbox/Mac/Downloads/tagseq_deseq_results.Rdata")
coldata=subset(coldata,coldata$samplingTimepoint == "2")

#to include afternoon comment out next line
coldata=subset(coldata,coldata$time_day == "morning")
coldata$treatment=factor(coldata$treatment,levels=c("heat","control","heatSwitched"))
rownames(coldata)=coldata$sampleNumber

keeps=rownames(coldata)
counts=counts[,keeps]

#to include afternoon add time_day to deseq formula
dds = DESeqDataSetFromMatrix(countData=counts, colData=coldata, 
                             design=~colony+age+treatment)
dds=DESeq(dds)

resultsNames(dds)

cvh=results(dds,name="treatment_control_vs_heat")
cvh=data.frame("cvh"=cvh$log2FoldChange,"gene"=rownames(cvh))

svh=results(dds,name="treatment_heatSwitched_vs_heat")
svh=data.frame("svh"=svh$log2FoldChange,"gene"=rownames(svh))

all=full_join(cvh,svh,by="gene")

#color by basemean GBM --------
load("~/Dropbox/project/Projects/stsw/mdRAD/mdRAD_deseq_results.Rdata")
coldata=subset(coldata,coldata$samplingTimepoint == "1"|coldata$samplingTimepoint == "2")
coldata$time_day <- factor(coldata$time_day, levels = c("morning","afternoon"))
coldata$samplingTimepoint=factor(coldata$samplingTimepoint,levels=c("1","2"))
coldata$treatment=factor(coldata$treatment,levels=c("control","heat","heatSwitched"))

rownames(coldata)=coldata$sampleNumber
keeps=rownames(coldata)
counts=counts[,keeps]
dds = DESeqDataSetFromMatrix(countData=counts, colData=coldata,
                             design=~1)
dds=DESeq(dds)
bm=results(dds,name = "Intercept")
bm = data.frame("bm"=bm$baseMean,"gene"=rownames(bm))

all=full_join(all,bm,by="gene")

#--------------------------------------------------------------------------------------------
#methylation lfc heat

load("~/Dropbox/project/Projects/stsw/mdRAD/mdRAD_deseq_results.Rdata")
coldata=subset(coldata,coldata$samplingTimepoint == "1"|coldata$samplingTimepoint == "2")
coldata$time_day <- factor(coldata$time_day, levels = c("morning","afternoon"))
coldata$samplingTimepoint=factor(coldata$samplingTimepoint,levels=c("1","2"))
coldata=subset(coldata,coldata$treatment == "heat")

rownames(coldata)=coldata$sampleNumber
keeps=rownames(coldata)
counts=counts[,keeps]
dds = DESeqDataSetFromMatrix(countData=counts, colData=coldata,
                             design=~colony+age+time_day+samplingTimepoint)
dds=DESeq(dds)
mh=results(dds,name = "samplingTimepoint_2_vs_1")
mh = data.frame("mh"=mh$log2FoldChange,"gene"=rownames(mh))

all=full_join(all,mh,by="gene")

df=data.frame("cvh"=all$cvh,"svh"=all$svh,"mh"=all$mh,"bm"=all$bm,"gene"=all$gene)
df=na.omit(df)


#load colors
load("~/Dropbox/project/Projects/stsw/mdRAD/ixn_term_investigation/allkME.Rdata")

allkME = allkME %>% select(magenta)
allkME$gene=rownames(allkME)

df=full_join(df,allkME,by="gene")


#ordisurf-----------------------
#basemean
gbm_bm=df$bm
os=ordisurf(df[,1:2]~gbm_bm)
plot(df[,1:2],asp=1,bty="n",col=rgb(0,0,0,0.3))
plot(os,pch=".",col="red")

summary(os)

#ordisurf for ggplot
extract.xyz <- function(obj) {
  xy <- expand.grid(x = obj$grid$x, y = obj$grid$y)
  xyz <- cbind(xy, c(obj$grid$z))
  names(xyz) <- c("x", "y", "z")
  return(xyz)
}

bm.contour.vals <- extract.xyz(obj = os)

#gbm lfc
gbm_lfc=df$mh
os=ordisurf(df[,1:2]~gbm_lfc)
plot(df[,1:2],asp=1,bty="n",col=rgb(0,0,0,0.3))
plot(os,pch=".",col="red")

summary(os)

#ordisurf for ggplot
extract.xyz <- function(obj) {
  xy <- expand.grid(x = obj$grid$x, y = obj$grid$y)
  xyz <- cbind(xy, c(obj$grid$z))
  names(xyz) <- c("x", "y", "z")
  return(xyz)
}

lfc.contour.vals <- extract.xyz(obj = os)


#color palette test for final fig
#hex col = bm; contour = lfc
a=ggplot(df, aes(x=cvh,y=svh,z=log10(bm)))+
  stat_summary_hex()+
  theme_classic()+
  #guides(fill=guide_legend(title="New Legend Title"))
  scale_fill_viridis(option="rocket",name=paste0("GBM", "\n" ,"basemean"))+
  stat_contour(data=bm.contour.vals,aes(x, y, z = z,color=..level..))+
  xlab("GE control vs heat")+
  ylab("GE switch vs heat")+
  theme(legend.position = "none")+
  scale_color_distiller(palette = "Purples",name=paste0("GBM", "\n" ,"basemean"))

s=lm(svh~cvh,data=df)
summary(s)


b=ggplot(df, aes(x=cvh,y=svh,z=gbm_lfc))+
  stat_summary_hex()+
  theme_classic()+
  #guides(fill=guide_legend(title="New Legend Title"))
  scale_fill_viridis(option = "mako", name=paste0("GBM", "\n" ,"lfc"))+
  stat_contour(data=lfc.contour.vals,aes(x, y, z = z,color=..level..))+
  xlab("GE control vs heat")+
  ylab("GE switch vs heat")+
  theme(legend.position = "none")+
  scale_color_distiller(palette = "Greens",name=paste0("GBM", "\n" ,"lfc"))

#grid.arrange(a,b,nrow=1)


c=ggplot(df, aes(x=cvh,y=svh,z=magenta))+
  stat_summary_hex()+
  theme_classic()+
  scale_fill_distiller(palette = "RdPu",name=paste0("magenta", "\n" ,"membership"),trans="reverse",guide=guide_colorbar(reverse = TRUE))+
  #stat_contour(data=bm.contour.vals,aes(x, y, z = z,color=..level..))+
  xlab("GE control vs heat")+
  ylab("GE switch vs heat")+
  theme(legend.position = "none")+
  scale_color_distiller(palette = "Purples",name=paste0("GBM", "\n","basemean"))


grid.arrange(a,b,c,nrow=1)



