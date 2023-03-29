#DAPC & PCA--------

library(DESeq2)
library(tidyverse)
library(ggplot2)
library(viridis)
library(vegan)
library(ggpmisc)
library(gridExtra)
rm(list=ls())

load("~/Desktop/tagseq_counts_coldata.Rdata")
coldata=subset(coldata,coldata$samplingTimepoint == "2"|coldata$samplingTimepoint=="1")
coldata=subset(coldata,coldata$treatment == "control"|coldata$treatment=="heat")

coldata$lab1=ifelse(grepl("morning",coldata$time_day)&grepl("1",coldata$samplingTimepoint)&grepl("control",coldata$treatment),"cm1",
                    ifelse(grepl("afternoon",coldata$time_day)&grepl("1",coldata$samplingTimepoint)&grepl("control",coldata$treatment),"ca1",
                           ifelse(grepl("morning",coldata$time_day)&grepl("2",coldata$samplingTimepoint)&grepl("control",coldata$treatment),"cm2",
                                  ifelse(grepl("afternoon",coldata$time_day)&grepl("2",coldata$samplingTimepoint)&grepl("control",coldata$treatment),"ca2",
                                         ifelse(grepl("morning",coldata$time_day)&grepl("1",coldata$samplingTimepoint)&grepl("heat",coldata$treatment),"hm1",
                                                ifelse(grepl("afternoon",coldata$time_day)&grepl("1",coldata$samplingTimepoint)&grepl("heat",coldata$treatment),"ha1",
                                                       ifelse(grepl("morning",coldata$time_day)&grepl("2",coldata$samplingTimepoint)&grepl("heat",coldata$treatment),"hm2","ha2")))))))

#coldata=subset(coldata,coldata$time_day == "afternoon")
coldata=subset(coldata,coldata$samplingTimepoint=="1")

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

datExpr=t(vsd)


#PCA---------------------------
library(vegan)
dds.pca=capscale(dist(datExpr,method = "manhattan")~coldata$lab1)
plot(dds.pca)
scores.vegan2 <- as.data.frame(dds.pca$CA$u)

ggplot(scores.vegan2, aes(MDS1, MDS2, color = coldata$lab1))+
  #stat_ellipse()+
  geom_point(size=3)+
  theme_classic()+
  #coord_fixed()+
  #ggtitle("Time point 2")+
  #theme(legend.position = "none")
  scale_color_manual(values=c("darkorange","pink","darkorange3","red"))


#DAPC: by treatment ---------------------------
library(adegenet)
dp=dapc(t(vsd),coldata$treatment,n.pca=4, n.da=1) 
# "predicting" the same data (for plotting)
pred=predict.dapc(dp,newdata=(t(vsd))) 

# measures of stress for each sample:
pred$ind.scores

# plotting "ghosts" 
plot(density(pred$ind.scores),col="white",bty="n",yaxt="n",ylab="",xlab="discriminant function",xlim=c(-4,4),ylim=c(0,1.8),main="",mgp=c(2.3,1,0))
polygon(density(pred$ind.scores[coldata$lab1=="cm1"]),border="black",col=alpha("pink",0.4))
polygon(density(pred$ind.scores[coldata$lab1=="hm1"]),border="black",col=alpha("red",0.4))
polygon(density(pred$ind.scores[coldata$lab1=="ca1"]),border="black",col=alpha("darkorange",.4))
polygon(density(pred$ind.scores[coldata$lab1=="ha1"]),border="black",col=alpha("darkorange3",.4))
# polygon(density(pred$ind.scores[coldata$lab1=="cm2"]),border="black",col=alpha("darkolivegreen2",.4))
# polygon(density(pred$ind.scores[coldata$lab1=="hm2"]),border="black",col=alpha("darkolivegreen",.4))
# polygon(density(pred$ind.scores[coldata$lab1=="ca2"]),border="black",col=alpha("cyan",.4))
# polygon(density(pred$ind.scores[coldata$lab1=="ha2"]),border="black",col=alpha("cyan4",.4))

legend(-3.9, 1.3, legend=c("cm1","hm1","ca1","ha1","cm2","hm2","ca2","ha2"),
       col=alpha(c("pink","red","darkorange","darkorange3","darkolivegreen2","darkolivegreen","cyan","cyan4"),0.6), pch = 19, cex=0.8,
       title="", text.font=4, bg=NA,box.lty = 0)
