#discriminant analysis of principal components (DAPC)
#fig 3b

library(DESeq2)
library(tidyverse)
library(ggplot2)
library(viridis)
library(vegan)
library(ggpmisc)
library(gridExtra)
rm(list=ls())

load("~/Desktop/tagseq_counts_coldata.Rdata")
coldata=subset(coldata,coldata$samplingTimepoint == "2")
coldata=subset(coldata,coldata$treatment == "control"|coldata$treatment=="heat")
coldata=subset(coldata,coldata$time_day == "morning")

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

datExpr=t(vsd)

#DAPC: by treatment ---------------------------
library(adegenet)
dp=dapc(t(vsd),coldata$treatment,n.pca=4, n.da=1) 
# "predicting" the same data (for plotting)
pred=predict.dapc(dp,newdata=(t(vsd))) 

# measures of stress for each sample:
pred$ind.scores

# plotting "ghosts" 
plot(density(pred$ind.scores),col="white",bty="n",yaxt="n",ylab="",xlab="discriminant function",ylim=c(0,1),xlim=c(-10,10),main="",mgp=c(2.3,1,0))
polygon(density(pred$ind.scores[coldata$treatment=="heat"]),col=alpha("brown2",0.4),border="black")
polygon(density(pred$ind.scores[coldata$treatment=="control"]),col=alpha("steelblue2",0.4),border="black")

#plot legend

legend(-3,1, legend=c("control","heat"),
       col=alpha(c("steelblue2","brown2"),0.4), pch = 19, cex=0.8,
       title="", text.font=1, bg=NA,box.lty = 0)
