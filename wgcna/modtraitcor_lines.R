#make plots for reviewer 1
rm(list=ls())
library(ggplot2)
library(tidyverse)

#line plot, treatment change for each module
load("~/Dropbox/project/Projects/stsw/mdRAD/ixn_term_investigation/tagseq_wgcna/wgcna_input.Rdata")
load("~/Dropbox/project/Projects/stsw/mdRAD/ixn_term_investigation/tagseq_wgcna/networkdata_signed.RData")
load("~/Dropbox/project/Projects/stsw/submission_materials/additional_analyses/moduleTraitCor.Rdata")

rownames(moduleTraitCor)=c("brown","greenyellow","black","yellow","blue","grey60","salmon","magenta","pink")

moduleTraitCor1 = as.data.frame(moduleTraitCor)
moduleTraitCor1$module = rownames(moduleTraitCor1)
moduleTraitCor1= moduleTraitCor1 %>% pivot_longer(cols=c("c1","h1","s1","c2","h2","s2"),
                    names_to='treat',
                    values_to='R')
moduleTraitCor1=separate(moduleTraitCor1, treat, into = c("treatment", "timepoint"), sep = 1)

#assign colors
moduleTraitCor1$colors= ifelse(moduleTraitCor1$treatment == "c", "steelblue",
                               ifelse(moduleTraitCor1$treatment =="h","coral2","orchid4"))

#make nicer legend labels
moduleTraitCor1$treatment= ifelse(moduleTraitCor1$treatment == "c", "control",
                               ifelse(moduleTraitCor1$treatment =="h","heat","switch"))

                              
ggplot(data = moduleTraitCor1, aes(timepoint, R,group=treatment)) +
  theme_bw()+
  geom_line(aes(color=treatment)) +
  geom_point(aes(color=treatment),size=2) + 
  ylab("expression level")+
  #theme(legend.position="none")+
  geom_hline(yintercept = 0,color="gray",linetype="dashed")+
  scale_color_manual(values=c("steelblue","coral2","orchid4"))+
  facet_wrap(~ module)

