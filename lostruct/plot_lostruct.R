library(lostruct)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(car)

#plot chromosome-level PCA
pca.full<-read.csv(file='~/Documents/StA/popgen/lostruct/lostruct_PCA.csv',h=T)
info<-read.table("~/Documents/StA/popgen/lostruct/allpops_outgroup_nodupes_pruned.fam",h=F,sep = "")
info<-data.frame(info[-c(307:310),])#remove commodus samples
info$region<-"Kauai"
info[grep("Hilo",info$V1),]$region<-"Hawaii"
info[grep("Oahu",info$V1),]$region<-"Oahu"
info[grep("Aus",info$V1),]$region<-"Australia"

#add sample info to PCA data frame
pca.full$pop<-rep(info$V1,14)
pca.full$region<-rep(info$region,14)
pca.full$sample<-rep(info$V2,14)
pca.full$fw<-rep(info$V6,14)

pca.plot<-ggplot(pca.full[!is.na(pca.full$chr),],aes(x=PC1,y=PC2))+
  theme_bw()+
  geom_point(size=0.75,alpha=0.75,aes(colour=region))+#facet_grid(.~chr,scales='free')+
  facet_grid(.~chr,scales='free')+
  theme(panel.grid=element_blank(),axis.text=element_blank(),
        legend.position='none',
        axis.ticks=element_blank(),
        strip.background=element_rect(fill='#eeeeee'))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  theme(strip.text.x = element_text(margin = margin(0,0,0,0, "cm")))

ggsave('chromosome_PCA.png',plot=pca.plot,
       dpi=600,height=1.4,width=10)


# now plot MDS in sliding windows across all chrs

### read existing files
mds.full<-read.csv("~/Documents/StA/popgen/lostruct/lostruct_MDS.csv")

colour_scheme<-c("#aaaaaa",'#5e82b5','#a060db','#c24a60')
mds.plot<-ggplot(mds.full,aes(x=bp/1e+06,y=value,colour=col))+
  geom_hline(yintercept=0,linewidth=0.5,colour='#888888',linetype='dashed')+
  facet_grid(variable~chr,scales='free')+
  theme_bw()+geom_point(size=0.1)+
  geom_smooth(method='loess',span=0.2,colour='black',linewidth=0.5,se=FALSE)+
  scale_colour_manual(values=colour_scheme)+
  theme(panel.grid=element_blank(),legend.position='none',
        axis.title=element_blank(),axis.text=element_blank(),
        axis.ticks=element_blank(),strip.background=element_rect(fill='#eeeeee'))+
  ggtitle("Evidence for massive inversions on at least 11 chromosomes")

library(patchwork)
mds.plot+pca.plot+plot_layout(nrow=2,heights=c(1.5,1))

ggsave('inversion_MDS_PCA.png',plot=mds.plot+pca.plot+plot_layout(nrow=2,heights=c(1.5,0.75)),
       dpi=600,height=5,width=10)


