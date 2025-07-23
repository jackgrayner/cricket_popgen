#RNA PCA

library(tidyr)
library(tidyverse)
library(ggbeeswarm)
library(plyr)
library(viridis)
library(patchwork)
library(ggrepel)
library(DESeq2)
library(edgeR)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(variancePartition)


setwd('~/Documents/StA/popgen/scripts/SVs_gene_expression/')

#first, read results of PCA from SNP calling
pca.all<-read.csv("./pca_allchr.csv")

pca.plot<-ggplot(pca.all,aes(x=PC1,y=PC2))+geom_point()+theme_bw()+
  theme(panel.grid.minor = element_blank())+
  facet_wrap(.~chr)+
  xlab(paste("PC1"))+
  ylab(paste("PC2"))

#define pca.allyotypes based on PCA clustering in the above plot
pca.all <- pca.all %>%
  mutate(Ninv = case_when(
    chr==1 & PC1 > 0 & PC1 < 0.12~1,
    chr==1 & PC1 > 0~0,
    chr==1 & PC1 < 0~2,
    
    chr==4 & PC1 < -0.2~2,
    chr==4 & PC1 < 0~1,
    chr==4 & PC1 > 0~0,
    
    chr==5 & PC1 < -0.2~2,
    chr==5 & PC1 < 0~1,
    chr==5 & PC1 > 0~0,
    
    chr==6 & PC1 > 0.2~2,
    chr==6 & PC1 > 0~1,
    chr==6 & PC1 < 0~0,
    
    chr==7 & PC1 > 0~1,
    chr==7 & PC1 < 0~0,
    
    chr==10 & PC1 < -0.1~2,
    chr==10 & PC1 < 0.1~1,
    chr==10 & PC1 > 0.1~0,
    
    chr==12 & PC1 > 0~2,
    chr==12 & PC1 > -0.1~1,
    chr==12 & PC1 < -0.1~0,
    
    chr==13 & PC1 > 0~2,
    chr==13 & PC1 > -0.2~1,
    chr==13 & PC1 < -0.2~0,
    
    TRUE~NA_real_
  ))

pca.plot<-ggplot(pca.all,aes(x=PC1,y=PC2,fill=Ninv))+geom_point(size=2,shape=21,colour='#555555')+theme_bw()+
  theme(panel.grid.minor = element_blank(),legend.position='none',axis.ticks=element_blank(),
        axis.text=element_blank(),panel.grid=element_blank(),
        panel.background = element_rect(fill='white',colour='#aaaaaa',linewidth=0.25),
        plot.title=element_text(hjust = 0.5),
        strip.background = element_rect(fill="#eeeeee"),
        strip.text = element_text(size = 10, margin = margin()))+
  facet_wrap(~chr,ncol=5)+scale_fill_viridis(discrete=FALSE)+
  xlab(paste("PC1"))+
  ylab(paste("PC2"))

#now run DE analysis

pheno<-read.csv('sample_info.csv',row.names=1)
cts<-read.csv('Toc_e_gene_count_matrix.csv',h=T,row.names="gene_id")
cts<-cts[-nrow(cts),]

#remove genes not expressed at >1 count per mill. in at least 3 samples
cts<-cts[rowSums(cpm(cts)>1)>=3, ]
nrow(cts)

#run variance partitioning of gene expression counts into predictors and SVs
form <- ~  (1 |Fw_geno) + (1 | Treatment) + (1 | Sex) + (Chr1_SV) + (Chr4_SV) + (Chr5_SV) + (Chr6_SV) + (Chr7_SV) +
  (Chr10_SV) + (Chr12_SV) + (Chr13_SV) 
varPart <- fitExtractVarPartModel(cts, form, pheno,n=5000)
vp <- sortCols(varPart)
vp<-plotVarPart(vp)$data

vp$variable<-factor(vp$variable,levels=c("Chr1_SV","Chr4_SV","Chr5_SV","Chr6_SV","Chr7_SV","Chr10_SV",
                                         "Chr12_SV","Chr13_SV","Fw_geno","Sex","Treatment","Residuals"))
vp$var1<-vp$variable
vp[vp$variable %in% c("Fw_geno","Sex","Treatment"),]$var1<-NA
g2<-ggplot(vp[!vp$variable=="Residuals" & !is.na(vp$variable),],aes(x=variable,y=value,fill=var1))+
  theme_bw()+theme(panel.grid=element_blank())+
  theme(axis.text.x=element_text(angle = 45,vjust=1, hjust=1,size=8),axis.title.y=element_blank(),
        panel.background = element_rect(fill='white',colour='#aaaaaa',linewidth=0.25))+
  geom_quasirandom(size=0.5,alpha=0.5,colour='#888888')+
  #geom_violin()+
  geom_boxplot(colour='black',width=0.5,alpha=0.75,outlier.shape=NA)+
  theme(legend.position='none')+ylab("Variance explained (%)")
g2<-g2+coord_flip(ylim=c(0,45)) 


#DE analysis
results_list<-list()
for (i in seq_along(c(1,4,5,6,7,10,12,13))) {
  chr<-c(1, 4, 5, 6, 7, 10, 12, 13)[i]
  pheno1<-pheno[, c(1, 2, 3, i + 3)]
  colnames(pheno1)[4]<-'chr_SV'
  
  dds <- DESeq(DESeqDataSetFromMatrix(
    countData=cts,
    colData=pheno1,
    design=~Sex+chr_SV+Treatment+Fw_geno))
  res<-results(dds, name = "chr_SV", alpha = 0.05)
  results_list[[paste0("chr",chr)]]<-data.frame(res,chr.sv=chr)
}
res.all<-do.call(rbind,results_list)
res.all$gene<-gsub("chr[0-9]*\\.","",rownames(res.all))
res.all<-res.all[!is.na(res.all$padj),]

#get gene positions and annotations
tocgff<-read.table("~/Documents/StA/Sw_WGS/TOC.asm.scaffold.gene.gff3")
tocgff$V1<-as.integer(gsub("scaffold_","",tocgff$V1))
tocgff<-tocgff[tocgff$V3=="gene",]
tocgff$gene<-gsub(";","",substr(tocgff$V9,4,16))

#run this to get annotations - but exclude from main analysis
#gene_annos<-read.table("~/Documents/StA/Sw_WGS/TOC.asm.scaffold.gene.SWISSPROT.blastp.top_hit.txt")
#colnames(gene_annos)<-c("gene","score","anno")
#gene_annos$gene<-gsub(".t[0-9]*","",gene_annos$gene)
#tocgff<-gtf[tocgff$gene %in% res.all$gene,]
#tocgff<-merge(tocgff,gene_annos,by='gene')

res.all<-merge(res.all,tocgff,by='gene')
res.all<-res.all[res.all$V1<15,]#get rid of genes on unplaced scaffolds

library(viridis)
dges<-ggplot(res.all[!is.na(res.all$pvalue) & !is.na(res.all$chr.sv) &
                           !is.na(res.all$V1),],
       aes(x=V4,y=(-log10(pvalue)),colour=factor(chr.sv)))+
  theme_minimal()+theme(axis.text.x=element_blank(),panel.grid=element_blank(),axis.title.x=element_blank(),
                   panel.background = element_rect(fill='white',colour='white',linewidth=0.25),
                   axis.ticks=element_blank(),
                   legend.position='top',strip.background = element_rect(fill='white',colour='white'),
                   legend.title=element_blank())+
                     guides(colour = guide_legend(override.aes = list(size=2,alpha=1)))+
  geom_point(size=0.5,alpha=1)+facet_grid(.~V1,scales = 'free',switch='both')+
  ylab("-log10(P)")#+
  scale_color_viridis(discrete=TRUE,option="E")

ggsave('DE_effects_inversions_genome.png',plot=dges,dpi=600,height=4,width=9)


