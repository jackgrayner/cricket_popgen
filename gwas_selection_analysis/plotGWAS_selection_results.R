library(dplyr)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(topGO)
library(GenomicRanges)

# Subsetting - run only first time to avoid reading in large files, and add useful columns

# Oahu<-read.table("Oahu_Fw.assoc.txt",h=T) %>% mutate(Island="Oahu",Padj=p.adjust(p_lrt))
# Kauai<-read.table("Kauai_Fw.assoc.txt",h=T) %>% mutate(Island="Kauai",Padj=p.adjust(p_lrt))
# Hilo<-read.table("Hilo_Fw.assoc.txt",h=T) %>% mutate(Island="Hawaii",Padj=p.adjust(p_lrt))
# write.table(Oahu[Oahu$p_lrt<0.1,],"Oahu_Fw.assoc_P0.1.txt",row.names=FALSE,quote=FALSE)
# write.table(Kauai[Kauai$p_lrt<0.1,],"Kauai_Fw.assoc_P0.1.txt",row.names=FALSE,quote=FALSE)
# write.table(Hilo[Hilo$p_lrt<0.1,],"Hawaii_Fw.assoc_P0.1.txt",row.names=FALSE,quote=FALSE)

#######
#read/filter genome annotation, add functional annotations
######

#read genome annotation, reformat some columns, subset genes
tocgff<-read.table("TOC.asm.scaffold.gene.gff3") %>% 
  mutate(V1=as.integer(gsub("scaffold_","",V1)),gene=gsub("ID=.*;.*Name=","",V9)) %>% filter(V3=="gene")
#read functional annotations
gene_annos<-read.table("TOC.asm.scaffold.gene.SWISSPROT.blastp.top_hit.txt")
colnames(gene_annos)<-c("gene","score","anno")
gene_annos$gene<-gsub("\\.t.*","",gene_annos$gene)#remove transcript ID

#keep top hit per gene
gene_annos<-gene_annos[order(gene_annos$score,decreasing = TRUE),] %>% filter(!duplicated(gene_annos))

#remove accession info and spp ID from annotations for readability, then merge genome and functional annotation data frames
gene_annos$anno<-gsub(".*\\|","",gene_annos$anno)
gene_annos$anno<-gsub("_.*","",gene_annos$anno)#remove spp ID from annotation
tocgff.anno<-merge(tocgff,gene_annos,by='gene')

#create granges object
genes<-GRanges(seqnames=tocgff.anno$V1,#chr/scaffolds
  ranges=IRanges(start=tocgff.anno$V4, end=tocgff.anno$V5))#start and end of gene feature
mcols(genes)$gene <- tocgff.anno$gene#add gene ID (e.g., TOC.11T.00352)
mcols(genes)$anno <- tocgff.anno$anno#add annotation (eg., UBCD1)

window=50000 #50kb window used to identify nearby genes
#find nearby genes for sig. variants, output sig. vars with genes & annotations
annotation_function<-function(df){ 
  #make empty columns
  df$gene<-NA#to store gene IDs for nearby genes
  df$anno<-NA#to store annotations for nearby genes
  df$dist<-NA#to store distance to the nearby gene (0 if overlapping)
  df$within<-FALSE#whether the variant is within a gene feature
  
  #convert input into granges object
  variants<-GRanges(seqnames=df$chr,ranges=IRanges(start=df$ps,end=df$ps))
  #find overlaps between query/subject (i.e. variants within genes)
  variants_within<-findOverlaps(variants, genes)
  
  #store within (TRUE) and distance (0) for overlapping ranges
  df$within[queryHits(variants_within)]<-TRUE
  df$dist[queryHits(variants_within)]<-0
  
  #find nearby genes
  dist.nearest<-distanceToNearest(variants,genes,ignore.strand=TRUE)
  nearest_genes<-genes[subjectHits(dist.nearest)]
  distances<-mcols(dist.nearest)$distance
  #save results to variants df
  df$gene<-mcols(nearest_genes)$gene
  df$anno<-mcols(nearest_genes)$anno
  df$dist<-distances
  #note whether dist is < 50Kb
  df$within50kb<-df$dist<window
  return(df)
}

#custom theme for plots
cust.theme<-function(){
  theme_minimal() %+replace% 
    theme(panel.grid=element_blank(),plot.background=element_rect(fill='white',colour='white'),
          axis.ticks.y=element_blank(),legend.position='none',
          panel.spacing.x=unit(0.05, "lines"))
}


#######
#curly-wing GWAS - plot (only Chr2)
######

#read subset of association test results, keep only Chr2
All.Cw<-read.table("All_Cw_lmm.assoc_P0.1.txt",h=T) %>% filter(chr==2)
#annotate sig. SNPs
All.Cw$sig<-All.Cw$padj<0.05
All.Cw<-annotation_function(All.Cw)

All.Cw.plot<-ggplot(All.Cw,aes(x=ps,y=(-log10(p_lrt)),colour=sig))+
  cust.theme()+
  geom_point(size=0.5,alpha=1)+
  #plot gene annotations for SNPs that are significant and within 50kb of genes
  geom_text_repel(data=All.Cw[All.Cw$sig & All.Cw$within50kb,],
                  aes(x=ps,y=(-log10(p_lrt)),label=anno),
                  min.segment.length = 0.0001,size=2.5,
                  nudge_x = 0,nudge_y = 1.5,fontface = "italic",alpha=1)+
  scale_colour_manual(values=c("#d4cdcd","#8f3131"))+
  ggtitle("Curly-wing")+ylab("-log10(P)")+
  scale_y_continuous(expand=c(0,0))

rm(list="All.Cw")#remove from environment to save memory

#######
#flatwing GWAS - plot (only Chr1)
######

#read subset of association test results for each pop, keep only Chr1
kauai.fw<-read.table("Kauai_Fw.assoc_P0.1.txt",h=T) %>% filter(chr=="scaffold_1") %>% mutate(chr=1)
oahu.fw<-read.table("Oahu_Fw.assoc_P0.1.txt",h=T) %>% filter(chr=="scaffold_1") %>% mutate(chr=1)
#hilo.fw<-read.table("Hawaii_Fw.assoc_P0.1.txt",h=T) %>% filter(chr=="scaffold_1") #excl. as too few Fw samples

#store positions of top variants for plotting
topfw.kauai<-kauai.fw[order(kauai.fw$p_lrt),][1,]$ps
topfw.oahu<-oahu.fw[order(oahu.fw$p_lrt),][1,]$ps

#add annotations
kauai.fw<-annotation_function(kauai.fw)
oahu.fw<-annotation_function(oahu.fw)

#concatenate
All.Fw<-rbind(kauai.fw,oahu.fw)
rm(list=c("kauai.fw","oahu.fw"))#save memory

#define whether SNP is significant in either island
All.Fw$sig<-NA
All.Fw[All.Fw$Padj<0.05 & All.Fw$Island=="Kauai",]$sig<-"Kauai"
All.Fw[All.Fw$Padj<0.05 & All.Fw$Island=="Oahu",]$sig<-"Oahu"

All.Fw.plot<-ggplot(All.Fw,aes(x=ps/1e+6,y=(-log10(p_lrt)),colour=sig))+facet_grid(Island~.)+
  xlab("Chr1 pos. Mb")+cust.theme()+
  geom_vline(xintercept=253.6,linetype='dashed',linewidth=0.25,colour='black')+#add dsx location
  geom_point(size=0.5,alpha=1)+
  #plot gene annotations for SNPs that are significant and within 50kb of genes
  geom_text_repel(data=All.Fw[!is.na(All.Fw$sig)  & All.Fw$within50kb,],
                  aes(x=ps/1e+6,y=(-log10(p_lrt)),label=anno),colour='black',max.overlaps = 10,
                  min.segment.length = 0.0001,size=2,nudge_x = 0,nudge_y = 1.5,fontface = "italic",alpha=1)+
  scale_colour_manual(values=c("#00c08b","#c77cff","#d4cdcd"))+
  ggtitle("Flatwing")+ylab("-log10(P)")+
  scale_y_continuous(expand=c(0,0))#+coord_cartesian(xlim=(c(250,260)))

rm(list="All.Fw")

#######
#selection results (LFMM of genotype ~ fly attack rate)
######

flysel<-read.csv("LFMM_results.csv",h=T) %>% mutate(p_lrt=P)
flysel$sig<-flysel$p_lrt<quantile(flysel$p_lrt,0.0001)

#add gene info 
flysel<-annotation_function(flysel)

fly.sel.plot<-ggplot(flysel,aes(x=ps,y=(-log10(P)),colour=sig))+
  geom_vline(data=data.frame(chr=1),aes(xintercept=topfw.kauai),linetype='dashed',linewidth=0.35,colour='#00c08b')+
  geom_vline(data=data.frame(chr=1),aes(xintercept=topfw.oahu),linetype='dashed',linewidth=0.35,colour='#c77cff')+
  #plot gene annotations for SNPs that are outliers and within 50kb of genes
  geom_text_repel(data=flysel[flysel$sig & flysel$within50kb,],aes(x=ps,y=(-log10(P)),label=anno),
                  min.segment.length = 0.0001,size=2,nudge_x = 0,nudge_y = 1.5,
                  fontface = "italic",alpha=1,max.overlaps = 10)+
  cust.theme()+theme(axis.text.x=element_blank())+
  facet_grid(.~chr,space='free',scales='free',switch='both')+
  geom_point(size=0.5,alpha=1)+
  scale_colour_manual(values=c("#d4cdcd","#8f3131"))+
  ggtitle("Fly selection")+ylab("-log10(P)")+
  scale_y_continuous(expand=c(0,0),limits=c(0,7.3))

ggsave('fw_cw_fly.png',dpi=600,height=6,width=8.5,
       plot=(All.Fw.plot.zoomed+labs(tag="C")+All.Cw.plot.zoomed+labs(tag="D"))/
         fly.sel.plot+labs(tag="E")+plot_layout(nrow=3))


#now test overrepresentation of gene ontology categories among selection outliers

#read in gene ontology file
gaf<-read.table("~/Documents/StA/RNA/makeshift_gaf.tsv",h=F,sep='\t')
colnames(gaf)<-c("gene_ID","go_terms")

#test overrepresentation of gene ontology categories among selected loci
#using topGO functions

#create named list of P-values (names=gene IDs)
gene_list <- flysel$p_lrt
names(gene_list) <- flysel$gene
gene_list <- gene_list[!is.na(gene_list)]

#define function for genes of interest (0.01 percentile P outliers)
selection_function <- function(x) {
  return(x < quantile(flysel$P,0.0001))
}

#create list
gene2GO <- list()
for(i in 1:nrow(gaf)){
  #for each gene with GO annotations, split them by ; character, then create list (GOs per gene)
  if(!is.na(gaf[i,]$go_terms) & !gaf[i,]$go_terms=="") {
    go_list <- (unlist(strsplit(gaf[i,]$go_terms,";")))
    gene2GO[[gaf[i,]$gene_ID]] <- go_list
  }
}

#create topGO gene ontology objects
GOdata_BP <-  new("topGOdata",
                  ontology = "BP",#biological process
                  allGenes = gene_list,#full list of genes with GO terms
                  geneSel = selection_function,#define outliers
                  annot = annFUN.gene2GO,
                  gene2GO = gene2GO)
GOdata_MF <- new("topGOdata",
                 ontology = "MF",#molecular function
                 allGenes = gene_list,#full list of genes with GO terms
                 geneSel = selection_function,#define outliers
                 annot = annFUN.gene2GO,
                 gene2GO = gene2GO)

#run Fisher's exact tests
test_BP <- runTest(GOdata_BP, statistic = "fisher")
test_MF <- runTest(GOdata_MF, statistic = "fisher")
#generate results tables
results_BP <- GenTable(GOdata_BP, Fisher = test_BP, topNodes = 100) %>% mutate(Category="Biological Process")
results_MF <- GenTable(GOdata_MF, Fisher = test_MF, topNodes = 100) %>% mutate(Category="Molecular Function")

#concatenate BP and MF results, then subset terms with P < 0.01 and more than 1 outlier
all_results <- rbind(results_BP, results_MF)
all_results$Fisher<-as.numeric(all_results$Fisher)
all_results<-all_results[all_results$Fisher<0.01 & all_results$Significant>1,]
View(all_results)
write.csv(all_results,"fly_selection_outlier_GOoverrep.csv",row.names=FALSE)

