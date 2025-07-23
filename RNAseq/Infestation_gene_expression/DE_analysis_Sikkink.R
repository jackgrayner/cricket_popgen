library(DESeq2)
library(topGO)
library(edgeR)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(stringr)

setwd('~/Documents/StA/popgen/scripts/Infestation_gene_expression/')

gaf<-read.table("~/Documents/StA/RNA/makeshift_gaf.tsv",h=F,sep='\t')
#create gene ontology overrep. functions
run_topGO<-function(res.df){
  gene_list<-res.df$padj
  names(gene_list)<-rownames(res.df)
  gene_list<-gene_list[!is.na(gene_list)]
  selection_function<-function(x) {
    return(x<0.05)#test geneset defined as Padj < 0.05; no logFC threshold
  }
  gene2GO<-list()
  for(i in 1:nrow(gaf)) {
    gene_symbol<-gaf[i,1] 
    go_terms<-gaf[i,2] 
    
    if(!is.na(go_terms) & !go_terms=="") {
      go_list<-trimws(unlist(strsplit(go_terms,";")))
      go_list<-go_list[go_list!=""]
      gene2GO[[gene_symbol]]<-go_list
    }
  }
  #test biological processes
  GOdata_BP <-  new("topGOdata",
                    ontology = "BP",
                    allGenes = gene_list,
                    geneSel = selection_function,
                    annot = annFUN.gene2GO,
                    gene2GO = gene2GO)
  #test molecular functions
  GOdata_MF <- new("topGOdata",
                   ontology = "MF",
                   allGenes = gene_list,
                   geneSel = selection_function,
                   annot = annFUN.gene2GO,
                   gene2GO = gene2GO)
  
  #run fisher tests
  test_BP <- runTest(GOdata_BP, algorithm = "weight01", statistic = "fisher")
  test_MF <- runTest(GOdata_MF, algorithm = "weight01", statistic = "fisher")
  results_BP <- GenTable(GOdata_BP, Fisher = test_BP, topNodes = 100)
  results_MF <- GenTable(GOdata_MF, Fisher = test_MF, topNodes = 100)
  results_BP$Category <- "Biological Process"
  results_MF$Category <- "Molecular Function"
  
  all_results <- rbind(results_BP, results_MF)
  all_results$Fisher<-as.numeric(all_results$Fisher)
  all_results<-all_results[all_results$Fisher<0.01,]#output only GO terms with P < 0.01
  return(all_results)
}

setwd('~/Documents/StA/RNA')
pheno<-read.csv('Sikkink_pheno.csv')
rownames(pheno)<-pheno$ID
cts<-read.csv('Sikkink_gene_count_matrix.csv',h=T,row.names="gene_id")
summary(colnames(cts) %in% rownames(pheno))
cts<-cts[,order(colnames(cts))]
pheno<-pheno[order(rownames(pheno)),]
summary(colnames(cts) == rownames(pheno))#check samples in same order

#DE analysis
dds <- DESeq(DESeqDataSetFromMatrix(countData = cts,
  colData = pheno,
  design= ~ Pop+Trt))#nb. two populations in analysis
resultsNames(dds)
summary(results(dds, name="Trt_D4_vs_C",alpha=0.05))#d4 vs c is the focus as it represents mid-stage infection
summary(results(dds, name="Trt_D7_vs_C",alpha=0.01))
res.d4<-data.frame(results(dds, name="Trt_D4_vs_C",alpha=0.05))
res.d4$gene<-rownames(res.d4)
res.d7<-data.frame(results(dds, name="Trt_D7_vs_C",alpha=0.05))
res.d7$gene<-rownames(res.d7)

#GO overrep
res.d4.go<-run_topGO(res.d4)
res.d4.go$fold<-res.d4.go$Significant/res.d4.go$Annotated
View(res.d4.go)


#plot volcano and GO overrep
g.volcano<-ggplot(res.d4,aes(x=log2FoldChange,y=-log10(pvalue)))+
  theme_bw()+theme(panel.grid=element_blank(),legend.position='none')+
  geom_point(size=0.75,aes(colour=padj<0.05))+
  geom_vline(xintercept=0,linetype='dashed')+
  xlim(c(-12,12))+
  scale_colour_manual(values=c("#aaaaaa","#C74955"))
g.go<-ggplot(res.d4.go[res.d4.go$Fisher<3.5e-5,],aes(y=Term,x=-log10(Fisher)))+
  theme_minimal()+theme(panel.grid=element_blank(),legend.position='none',
                        axis.title.y=element_blank(),axis.line.x=element_line())+
  geom_col(colour='black')+xlab("Fold enrichment")

g.volcano+g.go
ggsave('Infest_d4_vol_go.png',plot=g.volcano+labs(tag="F",title="Response to infestation")+
         g.go+theme(axis.text.y=element_text(size=8))+labs(tag="G")+
         plot_layout(widths=c(1.25,1)),height=2.6,width=7.5)
