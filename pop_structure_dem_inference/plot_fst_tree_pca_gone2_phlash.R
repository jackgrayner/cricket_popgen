#plot pairwise fst
library(pheatmap)
library(ggplot2)
library(dplyr)
library(ggplotify)
fst<-read.csv("./fst_results/fst.csv",row.names = 1)

#create and plot numeric matrix
fst<-apply(fst,1,as.numeric)
rownames(fst)<-colnames(fst)
#removal diagonal
diag(fst)<-NA
fst.heatmap<-pheatmap(fst,cluster_cols = FALSE,cluster_rows = FALSE,gaps_row = c(2,7),
         gaps_col = c(2,7),show_rownames = TRUE,show_colnames = TRUE,
         border_color = '#555555',display_numbers = FALSE)
#convert to ggplot object for multi-panel plot
fst.heatmap<-as.ggplot(fst.heatmap)

#read in phlash Ne estimates
phlash<-data.frame(t(read.csv("./phlash_results/Nes_Oahu_Chr5_9_all.csv",header=FALSE)))
#read in phlash time intervals
phlash.time<-read.csv("./phlash_results/Nes_Oahu_Chr5_9_T.csv",h=F)
phlash$T<-phlash.time$V1
q.0.25<-function(x){quantile(x,0.025)}
q.0.975<-function(x){quantile(x,0.975)}

#each column is a different iteration, so take the median for each sampled point (row)
phlash<-data.frame(apply(phlash,2,as.numeric))
phlash$median<-(apply(phlash[,c(2:101)],1,median))
phlash$q5<-apply(phlash[,c(2:101)],1,q.0.25)
phlash$q95<-apply(phlash[,c(2:101)],1,q.0.975)

phlash<-ggplot(phlash,aes(x=(T),y=(median)))+
  theme_bw()+theme(panel.grid=element_blank(),axis.ticks=element_line())+
  geom_ribbon(aes(ymin = (q5),ymax=(q95)),fill="#dddddd")+
  geom_line(size=1.5)+scale_y_log10()+scale_x_log10()+
  #scale_y_continuous(breaks=seq(4,9,1))+
  xlab("Generation")+annotation_logticks()+ylab("Effective population size")

#read in GONE2 results and plot
all<-rbind(
  read.table("./gone2_results/Hilo_Church_Lawn_all_u0.01_GONE2_Ne",h=T) %>% mutate(pop="Hawaii.CL")%>% mutate(island="Hawaii"),
  read.table("./gone2_results/Hilo_UH_all_u0.01_GONE2_Ne",h=T) %>% mutate(pop="Hawaii.UH") %>% mutate(island="Hawaii"),
  read.table("./gone2_results/Kauai_Wailua_all_u0.01_GONE2_Ne",h=T) %>% mutate(pop="Kauai.Wailua") %>% mutate(island="Kauai"),
  read.table("./gone2_results/Kauai_North_all_u0.01_GONE2_Ne",h=T) %>% mutate(pop="Kauai.North") %>% mutate(island="Kauai"),
  read.table("./gone2_results/Kauai_Pono_Kai_all_u0.01_GONE2_Ne",h=T) %>% mutate(pop="Kauai.Pono_Kai") %>% mutate(island="Kauai"),
  read.table("./gone2_results/Oahu_BYU_all_u0.01_GONE2_Ne",h=T) %>% mutate(pop="Oahu.BYU") %>% mutate(island="Oahu"),
  read.table("./gone2_results/Oahu_Community_Center_all_u0.01_GONE2_Ne",h=T) %>% mutate(pop="Oahu.CC") %>% mutate(island="Oahu"),
  read.table("./gone2_results/Oahu_Kamilo_all_u0.01_GONE2_Ne",h=T) %>% mutate(pop="Oahu.KP") %>% mutate(island="Oahu")
)

gone2<-ggplot(all,aes(x=Generation,y=(Ne_diploids),colour=pop))+
  theme_bw()+theme(panel.grid=element_blank(),axis.ticks=element_line(),legend.title=element_blank())+
  geom_rect(xmin=42,xmax=66.5,ymin=0,ymax=50000,fill='#eeeeee',colour='#eeeeee',alpha=0.5)+
  geom_vline(xintercept=115.5,colour='darkred',linetype='dotted')+
  geom_line(size=1)+scale_x_continuous(breaks=seq(0,150,20))+scale_y_log10()+
  #scale_y_continuous(breaks=seq(1,4,1))+
  annotation_logticks(sides = 'l')  +ylab("Effective population size")

#read in fasttree results and plot tree
library(ggtree)
library(ggplot2)
library(dplyr)
library(ape)

tree<-read.tree("./fasttree_results/fasttree_tree_file.txt")
pops<-read.table("./fasttree_results//pops.txt")

#change pop labels
pops$V2<-factor(pops$V2)
levels(pops$V2)<-c("Cairns","Hawaii.CL","Hawaii.UH","Kauai.AS","Kauai.CG","Kauai.PK",
                   "Kauai.VL","Kauai.WC","Mission","Oahu.AC","Oahu.BYU","Oahu.CC","Oahu.KP","Outgroup")
pops$island<-"commodus"
pops[grep("Kauai",pops$V2),]$island<-"Kauai"
pops[grep("Oahu",pops$V2),]$island<-"Oahu"
pops[grep("Hawaii",pops$V2),]$island<-"Hawaii"
pops[grep("Mission",pops$V2),]$island<-"Australia"
pops[grep("Cairns",pops$V2),]$island<-"Australia"

drop<-tree$tip.label[!tree$tip.label %in% pops$V1]
tree<-drop.tip(tree,drop)
pops[pops$island=="Australia" | pops$island=="commodus_outgroup",]$V2<-"outgroup/Australia"

#root tree by Australian mainland oceanicus sample
tree<-root(tree,outgroup = 'A4')
p <- ggtree(tree,size=0.5,aes(colour=V2)) %<+% pops +
  geom_tippoint(aes(color = V2,shape=island), size = 1) +  
  theme_tree2() +  
  theme(legend.title = element_blank())+
  scale_shape_manual(values=c('plus','plus','circle','square','triangle'))
scale_colour_manual(values=c("#d98558","black","#7caf5f","#a254a8","#3089c6")) #+

p<-flip(p, 377, 313)
print(p)

#perform and plot PCA

library(ggExtra)
library(patchwork)
library(pcadapt)
library(ggplotify)

info <- read.table("./pca_files/allHawaii_nodupes_pruned.fam",h=F)
pops<-info$V1
filename <- read.pcadapt("./pca_files/allpops_outgroup_nodupes_chr5_9_14_pruned.bed", type = "bed")
auto <- pcadapt(input = filename, K = 9,LD.clumping = list(size = 500, thr = 0.1))
plot(auto, option = "screeplot")#k=9 looks reasonable
pca.df<-data.frame(auto$scores)
pca.df<-pca.df[c(1:300),]
pca.df$population<-info$V1
pca.df$ID<-info$V2
pca.df$ID <- sub(".*-", "", pca.df$ID)
pca.df$island <- sub("_.*", "", pca.df$population)

#add fake aus point
pca.df[,1]<-as.numeric(pca.df[,1])
pca.df[,2]<-as.numeric(pca.df[,2])
pcaplot<-ggplot(pca.df,aes(x=X1,y=X2))+
  geom_point(aes(shape=island,colour=population),size=1.7)+
  scale_shape_manual(values=c('circle','square','triangle','plus'))+
  theme(panel.grid=element_blank())+xlim(c(-0.1,0.15))+ylim(c(-0.15,0.1))+
  theme_minimal()+theme(legend.position='none')+
  theme(panel.grid=element_blank(),axis.text=element_blank())+
  xlab("PC1")+ylab("PC2")


#plot all
library(patchwork)

ggsave("~/Documents/StA/popgen/writeup/pop_structure_demo.svg",height=10,width=14,
       plot=(pcaplot+labs(tag="A")+p+labs(tag="B")+fst.heatmap+labs(tag="C")+plot_layout(widths=c(0.66,0.66,1)))/
         (phlash+labs(tag="D")+gone2+labs(tag="E")))

