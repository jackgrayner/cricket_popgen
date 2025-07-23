library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(viridis)
library(ape)
library(ggrepel)

setwd("~/Documents/StA/popgen/map/")

#create barplot of morph frequencies
morph<-read.table("~/Documents/StA/popgen/ordered_sample_info.txt",h=T)
morph$Cw_pheno<-factor(morph$Cw_pheno)
morph$Fw_pheno<-factor(morph$Fw_pheno)
morph$Sw_pheno<-factor(morph$Sw_pheno)
levels(morph$Cw_pheno)<-c("Straight","Curly")
levels(morph$Fw_pheno)<-c("Normal","Flat")
levels(morph$Sw_pheno)<-c("Long","Small")
summary(factor(morph$Population))

morph$Phenotype<-factor(paste(morph$Sw_pheno,morph$Cw_pheno,morph$Fw_pheno,sep="_"),
                    levels=c("Long_Curly_Flat","Long_Curly_Normal","Long_Straight_Flat",
                             "Small_Curly_Normal","Small_NA_Normal","Long_Straight_Normal"))

#plot morph frequencies
g.freqs<-ggplot(morph,aes(x=Population,fill=Phenotype))+theme_minimal()+scale_fill_viridis(discrete=TRUE,option='E')+
  geom_bar(position='fill',alpha=1)+theme(panel.grid=element_blank(),axis.text.x=element_text(angle=90,hjust = 1,vjust = 0.5))+
  theme(legend.position='right',legend.title=element_blank(),axis.title=element_blank(),
        legend.text = element_text(size = 8),legend.key.size = unit(0.45, "cm"))+
  guides(shape = guide_legend(override.aes = list(size = 0.5)))+guides(color = guide_legend(override.aes = list(size = 0.5)))

ggsave('~/Documents/StA/popgen/flytrapping/morph_bar.png',plot=g.freqs,dpi=600,height=2.25,width=3.5)

#make map
world <- ne_countries(scale = "medium", returnclass = "sf")

#import data 
locations <- read.csv("cw_mapdata_silentvars.csv")
locations$fontface<-'bold'
locations[locations$singing_remain=="N",]$fontface<-'italic'
Map_Populations_crickets<-ggplot(data = world) + 
  geom_sf(fill= 'antiquewhite',colour="black") +
  geom_point(data=locations,aes(x=longitude, y=latitude, fill=island),shape=21,size=1,alpha=1,colour='black')+
  geom_text_repel(data=locations,
                  aes(x=longitude, y=latitude, label=name1, colour=island), 
                  size=4,min.segment.length = unit(0, 'lines'),force = 50,
                  show.legend = FALSE,force_pull = 0.1) + 
  annotation_scale(location = 'bl', width_hint = 0.5) +
  scale_colour_manual(values=c("#f8766d","#00c08b","#c77cff"))+
  scale_fill_manual(values=c("#f8766d","#00c08b","#c77cff"))+
  annotation_north_arrow(location = 'bl', which_north = 'true', pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering) + 
  coord_sf(xlim = c(-160.580614, -154.278239), ylim = c(18.268486, 22.798264), expand = FALSE)+
  xlab('Longitude') + 
  ylab('Latitude') + 
  theme(panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "white",colour='white'))+
  theme(legend.position = 'none')+
  guides(colour = guide_legend(override.aes = list(size=3)))

ggsave('patch_map_morph.png',plot=grid.arrange(Map_Populations_crickets+labs(tag="A"),g.freqs+labs(tag="B"),
                                               widths=c(1.5,1)),height=6,width=8,dpi=600)

