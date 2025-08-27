library(ggbeeswarm)
library(ggplot2)
library(viridis)
library(lme4)
library(car)
library(rptR)
library(dplyr)

flies<-read.csv('~/Documents/StA/popgen/flytrapping/flytrapping.csv',h=T)
flies<-flies[flies$rain %in% c("dry","light rain"),]#get rid of trials during heavy rain

sum(flies[flies$playback=="song",]$flies_caught)
sum(flies[flies$playback=="control",]$flies_caught)

#relabel sites
flies$site<-gsub("K\\.","Kauai.",flies$site)
flies$site<-gsub("O\\.","Oahu.",flies$site)
flies$site<-gsub("H\\.","Hawaii.",flies$site)

#classify sites according to proportion of singing crickets
flies$wtnwpres<-"no_song"
flies[flies$site %in% c("Kauai.CG","Kauai.HC","Oahu.CC","Oahu.KP"),]$wtnwpres<-"low_song"
flies[flies$site %in% c("Hawaii.CL","Kauai.PK","Hawaii.UH","Oahu.BYU"),]$wtnwpres<-"high_song"
flies$wtnwpres<-factor(flies$wtnwpres,levels=c("no_song","low_song","high_song"))

#create binary variable (0 = 0 flies, 1 = >0 flies)
flies$fly.yn<-0
flies[flies$flies_caught>0,]$fly.yn<-1

#now estimate relative fitness of silent males

#take mean flies attracted per song trap
flies.sum.night<-flies %>% filter(playback=='song' & time_postsunset<120) %>% 
  group_by(site,island,date,wtnwpres) %>% 
  dplyr::summarise(meanflies=mean(flies_caught),fly.yn=mean(fly.yn))

#that's for 10 mins (Pten). crickets sing for 7.63 mins (Pfly)
#'expected arrivals' per 7.63 minutes is
Nten = flies.sum.night$meanflies
Nfly = Nten^0.76

#probability of  being infected per night is...
#(only 61% of attacks end in infestation)
Pinf = 1 - (1 - 0.61)^(0.76*Nfly)

#PInf = 1 - exp(-0.61*Nfly)
est_lifespan = 1/Pinf + 6
est_lifespan[est_lifespan>29]<-29
flies.sum.night$lifeexpectency<-est_lifespan

NS_relfitness_nw<-(flies.sum.night$lifeexpectency)/29
SS_relfitness_nw<-1.391#sexual selection estimate from Tanner 2019
NS_relfitness_nw*SS_relfitness_nw # = net rel. fitness of singing Nw males
NS_relfitness_fw<-29/(flies.sum.night$lifeexpectency)
SS_relfitness_fw<-(1-0.391)#Sexual selection estimate taken from Tanner 2019
flies.sum.night$fw.rel.fitness<-NS_relfitness_fw*SS_relfitness_fw # = net rel. fitness of Fw males

g.silent<-ggplot(flies.sum.night,aes(x=site,y=NS_relfitness_fw,colour=island,group=date))+
  geom_hline(yintercept=1,linetype='dashed')+
  theme_bw()+theme(panel.grid=element_blank(),axis.text.x=element_text(angle=45,hjust=1,vjust=1),
                   axis.title.x=element_blank(),
                   strip.background=element_rect(fill='#eeeeee'))+
  facet_grid(.~wtnwpres,scales='free',space='free')+geom_quasirandom(alpha=0.75,width=0.1)+
  ylab("Rel. fitness of silent males under fly selection")+
  stat_summary(aes(group=site,fill=island),colour='#555555',shape=21,linewidth=0.5)+
  scale_colour_manual(values=c("#f8766d","#00c08b","#c77cff"))+
  scale_fill_manual(values=c("#f8766d","#00c08b","#c77cff"))#+scale_y_continuous(limits=seq(0,1))

ggsave("~/Documents/StA/popgen/flytrapping/rel_fitness.png",plot=g.silent1,dpi=600,height=3,width=5.25)

#plot num. of flies caught over time, with lines coloured by island
g.time<-ggplot(flies,aes(x=time_postsunset,y=flies_caught))+
  geom_rect(fill='#fcf3f2',xmin=0,xmax=120,ymin=0,ymax=1.5)+coord_cartesian(ylim=c(0,1.5))+
  theme_bw()+theme(panel.grid=element_blank(),legend.position='none')+
  geom_smooth(method='loess',span=1,aes(group=date,colour=island),se=FALSE,linewidth=0.35)+
  geom_smooth(method='loess',span=1,colour='black',alpha=1,fill='#cccccc',size=1.5)+
  ylab("Flies caught per trap")+xlab('Mins. post-sunset')+
  scale_colour_manual(values=c("#f8766d","#00c08b","#c77cff"))

#calculate repeatability of mean likelihood of attracting a fly (Y/N) across sites 
rep1<-rpt(formula=(log(fly.yn+1)~(1|site)),grname="site",data=data.frame(flies.sum.night),datatype = "Gaussian")
qqnorm(resid(rep1$mod))#resids look reasonably normally distributed
qqline(resid(rep1$mod))
shapiro.test(resid(rep1$mod))#as above
plot(rep1$mod)#could be heteroscedasticity
leveneTest(log(fly.yn+1) ~ factor(site), data = data.frame(flies.sum.night))#P = 0.10. Probably borderline 
plot(rep1)

#test for island/site diffs in likelihood of attracting fly
glm2<-lm(log(fly.yn+1)~island/site,data=data.frame(flies.sum.night))
anova(glm2)
plot(glm2)


#### investigate influence of rain and temperature
#load coordinates
sites<-read.csv("~/Documents/StA/popgen/map/cw_mapdata_silentvars.csv")

#read fly data
flies<-read.csv('~/Documents/StA/popgen/flytrapping/flytrapping.csv',h=T)
flies<-flies[flies$rain %in% c("dry","light rain"),]

flies.clim.clim<-merge(flies.clim.sum.night,sites,by.x='site',by.y='name1')

#retrieve daily temp
library(nasapower)
flies.clim$date1<-paste0(substr(flies.clim$date,7,11),"-",substr(flies.clim$date,4,5),"-",substr(flies.clim$date,1,2))
get_temp <- function(latitude, longitude, date) {
  result <- get_power(
    community = "AG",
    lonlat = c(longitude, latitude),
    pars = c("T2M"),  # T2M = Temperature at 2 meters
    dates = date,
    temporal_api = "daily"
  )
  print(result$T2M[1])
  return(result$T2M[1])
}

#retrieve daily rain
get_rain <- function(latitude, longitude, date) {
  result <- get_power(
    community = "AG",
    lonlat = c(longitude, latitude),
    pars = "PRECTOTCORR", 
    dates = date,
    temporal_api = "daily"
  )
  print(result$PRECTOTCORR[1])
  return(result$PRECTOTCORR[1])
}

#apply functions
flies.clim$temperature <- mapply(get_temp, flies.clim$latitude, flies.clim$longitude, flies.clim$date1)
flies.clim$rain <- mapply(get_rain, flies.clim$latitude, flies.clim$longitude, flies.clim$date1)

glm.temp<-glmer(log(meanflies.clim+1)~temperature+rain+(1|site),data=flies.clim)
glm.temp<-lm(log(meanflies.clim+1)~temperature+rain+Island/site,data=flies.clim)
qqnorm(resid(glm.temp))
qqline(resid(glm.temp))
car::Anova(glm.temp)
summary(glm.temp)

ggplot(flies.clim,aes(x=log(rain+1),y=log(meanflies.clim+1),colour=site))+
  theme_bw()+geom_point()+geom_smooth(method='lm',se=FALSE)+
  ylab("Daily mean flies.clim captured")+
  xlab("Log2 daily precipitation")#+facet_grid(.~temperature>median(temperature))




