library(lme4)
library(lostruct)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(car)
setwd("~/Documents/StA/popgen/ld_lostruct_pca/")
{#once everything has been run, this can be closed
{#' Minimum Enclosing Circle
#'
#' Finds the minimum enclosing circle to a set of points,
#' "Skyum algorithm based on the convex hull", as implemented at 
#'    http://www.uni-kiel.de/psychologie/rexrepos/posts/diagBounding.html#minimum-enclosing-circle
#'
#' @param xy The 2D coordinates.
#' @param circ A circle as returned by \code{enclosing_circle}.
#' @param n The number of angles to draw the circle at.
#' @param plot.it Actually do plotting or just return points?
#' @return A named list with entries \code{ctr} (the coordinates of the center), \code{rad} (the radius),
#' \code{three} (the three points in xy that define the circle),
#' and \code{index} (the indices of those three points in xy).
enclosing_circle <- function (xy) { 
  # deal with NAs, which chull doesn't like.
  na.inds <- is.na(xy[,1]) | is.na(xy[,2])
  out <- getMinCircle(xy[!na.inds,]) 
  out$index <- match( out$index, cumsum(!na.inds) )
  return(out)
}

#' @describeIn enclosing_circle Plots a circle.
plot_circle <- function (circ,n=200, plot.it=TRUE, ...) {
  angles <- seq(0, 2*pi, length.out=n)
  circ.lines <- cbind(circ$ctr[1] + circ$rad*cos(angles),
                      circ$ctr[2] + circ$rad*sin(angles))
  if (plot.it) lines(circ.lines, ... )
  return(invisible(circ.lines))
}

#' Circle defined by three points
getCircleFrom3 <- function(xy) {
  stopifnot(is.matrix(xy), is.numeric(xy), nrow(xy) == 3, ncol(xy) == 2)
  
  aa <- xy[1,  ]
  bb <- xy[2,  ]
  cc <- xy[3,  ]
  y  <- xy[ , 2]
  
  xDeltaA <- bb[1] - aa[1]
  yDeltaA <- bb[2] - aa[2]
  xDeltaB <- cc[1] - bb[1]
  yDeltaB <- cc[2] - bb[2]
  xDeltaC <- cc[1] - aa[1]
  yDeltaC <- cc[2] - aa[2]
  
  ## check if the points are collinear: qr(xy)$rank == 1, or:
  ## determinant of difference matrix = 0, no need to use det()
  dMat <- rbind(c(xDeltaA, yDeltaA), c(xDeltaB, yDeltaB))
  if(isTRUE(all.equal(dMat[1,1]*dMat[2,2] - dMat[1,2]*dMat[2,1], 0, check.attributes=FALSE))) {
    ## define the circle as the one that's centered between the points
    rangeX <- range(c(aa[1], bb[1], cc[1]))
    rangeY <- range(c(aa[2], bb[2], cc[2]))
    ctr    <- c(rangeX[1] + 0.5*diff(rangeX), rangeY[1] + 0.5*diff(rangeY))
    rad    <- sqrt((0.5*diff(rangeX))^2 + (0.5*diff(rangeY))^2)
  } else {
    rad <- prod(dist(xy)) / (2 * abs(det(cbind(xy, 1))))  # circle radius
    v1  <- rowSums(xy^2)                    # first vector in the numerator
    v2x <- c( xDeltaB, -xDeltaC,  xDeltaA)  # 2nd vector numerator for Mx
    v2y <- c(-yDeltaB,  yDeltaC, -yDeltaA)  # 2nd vector numerator for My
    ctr <- c(t(v1) %*% v2y, t(v1) %*% v2x) / as.vector(2 * (t(y) %*% v2x))  # center
  }
  
  return(list(ctr=ctr, rad=rad))
}

#' Vertex that produces the circle with the maximum radius
getMaxRad <- function(xy, S) {
  stopifnot(is.matrix(xy), is.numeric(xy), nrow(xy) >= 2, ncol(xy) == 2)
  stopifnot(is.numeric(S), length(S) >= 2, length(S) <= nrow(xy))
  
  n    <- length(S)                    # number of points
  Sidx <- seq(along=S)                 # index for points
  rads <- numeric(n)                   # radii for all circles
  post <- (Sidx %% n) + 1              # next point in S
  prev <- Sidx[order(post)]            # previous point in S
  for(i in Sidx) {
    pts     <- rbind(xy[S[prev[i]], ], xy[S[i], ], xy[S[post[i]], ])
    rads[i] <- getCircleFrom3(pts)$rad  # circle radius
  }
  
  return(which.max(rads))
}

#' Check if the angle at B in triangle ABC exceeds 90 degrees
isBiggerThan90 <- function(xy) {
  stopifnot(is.matrix(xy), is.numeric(xy), nrow(xy) == 3, ncol(xy) == 2)
  d   <- dist(xy)
  dAB <- d[1]
  dAC <- d[2]
  dBC <- d[3]
  return((dAB^2 + dBC^2 - dAC^2) < 0)
}

#' Maximum pairwise distance between two 2D-points
getMaxPairDist <- function(xy) {
  stopifnot(is.matrix(xy), is.numeric(xy), ncol(xy) == 2, nrow(xy) >= 2)
  
  # 2D -> only convex hull is relevant
  H    <- chull(xy)      # convex hull indices (vertices ordered clockwise)
  pts  <- xy[H, ]        # points that make up the convex hull
  N    <- nrow(pts)                      # number of points on hull
  dMat <- dist(pts, method="euclidean")  # distance matrix
  idx  <- which.max(as.matrix(dMat))     # maximum distance
  i    <- (idx-1) %/% N+1                # column -> point 1
  j    <- (idx-1) %%  N+1                # row    -> point 2
  mPts <- H[c(i, j)]                     # rows with max distance
  dst  <- max(dMat)                      # max distance
  
  return(list(d=dst, idx=mPts))
}

#' Minimal enclosing circle
getMinCircle <- function(xy) {
  stopifnot(is.matrix(xy), is.numeric(xy), nrow(xy) >= 2, ncol(xy) == 2)
  
  H    <- chull(xy)      # convex hull indices (vertices ordered clockwise)
  hPts <- xy[H, ]        # points that make up the convex hull
  
  ## min circle may touch convex hull in only two points
  ## if so, it is centered between the hull points with max distance
  maxPD  <- getMaxPairDist(hPts)
  idx    <- maxPD$idx    # index of points with max distance
  rad    <- maxPD$d / 2  # half the distance -> radius
  rangeX <- c(hPts[idx[1], 1], hPts[idx[2], 1])
  rangeY <- c(hPts[idx[1], 2], hPts[idx[2], 2])
  ctr    <- c(rangeX[1] + 0.5*diff(rangeX), rangeY[1] + 0.5*diff(rangeY))
  
  ## check if circle centered between hPts[pt1Idx, ] and hPts[pt2Idx, ]
  ## contains all points (all distances <= rad)
  dst2ctr <- dist(rbind(ctr, hPts[-idx, ]))      # distances to center
  if(all(as.matrix(dst2ctr)[-1, 1] <= rad)) {    # if all <= rad, we're done
    tri <- rbind(hPts[idx, ], ctr)
    return( c( getCircleFrom3(tri), list(three=hPts[idx,], index=H[idx]) ) )
  }
  
  ## min circle touches hull in three points - Skyum algorithm
  S <- H                               # copy of hull indices that will be changed
  while(length(S) >= 2) {
    n    <- length(S)                # number of remaining hull vertices
    Sidx <- seq(along=S)             # index for vertices
    post <- (Sidx %% n) + 1          # next vertex in S
    prev <- Sidx[order(post)]        # previous vertex in S
    mIdx <- getMaxRad(xy, S)         # idx for maximum radius
    
    ## triangle where mIdx is vertex B in ABC
    Smax <- rbind(xy[S[prev[mIdx]], ], xy[S[mIdx], ], xy[S[post[mIdx]], ])
    
    ## if there's only two hull vertices, we're done
    if(n <= 2) { break }
    
    ## check if angle(ABC) is > 90
    ## if so, eliminate B - if not, we're done
    if(isBiggerThan90(Smax)) { S <- S[-mIdx] } else { break }
  }
  
  return( c( getCircleFrom3(Smax), list(three=Smax, index=S[c(prev[mIdx],mIdx,post[mIdx])]) ) )
}
}

plotMDS<-function(chr,chrnum){
  #learned i need to be much stricter with geno filtering (>0.9 for pca)
  chr<-read_tped("LD_lostruct_mgeno0.01.tped",chrom = paste0("scaffold_",chrnum))
  #chr<-read_tped("allpops_outgroup_nodupes_pruned.tped",chrom = paste0("scaffold_",chrnum))
  #bps<-subset(read.table("lostruct_indep_allpops_allsamples_maf_depth_gq_geno_filt_thin0.1.coords.txt",h=F),
  #            V1==paste0("scaffold_",chrnum)) [,2]
  #chr<-read_tped("~/Documents/StA/popgen/Toc_RAD.tped",chrom = paste0("scaffold_",chrnum))
  #chr<-read_tped("~/Documents/StA/popgen/Toc_RAD.tped",chrom = paste0("Scaffold_",chrnum))#for leeban data only
  chr<-data.frame(chr)
  #chr<-chr[, !colSums(is.na(chr)) > nrow(chr)/10]  
  chr <- chr %>% 
     mutate_at(c(1:ncol(chr)), ~replace_na(.,0))
  #eigenstuff.chr <- (eigen_windows(chr, win=100, k=2))
  #pc1.full[seq(chrnum*300,chrnum*300*300,by=1),]$PC1<-t(eigen_windows(chr, win=nrow(chr), k=2)[,c(4:303)])
  
  #replace NA with 0
  
  eigenstuff.chr <- (eigen_windows(chr, win = 50, k=2))
  windist <- pc_dist( eigenstuff.chr, npc=2 )
  
  fit2d <- cmdscale( windist, eig=TRUE, k=2 )
  na.inds <- is.na(eigenstuff.chr[1,] )
  chr.mds <- cmdscale( windist, eig=TRUE, k=4 )
  mds.coords <- fit2d$points
  colnames(mds.coords) <- paste("MDS coordinate", 1:ncol(mds.coords))
  #image(windist)
  #pairs( mds.coords, pch=20)
  
  mincirc <- enclosing_circle( mds.coords[,1:2] )
  mds.corners <- corners( mds.coords[,1:2], prop=.05 )
  corner.cols <- c("red","blue","purple")
  ccols <- rep("black",nrow(mds.coords))
  for (k in 1:ncol(mds.corners)) {
    ccols[ mds.corners[,k] ] <- corner.cols[k]
  }
  mds.coords<-data.frame(mds.coords)
  mds.coords$col<-ccols
  mds.coords$chr<-chrnum
  mds.coords$window<-c(1:nrow(mds.coords))
  #mds.coords$bp<-bps[seq(1,length(bps),by=100)][-1]
  return(mds.coords)
}

mds.full<-data.frame(MDS.coordinate.1=NA,MDS.coordinate.2=NA,col=NA,chr=NA,window=NA)#,bp=NA)
for (chrnum in c(1:14)){
  chr<-paste0('chr',chrnum)
  mds.full1<-plotMDS(chr,chrnum)
  mds.full<-rbind(mds.full,mds.full1)
  
}
mds.full<-mds.full[-1,]
library(reshape2)
mds.full1<-melt(mds.full,measure.vars = c('MDS.coordinate.1','MDS.coordinate.2'))
nrow(mds.full1)
mds.full1$variable1<-"MDS.2"
mds.full1[c(1:(nrow(mds.full1)/2)),]$variable1<-'MDS.1'
write.csv(file='allpops_lostruct_MDS.csv',mds.full1)

unique(mds.full1$col)
colour_scheme<-c('#5e82b5','#a060db','#c24a60')
mds.plot<-ggplot(mds.full1,aes(x=window,y=value))+
  #geom_hline(yintercept=0,linewidth=0.5,colour='#888888',linetype='dashed')+
  facet_grid(variable1~chr,scales='free')+
  theme_bw()+
  geom_point(size=0.1,colour="#888888")+
  #geom_point(data=mds.full1[!mds.full1$col=="black",],
   #          aes(x=window,y=value,colour=col),size=0.25)+
  #geom_smooth(method='loess',span=0.2,colour='black',linewidth=0.5,se=FALSE)+
  scale_colour_manual(values=colour_scheme)+
  xlab("Position")+ylab("MDS value")+
  theme(panel.grid=element_blank(),legend.position='none',
        axis.text=element_blank(),
        strip.text = element_text(margin = margin(0.05,0.05,0.05,0.05, "cm")),
        axis.ticks=element_blank(),strip.background=element_rect(fill='#eeeeee'))

setwd("/Users/jackrayner/Documents/StA/popgen/lostruct/aus_toc/")
pca.full=data.frame(sample=NA,PC1=NA,PC2=NA,chr=NA,region=NA)

plotPCA<-function(chr,chrnum){
  chr<-read_tped("LD_lostruct_mgeno0.01.tped",chrom = paste0("scaffold_",chrnum))
  chr<-data.frame(chr)
  chr <- chr %>% 
    mutate_at(c(1:ncol(chr)), ~replace_na(.,0))
  PCs<-(eigen_windows(chr, win=nrow(chr), k=2))
  PC1<-PCs[grep("PC_1",dimnames(PCs)[[2]])]
  PC2<-PCs[grep("PC_2",dimnames(PCs)[[2]])]
  
  return(data.frame(sample=c(1:ncol(chr)),PC1=PC1,PC2=PC2,chr=chrnum,
         region=c(rep("Hawaii",(ncol(chr)-7)),rep("Aus.",7))))
}

for (chrnum in c(1:14)){
  chr<-paste0('chr',chrnum)
  pca.temp<-plotPCA(chr,chrnum)
  pca.full<-rbind(pca.full,pca.temp)
  
}
pca.full<-pca.full[-1,]
}

# #plot aus samples
# pca.full<-read.csv(file='~/Documents/StA/popgen/lostruct/lostruct_PCA.csv',h=T)
# length(levels(factor(pca.full$sample)))
# nrow(pca.full)/306
info<-read.table("LD_lostruct_mgeno0.01.tfam",h=F,sep = "")
# #info<-data.frame(info[-c(301:304),])
# info<-data.frame(info[-c(307:310),])
info$region<-"Kauai"
info[grep("Hilo",info$V1),]$region<-"Hawaii"
info[grep("Oahu",info$V1),]$region<-"Oahu"
# info[grep("Aus",info$V1),]$region<-"Australia"
pca.full$pop<-rep(info$V1,14)
pca.full$region<-rep(info$region,14)
pca.full$sample<-rep(info$V2,14)
# pca.full$fw<-rep(info$V6,14)
# write.csv(file='lostruct_PCA.csv',pca.full[!is.na(pca.full$chr),])
# 
pca.full$region<-factor(pca.full$region,levels=c("Kauai","Oahu","Hawaii"))
pca.plot<-ggplot(pca.full[!is.na(pca.full$chr),],aes(x=PC1,y=PC2))+
  theme_bw()+
  #geom_point(size=0.5,aes(colour=scale(PC1)))+#facet_grid(.~chr,scales='free')+
  geom_point(size=0.75,alpha=0.5,aes(colour=region))+#facet_grid(.~chr,scales='free')+
  facet_grid(region~chr,scales='free')+
  theme(panel.grid=element_blank(),axis.text=element_blank(),
    legend.position='none',
    axis.ticks=element_blank(),
    strip.background=element_rect(fill='#eeeeee'),
    strip.text = element_text(margin = margin(0.05,0.05,0.05,0.05, "cm")))+
  #scale_colour_manual(values=c("#ce7728","#7caf5f","#a254a8","#3089c6"))+
  #scale_colour_manual(values=c("#7caf5f","#a254a8","#3089c6"))+
  #scale_colour_viridis(discrete = TRUE,option='A')+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  scale_colour_manual(values=c("#00c08b","#c77cff","#f8766d"))

ggsave('inversions_PCA.png',plot=pca.plot,
       dpi=600,height=1.4,width=10)

mds.plot+pca.plot+plot_layout(nrow=2,heights=c(2,3))

ggsave('inversions_MDS_PCA.png',
       plot=mds.plot+pca.plot+plot_layout(nrow=2,heights=c(2,3)),
       dpi=600,height=4,width=8)


