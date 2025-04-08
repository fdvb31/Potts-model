library(DescTools)
library(R2jags)
library(PottsUtils)
library(imager)
library(dplyr)
library(tidyr)
library(gasper)
library(combinat)
library(SpatialPack)
library("imagerExtra")
library("tvR")




#Function that simulates the Potts mode with an external field, with the modal grid calculated at the end
mode_gibbs_external <- function(lattice,beta_,nruns=1,q,colours,ncomplete,var,cols_acc,colours2,rgbs_noise){
  wid <- height(lattice)
  len <- width(lattice)
  temp <- rep(0,len*wid*ncomplete)
  runs <- array(temp,c(len,wid,ncomplete))
  neighbors <- getNeighbors(mask=matrix(1, len, wid), neiStruc=c(2,2,0,0))
  blocks <- getBlocks(mask=matrix(1, len, wid), nblock=2)
  pottscol <- as.matrix(rep(0,length(lattice)*q,dim=c(length(lattice),q)))
  dim(pottscol) <- c(length(lattice),q)
  for(k in 1:q){
    pottscol[,k] <- -(0.5/(var))*((rgbs_noise[,1]-colours[k,1])^2+(rgbs_noise[,2]-colours[k,2])^2+(rgbs_noise[,3]-colours[k,3])^2)
  }  
  for (i in 1:ncomplete){
    cols <- as.integer(lattice)
    psnr <- c(PSNR(cols_acc,cols))
    for (j in 1:nruns){
      print(j)
      cols<- rPotts1(nvertex=(len*wid), ncolor=q, neighbors=neighbors, blocks=blocks,
                     spatialMat=NULL,beta=beta_,external=pottscol,colors = cols)
      psnr <- c(psnr,PSNR(cols_acc,cols))
    }
    runs[,,i] <- cols
  }
  temp <- rep(0,len*wid)
  mode <- array(temp,c(len,wid))
  for (i in 1:len){
    for (j in 1:wid){
      mode[i,j] <- Mode(runs[i,j,])[1]
    }
  }
  return(mode)
}

#Makes sure both colour lists represent the same colour
colour_correct <- function(colours,ext){
  x <- matrix(,nrow = nrow(colours), ncol = 2)
  for (i in 1:nrow(colours)){
    x[i,]<-rbind(c(i,as.integer(rownames(colours[i,]))))
  }
  x <- x[order(x[,2]),]
  new_ext <- ext
  for (i in 1:length(ext)){
    new_ext[[i]] <- x[ext[[i]],1]
  }
  return(new_ext)
}

#Adds external noise to an image
get_external_noise <- function(filepath,mu=0,var=0.25){
  z <- load.image(filepath)
  zdata <- as.data.frame(z)
  height(z)
  zwider <- pivot_wider(zdata,names_from="cc",values_from="value")
  rgbs <- zwider[,c(3,4,5)]
  rgbs <- as.matrix(rgbs)
  temp <- rnorm(length(rgbs),rgbs,sqrt(var))
  rgbs_noise <- as.matrix(sapply(temp, function(i) sapply(i, function(j) min(max(j, 0), 1))))
  dim(rgbs_noise) <- c(dim(rgbs))
  return(list(z,zwider,rgbs,rgbs_noise))
}

#K-means clustering for a RGB image
kmeanG <- function(rgbs,rgbs_noise,q,z){
  kmean <- kmeans(rgbs_noise,centers =q,nstart=25)
  colours <- as.data.frame(kmean[2])
  clusters <- as.data.frame(kmean[1])
  ext <- clusters
  ext <- as.matrix(ext)
  dim(ext) <- c(width(z),height(z))
  kmean2 <- kmeans(rgbs,centers =q,nstart=25)
  colours2 <- as.data.frame(kmean2[2])
  clusters2 <- as.data.frame(kmean2[1])
  ext2 <- clusters2
  ext2 <- as.matrix(ext2)
  dim(ext2) <- c(width(z),height(z))
  dissmall <- 1000000
  for (x in permn(c(1:q))){
    dis <- (sqrt(sum((colours2[as.vector(x),]-colours)^2)))
    if (dis<dissmall){
      dissmall <- dis
      perm <- as.vector(x)
    }
  }
  colours2 <- colours2[perm,]
  ext2 <- colour_correct(colours2,ext2)
  return(list(ext, ext2,colours,colours2))
}

file.path3 <- file.path()#File Path is Inserted here
base <- get_external_noise(file.path3,0,0.15)
z <- base[[1]]
zwider <- base[[2]]
rgbs <- base[[3]]
rgbs_noise <- base[[4]]
base2 <- kmeanG(rgbs,rgbs_noise,3,z)
ext <- base2[[1]]
colours <- base2[[3]]
ext2 <- base2[[2]]
colours2 <- base2[[4]]
mode <- mode_gibbs_external(ext2,log(1+sqrt(3)),20,3,colours,20,0.15,ext2,colours2,rgbs_noise)
image(1:nrow(flipvert(mode)), 1:ncol(flipvert(mode)), as.matrix(flipvert(mode)), col=gray(as.matrix(colours2)),bty ="n",axes=F,frame.plot=F, xaxt='n', ann=FALSE, yaxt='n', asp=745/962)
Axis(side=1, labels=FALSE)
Axis(side=2, labels=FALSE)

par(mar=rep(0, 4), xpd = NA) 
image(1:nrow(flipvert(ext2)), 1:ncol(flipvert(ext2)), as.matrix(flipvert(ext2)), col=gray(as.matrix(colours2)),bty ="n",axes=F,frame.plot=F, xaxt='n', ann=FALSE, yaxt='n', asp=745/962)
Axis(side=1, labels=FALSE)
Axis(side=2, labels=FALSE)

par(mar=rep(0, 4), xpd = NA) 
image(1:nrow(flipvert(ext)), 1:ncol(flipvert(ext)), as.matrix(flipvert(ext)), col=gray(as.matrix(colours2)),bty ="n",axes=F,frame.plot=F, xaxt='n', ann=FALSE, yaxt='n', asp=745/962)
Axis(side=1, labels=FALSE)
Axis(side=2, labels=FALSE)
