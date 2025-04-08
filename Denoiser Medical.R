library(DescTools)
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

#Function that simulates the Potts mode with an external field, with the modal grid caculate after each iteration
mode_gibbs_external2 <- function(lattice,beta_,nruns=1,q,colours,ncomplete,var,cols_acc,colours2,rgb_noise){
  wid <- height(lattice)
  len <- width(lattice)
  temp <- rep(0,len*wid*ncomplete)
  runs <- array(temp,c(len,wid,ncomplete))
  neighbors <- getNeighbors(mask=matrix(1, len, wid), neiStruc=c(2,2,0,0))
  blocks <- getBlocks(mask=matrix(1, len, wid), nblock=2)
  pottscol <- as.matrix(rep(0,length(lattice)*q,dim=c(length(lattice),q)))
  dim(pottscol) <- c(length(lattice),q)
  for(k in 1:q){
    pottscol[,k] <- -(0.5/var)*((rgbs_noise[,1]-colours[k,1])^2+(rgbs_noise[,2]-colours[k,2])^2+(rgbs_noise[,3]-colours[k,3])^2)
  }
  cols <- as.integer(lattice)
  for (i in 1:ncomplete){
    runs[,,i] <- cols
  }
  snrs <- rep(SNR(cols_acc,cols),ncomplete+1)
  psnrs <- rep(PSNR(cols_acc,cols),ncomplete+1)
  for (j in 1:nruns){
    for (i in 2:ncomplete){
      runs[,,i]<- rPotts1(nvertex=(len*wid), ncolor=q, neighbors=neighbors, blocks=blocks,
                          spatialMat=NULL,beta=beta_,external=pottscol,colors = as.integer(runs[,,i]))
    }
    temp <- rep(0,len*wid)
    mode <- array(temp,c(len,wid))
    for (l in 1:len){
      for (w in 1:wid){
        mode[l,w] <- Mode(runs[l,w,])[1]
      }
    }
  }
  return(list(mode))
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

#Makes sure both grayscale lists represent the same colour
colour_correct_gray <- function(colours,ext,perm){
  x <- matrix(,nrow = nrow(colours), ncol = 2)
  for (i in 1:nrow(colours)){
    x[i,]<-rbind(c(i,as.integer(perm[i])))
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

#Transforms for creating cimgs
cimglonger <- function(zwider2,rgbs){
  zwider2[,c(3,4,5)] <- rgbs
  zlong2 <- pivot_longer(zwider2,c("1","2","3"),names_to="cc",values_to="value")
  zlong2 <- arrange(zlong2,cc)
  zlong2$cc <- as.integer(zlong2$cc)
  ztv <- as.cimg(zlong2)
  return(ztv)
}


#Matrix transformation for displaying images

flipvert <- function(matrix){
  ext2 <- as.data.frame(matrix)
  ext2 <- rev(ext2)
  ext2 <- as.matrix(ext2)
  return(ext2)
}

#Adds external noise to an grayscale image

noisegrayscale <- function(filepath,mu=0,var=0.25){
  z <- load.image(filepath)
  z <- resize(z,256,256)
  z <- grayscale(z)
  zdata <- as.data.frame(z)
  grays <- zdata[,3]
  grays <- as.matrix(grays)
  temp <- as.matrix(rnorm(length(grays),grays,sqrt(var)))
  graynoise <- as.matrix(apply(temp, 1, function(j) min(max(j, 0), 1)))
  dim(graynoise) <- c(dim(grays))
  return(list(z,grays,graynoise))
}

#Adds multicplicative external noise to an grayscale image

noisegrayscalemulti <- function(filepath,mu=0,var=0.25){
  z <- load.image(filepath)
  z <- resize(z,256,256)
  z <- grayscale(z)
  zdata <- as.data.frame(z)
  grays <- zdata[,3]
  grays[which(grays==0)] <- 0.00000000001
  grays <- as.matrix(grays)
  temp <- as.matrix(rlnorm(length(grays),log(grays),sqrt(var)))
  graynoise <- as.matrix(apply(temp, 1, function(j) min(max(j, 0), 1)))
  dim(graynoise) <- c(dim(grays))
  return(list(z,grays,graynoise))
}


#K-means clustering for a grayscale image

kmeansGray <- function(gray,graynoise,q,grayz){
  kmean <- kmeans(graynoise,centers =q,nstart=25)
  colours <- as.data.frame(kmean[2])
  clusters <- as.data.frame(kmean[1])
  ext <- clusters
  ext <- as.matrix(ext)
  dim(ext) <- c(width(grayz),height(grayz))
  kmean2 <- kmeans(grays,centers =q,nstart=25)
  colours2 <- as.data.frame(kmean2[2])
  clusters2 <- as.data.frame(kmean2[1])
  ext2 <- clusters2
  ext2 <- as.matrix(ext2)
  dim(ext2) <- c(width(grayz),height(grayz))
  dissmall <- 1000000
  for (x in permn(c(1:q))){
    dis <- (sqrt(sum((colours2[as.vector(x),]-colours)^2)))
    if (dis<dissmall){
      dissmall <- dis
      perm <- as.vector(x)
    }
  }
  colours2 <- as.matrix(colours2[perm,])
  ext2 <- colour_correct_gray(colours2,ext2,perm)
  return(list(ext, ext2,colours,colours2))
}

#Function that simulates the Potts mode with an external field for a grayscale image, with the modal grid calculated at the end

mode_gibbs_external_gray <- function(lattice,beta_,nruns=1,q,colours,ncomplete,var,cols_acc,colours2,noise){
  wid <- height(lattice)
  len <- width(lattice)
  temp <- rep(0,len*wid*ncomplete)
  runs <- array(temp,c(len,wid,ncomplete))
  neighbors <- getNeighbors(mask=matrix(1, len, wid), neiStruc=c(2,2,0,0))
  blocks <- getBlocks(mask=matrix(1, len, wid), nblock=2)
  pottscol <- as.matrix(rep(0,length(lattice)*q,dim=c(length(lattice),q)))
  dim(pottscol) <- c(length(lattice),q)
  for(k in 1:q){
    pottscol[,k] <- -(0.5/var)*(noise-as.matrix(colours)[k])^2
  }  
  for (i in 1:ncomplete){
    cols <- as.integer(lattice)
    psnr <- c(PSNR(cols_acc,cols))
    for (j in 1:nruns){
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



#Function that simulates the Potts mode with an external field for a grayscale image, with the modal grid caculate after each iteration

mode_gibbs_external_gray2 <- function(lattice,beta_,nruns=1,q,colours,ncomplete,var,cols_acc,colours2,noise){
  wid <- height(lattice)
  len <- width(lattice)
  temp <- rep(0,len*wid*ncomplete)
  runs <- array(temp,c(len,wid,ncomplete))
  neighbors <- getNeighbors(mask=matrix(1, len, wid), neiStruc=c(2,2,0,0))
  blocks <- getBlocks(mask=matrix(1, len, wid), nblock=2)
  pottscol <- as.matrix(rep(0,length(lattice)*q,dim=c(length(lattice),q)))
  dim(pottscol) <- c(length(lattice),q)
  for(k in 1:q){
    pottscol[,k] <- -(0.5/var)*(noise-as.matrix(colours)[k])^2
  }  
  cols <- as.integer(lattice)
  for (i in 1:ncomplete){
    runs[,,i] <- cols
  }
  snrs <- rep(SNR(cols_acc,cols),ncomplete+1)
  psnrs <- rep(PSNR(cols_acc,cols),ncomplete+1)
  for (j in 1:nruns){
    for (i in 2:ncomplete){
      runs[,,i]<- rPotts1(nvertex=(len*wid), ncolor=q, neighbors=neighbors, blocks=blocks,
                          spatialMat=NULL,beta=beta_,external=pottscol,colors = as.integer(runs[,,i]))
    }
    temp <- rep(0,len*wid)
    mode <- array(temp,c(len,wid))
    for (l in 1:len){
      for (w in 1:wid){
        mode[l,w] <- Mode(runs[l,w,])[1]
      }
    }
    psnrs <- rbind(psnrs,c(apply(runs,3, function(x) PSNR(cols_acc,x)),PSNR(cols_acc,mode)))
    snrs <- rbind(snrs,c(apply(runs,3, function(x) SNR(cols_acc,x)),SNR(cols_acc,mode)))
  }
  return(list(mode,psnrs,snrs))
}

file.path4 <- file.path()#File Path is Inserted here

#Calculate the PSNR, SNR, and Entropy for different values of var and beta for the Potts Model
mat1 <- matrix(rep(0,25),nrow=5)
mat2 <- matrix(rep(0,25),nrow=5)
mat3 <- matrix(rep(0,25),nrow=5)
i_count <- 1
j_count <- 1
for (i in c(0.05,0.1,0.15,0.2,0.25)){
  graybase <- noisegrayscalemulti2(file.path4,0,var=i^2)
  zgray <- graybase[[1]]
  grays <- graybase[[2]]
  graynoise <- graybase[[3]]
  
  graybase2 <- kmeansGray(grays,graynoise,4,zgray)
  ext <- graybase2[[1]]
  colours <- graybase2[[3]]
  ext2 <- graybase2[[2]]
  colours2 <- graybase2[[4]]
  for (j in c(0.1,0.5,1,2.5,5)){
    print(i_count)
    print(j_count)
    gibbs <- mode_gibbs_external_gray(ext,j,15,4,log(colours),25,i^2,ext2,colours2,log(graynoise))
    mat1[i_count,j_count] <- PSNR(gibbs,ext2)
    mat2[i_count,j_count] <- SNR(gibbs,ext2)
    mat3[i_count,j_count] <- Entropy(gibbs)
    j_count <- j_count+1
  }
  i_count <- i_count+1
  j_count <- 1
}

#Calculate the PSNR, SNR, and Entropy for different values of var for the other denoising aglorithms
mat4 <- matrix(rep(0,30),nrow=5)
mat5 <- matrix(rep(0,30),nrow=5)
mat6 <- matrix(rep(0,30),nrow=5)
i_count <- 1
for (i in c(0.01,0.02,0.05,0.1,0.15)){
  print(i_count)
  graybase <- noisegrayscale(file.path4,0,var=i^2)
  zgray <- graybase[[1]]
  grays <- graybase[[2]]
  graynoise <- graybase[[3]]
  
  graybase2 <- kmeansGray(grays,graynoise,4,zgray)
  ext <- graybase2[[1]]
  colours <- graybase2[[3]]
  ext2 <- graybase2[[2]]
  colours2 <- graybase2[[4]]
  
  noisy <- matrix(as.numeric(as.matrix(colours[ext,])))
  grays <- noisy
  newgray <- as.data.frame(zgray)
  newgray$value <- noisy
  newgray <- as.cimg(newgray)
  
  clean <- matrix(as.numeric(as.matrix(colours2[ext2,])))
  cleangray<- as.data.frame(zgray)
  cleangray$value <- clean
  cleangray <- as.cimg(cleangray)
  zgray <- newgray
  tv <- denoise2(grays,niter=100)
  dct <- DenoiseDCT(zgray,i)
  zgrayM <- as.matrix(zgray)
  cleangrayM <- as.matrix(cleangray)
  Median <- denoise(zgrayM,type="median")
  Lee <- denoise(zgrayM,type="Lee")
  Kuan <- denoise(zgrayM,type="Kuan",looks=5)
  Nathan <- denoise(zgrayM,type="Nathan")
  mat4[i_count,] <- c(PSNR(clean,tv),PSNR(dct,cleangray),PSNR(Median,cleangrayM),PSNR(Lee,cleangrayM),PSNR(Kuan,cleangrayM),PSNR(Nathan,cleangrayM))
  mat5[i_count,] <- c(SNR(clean,tv),SNR(dct,cleangray),SNR(Median,cleangrayM),SNR(Lee,cleangrayM),SNR(Kuan,cleangrayM),SNR(Nathan,cleangrayM))
  mat6[i_count,] <- c(Entropy(tv),Entropy(dct),Entropy(Median),Entropy(Lee),Entropy(Kuan),Entropy(Nathan))
  i_count <- i_count+1
}

