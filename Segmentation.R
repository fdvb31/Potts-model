library(PottsUtils)
library(imager)
library(dplyr)
library(tidyr)
library(gasper)
library(combinat)
library("imagerExtra")
library("tvR")

#Function for loading and k-means clustering an image
imageloader <- function(filepath,q){
  #Loads and transform image to dataframe
  z <- load.image(filepath)
  plot(z)
  zdata <- as.data.frame(z)
  height(z)
  zwider <- pivot_wider(zdata,names_from="cc",values_from="value")
  #Extracts Colour values
  rgbs <- zwider[,c(3,4,5)]
  rgbs <- as.matrix(rgbs)
  #Applies k-Means to the colours with q clusters
  kmean2 <- kmeans(rgbs,centers =q,nstart=25)
  colours2 <- as.data.frame(kmean2[2])
  clusters2 <- as.data.frame(kmean2[1])
  ext2 <- clusters2
  ext2 <- as.matrix(ext2)
  dim(ext2) <- c(width(z),height(z))
  #Returns the matrix of group of each vertex with the colour that each group represents
  return(list(ext2,colours2))
}

#Function for loading and segmenting an image
segmenting <- function(lattice,nruns,beta_,var,q){
  #Gets the necessary grid,block and neighbourhood structure
  wid <- height(lattice)
  len <- width(lattice)
  temp <- rep(0,len*wid)
  neighbors <- getNeighbors(mask=matrix(1, len, wid), neiStruc=c(2,2,0,0))
  blocks <- getBlocks(mask=matrix(1, len, wid), nblock=2)
  pottscol <- as.matrix(rep(0,length(lattice)*q,dim=c(length(lattice),q)))
  dim(pottscol) <- c(length(lattice),q)
  #Calculates the external field
  for(k in 1:q){
    pottscol[,k] <- -(0.5/(var))*((rgbs[,1]-colours2[k,1])^2+(rgbs[,2]-colours2[k,2])^2+(rgbs[,3]-colours2[k,3])^2)
  }
  cols <- as.integer(lattice)
  #Runs the Pott model nruns time
    for (j in 1:nruns){
      cols <- as.integer(cols)
      print(j)
      cols<- rPotts1(nvertex=(len*wid), ncolor=q, neighbors=neighbors, blocks=blocks,
                     spatialMat=NULL,,external=pottscol,beta=beta_,colors = cols)  }
    cols <- as.matrix(cols)
    dim(cols) <- c(width(z),height(z))
  #Returns the final grid
}
