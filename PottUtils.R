install.packages("PottsUtils")
install.packages("R2jags")
install.packages("MASS")
install.packages("mvtnrom")
library(PottsUtils)
library(?R2jags)
library(imager)
library(dplyr)
library(tidyr)
library(MASS)
library(mvtnorm)

ibrary(grid)

cols <- c("1"="cornflowerblue","2"="lightblue","3"="white")

#Function for creating a len  x wid simulate of the Potts model using the "PottsUtil" package
gibbs <- function(len,wid,q,beta_,nruns){
  neighbors <- getNeighbors(mask=matrix(1, len, wid), neiStruc=c(2,2,0,0))
  blocks <- getBlocks(mask=matrix(1, len, wid), nblock=2)
  x <- BlocksGibbs(n=nruns, nvertex=(len*wid), ncolor=q, neighbors=neighbors, blocks=blocks,
                   spatialMat=NULL, beta=beta_)
  y <- x[,ncol(x)]
  dim(y) <- c(len,wid)
  return(y)
}
par(mfrow=c(1,1))
y <- gibbs(500,500,3,(log(1+sqrt(3))),200)
par(mar=rep(0, 4), xpd = NA) 
image(1:nrow(y), 1:ncol(y), as.matrix(y), col=cols,bty ="n",axes=F,frame.plot=F, xaxt='n', ann=FALSE, yaxt='n', asp=745/962)
Axis(side=1, labels=FALSE)
Axis(side=2, labels=FALSE)