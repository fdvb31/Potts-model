install.packages("igraph")
install.packages("IsingSampler")
install.packages("bayess")
install.packages("sigmoid")
library(bayess)
library(igraph)
library(IsingSampler)
library(sigmoid)
par(pty="s")
library(grid)


cols <- c('1'="Pink","2"="Purple")
par(mfrow=c(1,1))

#Function for creating a random 2 state "len" by "wid" grid
set_matrix2 <- function(len,wid){
  ndata <- len*wid
  matrix(data=sample(c(1,-1),ndata,replace=T),len,wid)
}


#Function for Gibbs Sampler for Potts Model
gibb_lattice <- function(len,wid,beta,nruns){
  lattice <- set_matrix2(len,wid)
  lattices <- list(rep(0,nruns+1))
  lattices[[1]] <- lattice
  for (a in 1:nruns){
    for (i in 1:len){
      for (j in 1:wid){
        #Set Periodic Boundary conditions
        up <- ifelse(i == 1, len, i - 1)
        down <- ifelse(i == len, 1, i + 1)
        left <- ifelse(j == 1, wid, j - 1)
        right <- ifelse(j == wid, 1, j + 1)
        self <- lattice[i,j]
        neighbours <- c(lattice[up,j],lattice[down,j],lattice[i,left],lattice[i,right])
        #Calculate Enegery for the state being in 1 or -1
        E1 <- sum(1==neighbours)
        E2 <- sum(-1==neighbours)
        E_change <- E1-E2
        #Updates the vertex [i,j] using Gibbs sampler
        Prob <- exp(beta*E_change)
        Prob <- Prob/(1+Prob)
        U <- runif(1,0,1)
        if (Prob>U){
          lattice[i,j] <- 1
        } else{
          lattice[i,j] <- -1
        }
      }
    }
    lattices[[a+1]] <- lattice
  }
  return(lattice)
}

#Sample run with len x wid =200x200, beta=log(1+sqrt(2)), number of runs =200
g <- (gibb_lattice(200,200,log(1+sqrt(2)),200))

#display the grid
par(mar=rep(0, 4), xpd = NA) library(grid)
image(1:nrow(g), 1:ncol(g), as.matrix(g), col=cols,bty ="n",axes=F,frame.plot=F, xaxt='n', ann=FALSE, yaxt='n', asp=745/962)



#Run using the "isinghm" function 
image(1:200,1:200,isinghm(200,200,200,beta=0.1),
      col=cols,bty ="n",axes=F,frame.plot=F, xaxt='n', ann=FALSE, yaxt='n', asp=745/962)
