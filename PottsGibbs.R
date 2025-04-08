
cols <- c("1"="Pink","2"="Purple","3"="red")
cols <- c("1"="Pink","2"="Purple")
par(mar=rep(0, 4), xpd = NA) 
library(grid)


set_matrix <- function(len,wid,q){
  #Randomly creates a length by width matrix of q colours
  ndata <- len*wid
  matrix(data=sample(c(1:q),ndata,replace=T),len,wid)
}

potts_lattice_gibbs <- function(len,wid,beta,nruns,q){
  lattice <- set_matrix(len,wid,q) #Creates a grid to use as the starting conditions
  lattices <- list(rep(0,nruns+1))
  lattices[[1]] <- lattice
  for (a in 1:nruns){
    for (i in 1:len){
      for (j in 1:wid){
        #Computes the neighbours of vertex [i,j] a using a periodic boundary
        up <- ifelse(i == 1, len, i - 1)
        down <- ifelse(i == len, 1, i + 1)
        left <- ifelse(j == 1, wid, j - 1)
        right <- ifelse(j == wid, 1, j + 1)
        self <- lattice[i,j]
        neighbours <- c(lattice[up,j],lattice[down,j],lattice[i,left],lattice[i,right])
        energies <- c(rep(0,q))
        for (x in 1:q){
          #For each colour, q, computes the energy given the neighbours of vertex [i,j]
          energies[x] <- beta*sum(x==neighbours)
        }
        energies <- exp(beta*energies)
        Prob <- energies/sum(energies) #Normalises the probailities
        lattice[i,j] <- sample(c(1:q),1,replace=FALSE,prob=Prob) 
        #Updates the vertex with the calculated probailites for each colour q
        }
      }
  }
  image(1:nrow(lattice), 1:ncol(lattice), as.matrix(lattice), col=cols,bty ="n",axes=F,frame.plot=F, xaxt='n', ann=FALSE, yaxt='n', asp=745/962)
}



cols <- c("1"="purple","2"="pink","3"="red")
#potts_lattice(length,width,beta,#number_of_runs,q)
potts_lattice_gibbs(200,200,0.99,200,3)



