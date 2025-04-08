par(mar=rep(0, 4), xpd = NA) 
library(grid)
set_matrix <- function(len,wid,q){ #Assigns a value of {1,..,q} to every cell of the matrix
  ndata <- len*wid
  matrix(data=sample(c(1:q),ndata,replace=T),len,wid)
}

potts_lattice_MCMC <- function(len,wid,beta,nruns,q){
  lattice <- set_matrix(len,wid,q) #Creates a grid to use as the starting conditions
  for (a in 1:nruns){
    for (i in 1:len){
      for (j in 1:wid){
        #Computes the neighbours of vertex [i,j] using a periodic boundary
        up <- ifelse(i == 1, len, i - 1)
        down <- ifelse(i == len, 1, i + 1)
        left <- ifelse(j == 1, wid, j - 1)
        right <- ifelse(j == wid, 1, j + 1)
        self <- lattice[i,j]
        neighbours <- c(lattice[up,j],lattice[down,j],lattice[i,left],lattice[i,right])
        E_now <- sum(self==neighbours) #Computes the current energy of vertex [i,j]
        prop <- sample(c(1:q),size = 1)
        E_new <- sum(prop==neighbours) #Computes the energy of a randomly selected colour q
        E_change <- E_new-E_now
        #Updates the vertex [i,j] using a MCMC algorithm
        U <- runif(1,0,1)
        Prob <- exp(beta*E_change)
        if (Prob>U){
          lattice[i,j] <- prop
        }
      }
    }
  }
  #Displays the final grid
  image(1:nrow(lattice), 1:ncol(lattice), as.matrix(lattice), col=cols,bty ="n",axes=F,frame.plot=F, xaxt='n', ann=FALSE, yaxt='n', asp=745/962)
}


cols <- c("1"="purple","2"="pink","3"="red","4"="violet","5"="darkred")
#potts_lattice(length,width,beta,#number_of_runs,q)
potts_lattice_MCMC(200,200,log(1+sqrt(5)),200,5)


