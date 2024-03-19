#library(MCMCpack)
library(mvtnorm)

# Proposal (and density) for betas
update.betaj <- function(j,incidence,sigma2j,TN){
  
  n <- ncol(incidence) # number of variables
  
  colnames(incidence) <- rownames(incidence) <- colnames(data)
  
  A1j <- names(incidence[,j][incidence[,j] != 0])
  A0j <- colnames(incidence)[!colnames(incidence) %in% A1j]
  
  parentnodes <- which(incidence[,j]==1)
  lp <- length(parentnodes)
  
  # Design matrix
  if(lp == 0){
    betaj <- t(as.matrix(rep(0,n)))
    colnames(betaj) <- colnames(incidence)
  } else {
    mu = solve(TN[parentnodes,parentnodes]) %*% TN[parentnodes, j]
    Sigma = sigma2j * solve(TN[parentnodes,parentnodes])
    
    # Proposal for betas
    beta1j <- mvtnorm::rmvnorm(1, mu, Sigma)
    beta0j <- rep(0,length(A0j))
    betaj <- cbind(beta1j,t(beta0j))
    colnames(betaj) <- c(A1j,A0j)
    betaj <- betaj[,sort(colnames(betaj))]
  }
  
  return(betaj)
}

density.betaj <- function(j,incidence,values,sigma2j,TN){
  
  n <- ncol(incidence) # number of variables

  parentnodes <- which(incidence[,j]==1)
  lp <- length(parentnodes)
  
  # Design matrix
  if(lp == 0){
    logdbeta <- 0
  } else {
    mu = solve(TN[parentnodes,parentnodes]) %*% TN[parentnodes, j]
    Sigma = sigma2j * solve(TN[parentnodes,parentnodes])
    
    logdbeta <- mvtnorm::dmvnorm(values[which(values!=0)], mu, Sigma, log = TRUE)
  }
  
  return(logdbeta)
}


# Proposal (and density) for sigma2

update.sigma2j <- function(j,N,incidence,am=1,TN){
  
  n <- ncol(incidence)

  aw <- n + am + 1
  t <- am * (aw - n - 1) / (am + 1)
  
  parentnodes <- which(incidence[,j]==1)
  lp <- length(parentnodes)
  
  a <- (aw - n + lp + 1)/2
  b <- t/2
  
  aj <- N/2 + a
  
  # Design matrix
  if(lp == 0){
    bj <- 1/2*TN[j,j]
  } else {
    bj <- 1/2*(TN[j,j] - t(TN[parentnodes,j]) %*% solve(TN[parentnodes,parentnodes]) %*% TN[parentnodes,j])
  }
  
  # Proposal for sigma2
  sigmaj <- invgamma::rinvgamma(1, shape = aj, rate = bj)
  return(sigmaj)
}


density.sigma2j <- function(j,N,incidence,values,am=1,TN){
  
  n <- ncol(incidence)
  
  aw <- n + am + 1
  t <- am * (aw - n - 1) / (am + 1)

  parentnodes <- which(incidence[,j]==1)
  lp <- length(parentnodes)
  
  a <- (aw - n + lp + 1)/2
  b <- t/2
  
  aj <- N/2 + a
  
  # Design matrix
  if(lp == 0){
    bj <- 1/2*TN[j,j]
  } else {
    bj <- 1/2*(TN[j,j] - t(TN[parentnodes,j]) %*% solve(TN[parentnodes,parentnodes]) %*% TN[parentnodes,j])
  }
  # Calculating densities given sigma2 values
  logdsig <- invgamma::dinvgamma(values, shape = aj, rate = bj, log=TRUE)
  return(logdsig)
}


logdet <- function(x){
  return(2 * sum(log(diag(chol(x)))))
}


pMVN <- function(incidence,values,j,t,sigma2j,log=FALSE){
  
  parentnodes <- which(incidence[,j]==1)
  lp <- length(parentnodes)
  
  # Design matrix
  if(lp == 0){
    logdbeta <- 0
  } else {
    Lambda0 <- diag(rep(t, lp), lp)
    logdbeta <- mvtnorm::dmvnorm(values[which(values != 0)],  
                                 sigma = as.matrix(sigma2j*solve(Lambda0)), log = TRUE)
  }
  
  if (log == FALSE) logdbeta <- exp(logdbeta)
  
  return(logdbeta)
}


# Prior distribution for sigma^2
pIG <- function(x, a, b, log=FALSE){
  
  dens <- invgamma::dinvgamma(x, shape=a, rate=b, log=TRUE)
  
  if (log == FALSE) dens <- exp(dens)
  
  return(dens)
}

psigma2 <- function(x, a, b, log=FALSE){
  
  dens <- sum(invgamma::dinvgamma(x, shape=a, rate=b, log=TRUE))
  
  if (log == FALSE) dens <- exp(dens)
  
  return(dens)
}


loglikelihood <- function(j,data,incidence,betas,sigma2s,am=1){
  
  colnames(incidence) <- rownames(incidence) <- colnames(data)
  nobvs <- dim(data)[1]
  #n <- dim(data)[2]
  
  yj <- data[,j]
  
  A1j <- c(names(incidence[,j][incidence[,j] != 0]))
  #cat("Parents of ",j,"are", names(A1j))
  
  if (length(A1j) != 0){
    # Design matrix
    Xj <- as.matrix(data[,A1j])
    #colnames(Xj) <- c(colnames(data)[-j])
    Beta_j <- betas[which(betas != 0)]
    
    deltaXj <- yj - Xj %*% Beta_j
  } else {
    deltaXj <- yj
  }
  
  logl <- 1/2*(log(am) - log(nobvs + am)) - nobvs/2 * (log(2*pi) + log(sigma2s)) - 
    1/(2*sigma2s) * (t(deltaXj) %*% deltaXj - (1/(nobvs + am))*sum(deltaXj)^2)
  #(1/(n + am))*(t(deltaXj) %*% rep(1,nobvs))^2)
  
  return(logl)
}



