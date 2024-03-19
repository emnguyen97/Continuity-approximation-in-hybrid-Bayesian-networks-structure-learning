

# Implement Tabu search (TABU) greedy search from bnlearn and return adjacency matrix
runtabu <- function(data, datmat, score){
  rtabu <- bnlearn::tabu(data, score=score)
  tabumat <- as.matrix(get.adjacency(graph.edgelist(rtabu$arcs)))
  tabumat <- as.data.frame(tabumat)[,order(colnames(as.data.frame(tabumat)))]
  tabumat <- tabumat[order(rownames(tabumat)),]
  tabumat <- as.matrix(tabumat)
  
  # If the dimension of the estimated matrix does not match the number of nodes, fill in these empty links
  nvar <- ncol(data)
  num <- ncol(tabumat)
  if (num < nvar){
    tbadded <- setdiff(colnames(datmat),colnames(tabumat))
    n <- length(tbadded)
    fillc <- matrix(0, nvar-n, n)
    fillr <- matrix(0, n, nvar)
    tabumat <- cbind(tabumat,fillc)
    tabumat <- rbind(tabumat,fillr)
    colnames(tabumat)[(num+1):(num+n)] <- rownames(tabumat)[(num+1):(num+n)] <- tbadded
    tabumat <- tabumat[LETTERS[1:nvar],LETTERS[1:nvar]]
  }
  
  return(tabumat)
}

# Implement PC algorithm from bnlearn and return adjacency matrix
runhc <- function(data, datmat, score){
  rhc <- bnlearn::hc(as.data.frame(data), score=score)
  hcmat <- as.matrix(get.adjacency(graph.edgelist(rhc$arcs)))
  hcmat <- as.data.frame(hcmat)[,order(colnames(as.data.frame(hcmat)))]
  hcmat <- hcmat[order(rownames(hcmat)),]
  hcmat <- as.matrix(hcmat)
  
  # If the dimension of the estimated matrix does not match the number of nodes, fill in these empty links
  nvar <- ncol(data)
  num <- ncol(hcmat)
  if (num < nvar){
    tbadded <- setdiff(colnames(datmat),colnames(hcmat))
    n <- length(tbadded)
    fillc <- matrix(0, nvar-n, n)
    fillr <- matrix(0, n, nvar)
    hcmat <- cbind(hcmat,fillc)
    hcmat <- rbind(hcmat,fillr)
    colnames(hcmat)[(num+1):(num+n)] <- rownames(hcmat)[(num+1):(num+n)] <- tbadded
    hcmat <- hcmat[LETTERS[1:nvar],LETTERS[1:nvar]]
  }
  
  return(hcmat)
}

# Implement PMCMC from BiDAG and return adjacency matrix

modelpcore<-function(MCMCchain, p, pdag=FALSE, burnin=0.2, DBN=FALSE, nsmall=0, dyn=0, b=0) {
  
  varlabels<-colnames(MCMCchain[[1]])
  n<-nrow(MCMCchain[[1]])
  incidence<-matrix(rep(0, n*n), nrow=n, ncol=n)
  endstep<-length(MCMCchain)
  startstep<-max(as.integer(burnin*endstep),1)
  if (pdag) {
    cpdags<-lapply(MCMCchain[startstep:endstep],dagadj2cpadj)
    incidence[which(as.matrix(Reduce('+', cpdags)/(endstep-startstep+1))>p)]<-1
  } else {
    incidence[which(as.matrix(Reduce('+', MCMCchain[startstep:endstep])/(endstep-startstep+1))>p)]<-1
  }
  colnames(incidence)<-varlabels
  rownames(incidence)<-varlabels
  if(DBN) {
    incidence<-DBNcut(incidence,dyn=dyn,b=b)
    incidence.init<-DBNinit(incidence,dyn=dyn,b=b)
    incidence[1:(dyn+b),1:(dyn+b)]<-incidence.init
  }
  return(incidence)
}

runpartMCMC <- function(data,score){
  
  modelpcore<-function(MCMCchain, p, pdag=FALSE, burnin=0.2, DBN=FALSE, nsmall=0, dyn=0, b=0) {
    
    varlabels<-colnames(MCMCchain[[1]])
    n<-nrow(MCMCchain[[1]])
    incidence<-matrix(rep(0, n*n), nrow=n, ncol=n)
    endstep<-length(MCMCchain)
    startstep<-max(as.integer(burnin*endstep),1)
    if (pdag) {
      cpdags<-lapply(MCMCchain[startstep:endstep],dagadj2cpadj)
      incidence[which(as.matrix(Reduce('+', cpdags)/(endstep-startstep+1))>p)]<-1
    } else {
      incidence[which(as.matrix(Reduce('+', MCMCchain[startstep:endstep])/(endstep-startstep+1))>p)]<-1
    }
    colnames(incidence)<-varlabels
    rownames(incidence)<-varlabels
    if(DBN) {
      incidence<-DBNcut(incidence,dyn=dyn,b=b)
      incidence.init<-DBNinit(incidence,dyn=dyn,b=b)
      incidence[1:(dyn+b),1:(dyn+b)]<-incidence.init
    }
    return(incidence)
  }
  
  myScore<-BiDAG::scoreparameters(score, data)
  
  partfit<-BiDAG::partitionMCMC(myScore, stepsave = 1)
  
  result <- modelpcore(partfit$traceadd$incidence, 0.5)
  
  return(as.matrix(result))
}

# Summary of results given a list of experimental results
get_result <- function(results,truemat){
  com <- lapply(1:nsim, function(i) compareDAGs(results[[i]],truemat,cpdag = TRUE))
  results <- sapply(1:8, function (x) mean(as.numeric(lapply(1:nsim, function(i) com[[i]][x])),na.rm = TRUE)) # Compare to True DAGs
  names(results) <- names(com[[1]])
  return(results)
}

get_metric <- function(results,truemat,metric){
  com <- lapply(1:nsim, function(i) compareDAGs(results[[i]],truemat,cpdag = TRUE))
  results <- as.numeric(lapply(1:nsim, function(i) com[[i]][metric])) # Compare to True DAGs
  return(results)
}

MCMCres <- function(MCMC,truemat){
  compDags <- lapply(1:nsim, function(i) compareDAGs(MCMC[[i]][[2]],truemat,cpdag = TRUE))
  results <- sapply(1:8, function (x) mean(as.numeric(lapply(1:nsim, function(i) compDags[[i]][x])),na.rm = TRUE)) # Compare to True DAGs
  names(results) <- names(compDags[[1]])
  return(results)
}

get_metric_MC <- function(MCMC,truemat,metric){
  com <- lapply(1:nsim, function(i) compareDAGs(MCMC[[i]][[2]],truemat,cpdag = TRUE))
  results <- as.numeric(lapply(1:nsim, function(i) com[[i]][metric])) # Compare to True DAGs
  return(results)
}
