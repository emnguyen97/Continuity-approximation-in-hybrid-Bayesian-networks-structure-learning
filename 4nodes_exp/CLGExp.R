# Libraries

library(BiDAG)
library(bnlearn)
library(tibble)
#library(ggdag)
library(dagitty)
library(bnlearn)
library(reshape2)
library(igraph)
library(gRain)
library(dplyr)
library(ggplot2)

# Load the necessary functions
source('RAG/structureMCMC/combinations.R')
source('RAG/structureMCMC/scoretables.R')
source('RAG/structureMCMC/structurefns.R')
source('RAG/structureMCMC/samplefns.R')
source('RAG/structureMCMC/structureMCMC.R')
source('RAG/4nodes_exp/CombData.R')
source('RAG/4nodes_exp/DiscreteData.R')

set.seed(3636)

nrep = 200 #No of observations
n <- nvar <- 4 #No of variables
nsim = 100 # no of sims

# Empty dataset

dat <- matrix(NA,nrow=nrep,ncol = nvar)
dat <- as_tibble(dat)
colnames(dat) <- LETTERS[1:nvar]

# True DAG structure
model <- model2network("[A][C][B|A:C][D|A:C]")
#graphviz.plot(model)
datmat <- as.matrix(get.adjacency(graph.edgelist(arcs(model))))
datmat <- as.data.frame(datmat)[,order(colnames(as.data.frame(datmat)))]
datmat <- datmat[order(rownames(datmat)),]
datmat <- as.matrix(datmat)


# List of datasets to be used for inference
Comb.dats <- lapply(1:nsim, function(i) CombData(data=dat,nrep,nvar))
Comb_cd.dats <- lapply(1:nsim, 
                       function(i) CombData_cd(data=dat,nrep,nvar,muA=-3,muC=6))

#####CLG########

runMCMC <- function(data,blklist){
  
  # Choose maximum number of parents
  
  maxparents<-n-1 # Maximum number of parents to allow
  
  # Fill up a matrix with possible parents
  
  parenttable<-listpossibleparents(maxparents,c(1:n))
  tablelength<-nrow(parenttable[[1]]) # size of the table
  
  # Now need to score them!
  
  #scoretable<-scorepossibleparents(parenttable,n) 
  
  iterations<-10000 #number of iterations in the chain
  moveprobs<-1 # having length 1 disallows the new edge reversal move
  #if(!(length(moveprobs)==1)){print('Vector of move probabilities has the wrong length!')}
  
  stepsave<-1 #stepsave<-iterations/1000 #how often to save the result
  
  set.seed(97)
  
  startDAG<-matrix(numeric(n*n),nrow=n) # starting DAG is empty say
  
  revallowed<-1 # allow standard edge reversals
  
  example<-structureMCMC(n,data,startDAG,iterations,stepsave,maxparents,parenttable,scoretable=NULL,revallowed,moveprobs,blacklist = blklist,scoretype = "bic-cg")
  
  uniqstructure <- unique(example$incidence)
  
  freq <- tabulate(match(example$incidence, uniqstructure))
  
  bgefreq <- matrix(NA,nrow=length(uniqstructure),ncol = 4)
  colnames(bgefreq) <- c("Indx", "StructureStr", "DAGscores", "Frequency")
  bgefreq <- as.data.frame(bgefreq)
  bgefreq$Frequency <- freq
  bgefreq$Indx <- match(uniqstructure,example$incidence)
  bgefreq$DAGscores <- sapply(bgefreq$Indx, function(x) example$DAGlogscore[[x]])
  bgefreq$StructureStr <- lapply(bgefreq$Indx, function(x) as.vector(example$incidence[[x]]))
  
  DAGscores<-unlist(example$DAGlogscore)
  maxDAG<-example$incidence[which(DAGscores==max(DAGscores))][[1]]
  colnames(maxDAG) <- rownames(maxDAG) <- colnames(datmat)
  
  return(list(summary=bgefreq, DAGlogscores=DAGscores, DAG=maxDAG))
}

#coda::traceplot(as.mcmc(unlist(example$DAGlogscore)[-c(1:100)]))

library(parallel)
library(doParallel)

multiResultClass <- function(result=NULL)
{
  me <- list(
    result = result
  )
  
  ## Set the name for the class
  class(me) <- append(class(me),"multiResultClass")
  return(me)
}
strcMCMC <- list()
#cl <- parallel::makeCluster(4)
#doParallel::registerDoParallel(cl)
numCores <- detectCores()
cl <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl) 

blklist <- matrix(0,4,4)
blklist[2,1] <- blklist[2,3] <- blklist[4,1] <- blklist[4,3] <- 1
colnames(blklist) <- rownames(blklist) <- colnames(datmat)

strcMCMC <- foreach(i = seq_len(100), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Comb.dats[[i]])
  result$result <- runMCMC(data,blklist)
  return(result)
}
stopCluster(cl)
save(strcMCMC, file = "4nodes_exp/4nodes-Sdc.RData")

############### Scd #####################################

multiResultClass <- function(result=NULL)
{
  me <- list(
    result = result
  )
  
  ## Set the name for the class
  class(me) <- append(class(me),"multiResultClass")
  return(me)
}
strcMCMC <- list()
#cl <- parallel::makeCluster(4)
#doParallel::registerDoParallel(cl)
numCores <- detectCores()
cl <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl) 

blklist_cd <- matrix(0,4,4)
blklist_cd[1,2] <- blklist_cd[3,2] <- blklist_cd[1,4] <- blklist_cd[3,4] <- 1
colnames(blklist_cd) <- rownames(blklist_cd) <- colnames(datmat)

strcMCMC <- foreach(i = seq_len(100), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Comb_cd.dats[[i]])
  result$result <- runMCMC(data,blklist_cd)
  return(result)
}
stopCluster(cl)
save(strcMCMC, file = "4nodes/4nodes-Scd.RData")


#########################################################
# Extracting frequency table for unique structures
Freqtab <- lapply(1:nsim, function(i) strcMCMC[[i]][[1]])

compDags <- lapply(1:nsim, function(i) compareDAGs(strcMCMC[[i]][[3]],datmat))

options(digits = 4)

sapply(1:8, function (x) mean(as.numeric(lapply(1:nsim, function(i) compDags[[i]][x])),na.rm = TRUE)) # Compare to True DAGs

# Getting frequency ratios

getIdx <- function(n, freqtable, k){
  which(sapply(1:n, function(i) identical(freqtable[[k]]$StructureStr[[i]], as.vector(datmat))))
}

# Get the frequency for the correct dag
get_freq <- function(freqtable, k, idx){
  freqtable[[k]]$Frequency[idx]
}

F_ratio <- function(freqtable){
  # Get the index of correct dags for each of the datasets
  c_indx <- sapply(1:nsim, function(x) getIdx(nrow(freqtable[[x]]),freqtable, x))
  
  # Get the frequency of correct dags for each of the datasets
  c_freq <- sapply(1:nsim, function(x) get_freq(freqtable, x, c_indx[[x]])) 
  
  c_freq_unlist <- as.numeric(sapply(c_freq, function(s) if (length(s) == 0) 0 else paste(s, collapse = " ")))
  
  sum(c_freq_unlist)/sum(10000 - c_freq_unlist)
}

# Extracting frequency table for unique structures
F_ratio(Freqtab)

