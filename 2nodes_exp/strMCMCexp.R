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

nrep = 200
n = nvar = 2
nsim = 100

# True DAG structure
newmodel <- model2network("[A][B|A]")
graphviz.plot(newmodel)
truemat <- as.matrix(get.adjacency(graph.edgelist(arcs(newmodel))))

# Empty dataset
dat <- matrix(NA,nrow=nrep,ncol = nvar)
dat <- as_tibble(dat)
colnames(dat) <- LETTERS[1:nvar]

# Function to generate combination data
rbern <- function(n, p, clg=FALSE){
  sims <- sample(0:1, size = n, replace = TRUE, prob = c(1-p, p))
  if (clg == TRUE){
    sims <- sample(c("a","b"), size = n, replace = TRUE, prob = c(1-p, p))
  }
  return(sims)
}

CombData <- function(data, beta, nrep, clg=FALSE){
  data$A <- rbern(nrep, p=0.5) #  has no parents
  data$B <- rnorm(nrep, mean = data$A*beta + 1, sd = 1)
  if (clg == TRUE){
    data$A <-as.factor(data$A)
    levels(data$A) <- c("a","b")
  }
  return(data)
}

# List of datasets to be used for inference
Comb0.05_ds <- lapply(1:nsim, function(i) CombData(data=dat,beta=0.05,nrep,clg=TRUE))
Comb0.1_ds <- lapply(1:nsim, function(i) CombData(data=dat,beta=0.1,nrep,clg=TRUE))
Comb0.5_ds <- lapply(1:nsim, function(i) CombData(data=dat,beta=0.5,nrep,clg=TRUE))
Comb1_ds <- lapply(1:nsim, function(i) CombData(data=dat,beta=1,nrep,clg=TRUE))
Comb1.5_ds <- lapply(1:nsim, function(i) CombData(data=dat,beta=1.5,nrep,clg=TRUE))
Comb2_ds <- lapply(1:nsim, function(i) CombData(data=dat,beta=2,nrep,clg=TRUE))
Comb2.5_ds <- lapply(1:nsim, function(i) CombData(data=dat,beta=2.5,nrep,clg=TRUE))
Comb3_ds <- lapply(1:nsim, function(i) CombData(data=dat,beta=3,nrep,clg=TRUE))
Comb5_ds <- lapply(1:nsim, function(i) CombData(data=dat,beta=5,nrep,clg=TRUE))
Comb10_ds <- lapply(1:nsim, function(i) CombData(data=dat,beta=10,nrep,clg=TRUE))


#####S_dc########

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
  colnames(maxDAG) <- rownames(maxDAG) <- colnames(truemat)
  
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

blklist <- matrix(0,2,2)
blklist[2,1] <- 1
colnames(blklist) <- rownames(blklist) <- colnames(truemat)

Comb0.05 <- foreach(i = seq_len(100), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Comb0.05_ds[[i]])
  result$result <- runMCMC(data,blklist)
  return(result)
}

Comb0.1 <- foreach(i = seq_len(100), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Comb0.1_ds[[i]])
  result$result <- runMCMC(data,blklist)
  return(result)
}

Comb0.5 <- foreach(i = seq_len(100), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Comb0.5_ds[[i]])
  result$result <- runMCMC(data,blklist)
  return(result)
}

Comb1 <- foreach(i = seq_len(100), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Comb1_ds[[i]])
  result$result <- runMCMC(data,blklist)
  return(result)
}

Comb1.5 <- foreach(i = seq_len(100), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Comb1.5_ds[[i]])
  result$result <- runMCMC(data,blklist)
  return(result)
}

Comb2 <- foreach(i = seq_len(100), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Comb2_ds[[i]])
  result$result <- runMCMC(data,blklist)
  return(result)
}

stopCluster(cl)
save.image(file="RAG/2nodes_exp/2nodes-Sdc.RData")

#########################################################
# Extracting frequency table for unique structures
Freqtab <- lapply(1:nsim, function(i) Comb2[[i]][[1]])

compDags <- lapply(1:nsim, function(i) compareDAGs(Comb2[[i]][[3]],truemat))

options(digits = 4)

sapply(1:8, function (x) mean(as.numeric(lapply(1:nsim, function(i) compDags[[i]][x])),na.rm = TRUE)) # Compare to True DAGs

# Getting frequency ratios

getIdx <- function(n, freqtable, k){
  which(sapply(1:n, function(i) identical(freqtable[[k]]$StructureStr[[i]], rep(0,4))))
}

# Get the frequency for the correct dag
get_freq <- function(freqtable, k, idx){
  freqtable[[k]]$Frequency[idx]
}

F_ratio <- function(freqtable){
  # Get the index of incorrect dags for each of the datasets
  c_indx <- sapply(1:nsim, function(x) getIdx(nrow(freqtable[[x]]),freqtable, x))
  
  # Get the frequency of incorrect dags for each of the datasets
  c_freq <- sapply(1:nsim, function(x) get_freq(freqtable, x, c_indx[[x]])) 
  
  c_freq_unlist <- as.numeric(sapply(c_freq, function(s) if (length(s) == 0) 0 else paste(s, collapse = " ")))
  
  sum(c_freq_unlist)/sum(10000 - c_freq_unlist)
}

# Extracting frequency table for unique structures
F_ratio(Freqtab)

