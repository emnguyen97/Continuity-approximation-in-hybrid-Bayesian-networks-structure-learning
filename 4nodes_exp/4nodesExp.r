# Libraries

library(BiDAG)
library(bnlearn)
library(tibble)
library(ggdag)
library(dagitty)
library(bnlearn)
library(reshape2)
library(igraph)
library(gRain)
library(dplyr)
library(ggplot2)

# Load the necessary functions

source('4nodes_exp/ContData.R')
source('4nodes_exp/CombData.R')
source('4nodes_exp/DiscreteData.R')
source('4nodes_exp/Inference.R')
source('4nodes_exp/ContDiscretised.r')
source('4nodes_exp/CombDiscretised.r')
source('4nodes_exp/partMCMC.R')

set.seed(3636)

nrep = 200 #No of observations
nvar = 4 #No of variables
nsim = 100 # To be changed to 100

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
datasets <- lapply(1:nsim, function(i) ContData(data=dat,nrep,nvar))
Comb.dats <- lapply(1:nsim, function(i) CombData(data=dat,nrep,nvar))
Comb_cd.dats <- lapply(1:nsim, 
                       function(i) CombData_cd(data=dat,nrep,nvar,muA=-3,muC=6))
dcrt.dats <- lapply(1:nsim, function(i) DiscreteData(data=dat,nrep))



# Running all combinations

type <- c("Continuous", "Combination")

discretise <- c("2l", "4l", "2eq", "4eq")

model <- c("bge", "bdecat")


for (d in discretise){
  for (m in model) {
    result <- partMCMC(type, d, m, nvar, nsim) 
    save.image(paste0("4nodes_exp/",type,d,m,".RData"))
  }
}

for (d in discretise){
  for (m in model) {
    result <- partMCMC(type, d, m, nvar, nsim) 
    save.image(paste0("4nodes_exp/",type,d,m,".RData"))
  }
}


Dcrt <- c("None")
Cat <- c("Continuous", "Combination", "Categorical")

for (c in Cat) {
  for (d in Dcrt){
    result <- partMCMC(c, d, "bge", nvar, nsim) 
    save.image(paste0("4nodes_exp/",c,d,"bge",".RData"))
  }
}


Catbde <- partMCMC("Categorical", "None", "bdecat", nvar, nsim) 
save(Catbde, file = paste0("4nodes_exp/","Categorical","None","bdecat",".RData"))

Catbge <- partMCMC("Categorical", "None", "bge", nvar, nsim) 
save(Catbge, file = paste0("4nodes_exp/","Categorical","None","bge",".RData"))


# Code to compare obtained dags with true dag
l <- length(result)
options(digits = 4)
# Compare to True DAGs
sapply(1:8, function (x) mean(as.numeric(lapply((nsim+1):l, function(i) result[[i]][x])))) 
options(digits = 5)
mean(as.numeric(lapply(1:nsim, function (x) nrow(result[[x]][[1]])))) # Avg no. of unique DAGs
min(as.numeric(lapply(1:nsim, function (x) nrow(result[[x]][[1]]))))
max(as.numeric(lapply(1:nsim, function (x) nrow(result[[x]][[1]]))))


# Extracting frequency table for unique structures
Freq_tab <- lapply(1:nsim, function(i) result[[i]][[1]])

compDAGs <- lapply(1:nsim, function(i) compareDAGs(matrix(Freq_tab[[i]]$StructureStr[[which.max(Freq_tab[[i]]$DAGscores)]],nrow = 4),datmat))
l <- length(compDAGs)

sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs[[i]][x])))) # Compare to True DAGs

# Getting frequency ratios
f_ratio <- list()

getIdx <- function(n, freqtable, k){
  which(sapply(1:n, function(i) identical(freqtable[[k]]$StructureStr[[i]], as.vector(datmat))))
}

fratio <- function(freqtable, k, idx){
  freqtable[[k]]$Frequency[idx]/(100000-freqtable[[k]]$Frequency[idx])
}


for (k in 1:nsim) {
  f_ratio[[k]] <- fratio(Freq_tab, k, getIdx(nrow(Freq_tab[[k]]), Freq_tab, k))
}

mean(unlist(f_ratio))



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
  
  sum(c_freq_unlist)/sum(100000 - c_freq_unlist)
}

# Extracting frequency table for unique structures
Freq_tab <- lapply(1:nsim, function(i) result[[i]][[1]])
F_ratio(Freq_tab)


#####CLG########

runMCMC <- function(i){
  
  data <- as.data.frame(Comb.dats[[i]])
  
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
  
  blklist <- matrix(0,4,4)
  blklist[2,1] <- blklist[2,3] <- blklist[4,1] <- blklist[4,3] <- 1
  colnames(blklist) <- rownames(blklist) <- colnames(data)
  
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

strcMCMC <- foreach(i = seq_len(100), .combine = 'c') %dopar% {
  result <- multiResultClass()
  result$result <- runMCMC(i)
  return(result)
}
stopCluster(cl)
save(strcMCMC, file = "/Users/emmanguyen/Documents/BNProject/4nodes/Sdc.RData")

