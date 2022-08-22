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


