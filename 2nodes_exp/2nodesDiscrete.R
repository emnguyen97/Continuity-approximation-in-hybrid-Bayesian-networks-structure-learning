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

source('4nodes_exp/Inference.R')

set.seed(3636)

nrep = 200
nvar = 2
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
rbern <- function(n, p){
  sims <- sample(0:1, size = n, replace = TRUE, prob = c(1-p, p))
  return(sims)
}

DiscreteData <- function(data, nrep, beta){
  data$A <- rbern(nrep, 0.1)
  for (i in 1:nrep){
    if(data$A[i] == 0) {
      data$B[i] <- rbern(1, 1-beta/2-0.5)
    } else if(data$A[i] == 1) {
      data$B[i] <- rbern(1, beta/2 + 0.5)
    } 
  }
  data <- data %>% mutate_if(sapply(data, is.numeric), as.factor)
  return(data)
}

# List of datasets to be used for inference
Dcrt0.1_ds <- lapply(1:nsim, function(i) DiscreteData(data=dat,beta=0.1,nrep))
Dcrt0.25_ds <- lapply(1:nsim, function(i) DiscreteData(data=dat,beta=0.25,nrep))
Dcrt0.4_ds <- lapply(1:nsim, function(i) DiscreteData(data=dat,beta=0.4,nrep))
Dcrt0.6_ds <- lapply(1:nsim, function(i) DiscreteData(data=dat,beta=0.6,nrep))
Dcrt0.75_ds <- lapply(1:nsim, function(i) DiscreteData(data=dat,beta=0.75,nrep))
Dcrt0.9_ds <- lapply(1:nsim, function(i) DiscreteData(data=dat,beta=0.9,nrep))


# Partition MCMC
Dcrt0.1_bdecat <- lapply(1:nsim, function(i) Inference("bdecat", Dcrt0.1_ds[[i]]))
Dcrt0.25_bdecat <- lapply(1:nsim, function(i) Inference("bdecat", Dcrt0.25_ds[[i]]))
Dcrt0.4_bdecat <- lapply(1:nsim, function(i) Inference("bdecat", Dcrt0.4_ds[[i]]))
Dcrt0.6_bdecat <- lapply(1:nsim, function(i) Inference("bdecat", Dcrt0.6_ds[[i]]))
Dcrt0.75_bdecat <- lapply(1:nsim, function(i) Inference("bdecat", Dcrt0.75_ds[[i]]))
Dcrt0.9_bdecat <- lapply(1:nsim, function(i) Inference("bdecat", Dcrt0.9_ds[[i]]))


# Extracting frequency table for unique structures
Freq_dcrt0.1 <- lapply(1:nsim, function(i) Dcrt0.1_bdecat[[i]][[1]])
Freq_dcrt0.25 <- lapply(1:nsim, function(i) Dcrt0.25_bdecat[[i]][[1]])
Freq_dcrt0.4 <- lapply(1:nsim, function(i) Dcrt0.4_bdecat[[i]][[1]])
Freq_dcrt0.6 <- lapply(1:nsim, function(i) Dcrt0.6_bdecat[[i]][[1]])
Freq_dcrt0.75 <- lapply(1:nsim, function(i) Dcrt0.75_bdecat[[i]][[1]])
Freq_dcrt0.9 <- lapply(1:nsim, function(i) Dcrt0.9_bdecat[[i]][[1]])


# Compare with true DAG
compDAGs_dcrt0.1 <- lapply(1:nsim, function(i) compareDAGs(matrix(Freq_dcrt0.1[[i]]$StructureStr[[which.max(Freq_dcrt0.1[[i]]$DAGscores)]],2,2),truemat))

compDAGs_dcrt0.25 <- lapply(1:nsim, function(i) compareDAGs(matrix(Freq_dcrt0.25[[i]]$StructureStr[[which.max(Freq_dcrt0.25[[i]]$DAGscores)]],2,2),truemat))

compDAGs_dcrt0.4 <- lapply(1:nsim, function(i) compareDAGs(matrix(Freq_dcrt0.4[[i]]$StructureStr[[which.max(Freq_dcrt0.4[[i]]$DAGscores)]],2,2),truemat))

compDAGs_dcrt0.6 <- lapply(1:nsim, function(i) compareDAGs(matrix(Freq_dcrt0.6[[i]]$StructureStr[[which.max(Freq_dcrt0.6[[i]]$DAGscores)]],2,2),truemat))

compDAGs_dcrt0.75 <- lapply(1:nsim, function(i) compareDAGs(matrix(Freq_dcrt0.75[[i]]$StructureStr[[which.max(Freq_dcrt0.75[[i]]$DAGscores)]],2,2),truemat))

compDAGs_dcrt0.9 <- lapply(1:nsim, function(i) compareDAGs(matrix(Freq_dcrt0.9[[i]]$StructureStr[[which.max(Freq_dcrt0.9[[i]]$DAGscores)]],2,2),truemat))


l <- length(compDAGs_dcrt0.1)
# b = 0.1
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_dcrt0.1[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Dcrt0.1_bdecat[[x]][[1]])))) # Avg no. of unique DAGs
# Number of correctly identified dags
length(which(sapply(1:nsim, function(i) identical(Dcrt0.1_bdecat[[i]][[2]],truemat)) == TRUE))


# b = 0.25
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_dcrt0.25[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Dcrt0.25_bdecat[[x]][[1]])))) # Avg no. of unique DAGs
# Number of correctly identified dags
length(which(sapply(1:nsim, function(i) identical(Dcrt0.25_bdecat[[i]][[2]],truemat)) == TRUE))

# b = 0.4
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_dcrt0.4[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Dcrt0.4_bdecat[[x]][[1]])))) # Avg no. of unique DAGs
# Number of correctly identified dags
length(which(sapply(1:nsim, function(i) identical(Dcrt0.4_bdecat[[i]][[2]],truemat)) == TRUE))

# b = 0.6
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_dcrt0.6[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Dcrt0.6_bdecat[[x]][[1]])))) # Avg no. of unique DAGs
# Number of correctly identified dags
length(which(sapply(1:nsim, function(i) identical(Dcrt0.6_bdecat[[i]][[2]],truemat)) == TRUE))

# b = 0.75
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_dcrt0.75[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Dcrt0.75_bdecat[[x]][[1]])))) # Avg no. of unique DAGs
# Number of correctly identified dags
length(which(sapply(1:nsim, function(i) identical(Dcrt0.75_bdecat[[i]][[2]],truemat)) == TRUE))

# b = 0.9
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_dcrt0.9[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Dcrt0.9_bdecat[[x]][[1]])))) # Avg no. of unique DAGs
# Number of correctly identified dags
length(which(sapply(1:nsim, function(i) identical(Dcrt0.9_bdecat[[i]][[2]],truemat)) == TRUE))


DiscreteData.num <- function(data, nrep, beta){
  data$A <- rbern(nrep, 0.5)
  for (i in 1:nrep){
    if(data$A[i] == 0) {
      data$B[i] <- rbern(1, 1-beta)
    } else if(data$A[i] == 1) {
      data$B[i] <- rbern(1, beta)
    } 
  }
  return(data)
}

# List of datasets to be used for inference
Dcrt0.1_ds.n <- lapply(1:nsim, function(i) DiscreteData.num(data=dat,beta=0.1,nrep))
Dcrt0.25_ds.n <- lapply(1:nsim, function(i) DiscreteData.num(data=dat,beta=0.25,nrep))
Dcrt0.4_ds.n <- lapply(1:nsim, function(i) DiscreteData.num(data=dat,beta=0.4,nrep))
Dcrt0.6_ds.n <- lapply(1:nsim, function(i) DiscreteData.num(data=dat,beta=0.6,nrep))
Dcrt0.75_ds.n <- lapply(1:nsim, function(i) DiscreteData.num(data=dat,beta=0.75,nrep))
Dcrt0.9_ds.n <- lapply(1:nsim, function(i) DiscreteData.num(data=dat,beta=0.9,nrep))

# Partition MCMC
Dcrt0.1_bge <- lapply(1:nsim, function(i) Inference("bge", Dcrt0.1_ds.n[[i]]))
Dcrt0.25_bge <- lapply(1:nsim, function(i) Inference("bge", Dcrt0.25_ds.n[[i]]))
Dcrt0.4_bge <- lapply(1:nsim, function(i) Inference("bge", Dcrt0.4_ds.n[[i]]))
Dcrt0.6_bge <- lapply(1:nsim, function(i) Inference("bge", Dcrt0.6_ds.n[[i]]))
Dcrt0.75_bge <- lapply(1:nsim, function(i) Inference("bge", Dcrt0.75_ds.n[[i]]))
Dcrt0.9_bge <- lapply(1:nsim, function(i) Inference("bge", Dcrt0.9_ds.n[[i]]))

# Extracting frequency table for unique structures
Freq_dcrt0.1_bge <- lapply(1:nsim, function(i) Dcrt0.1_bge[[i]][[1]])
Freq_dcrt0.25_bge <- lapply(1:nsim, function(i) Dcrt0.25_bge[[i]][[1]])
Freq_dcrt0.4_bge <- lapply(1:nsim, function(i) Dcrt0.4_bge[[i]][[1]])
Freq_dcrt0.6_bge <- lapply(1:nsim, function(i) Dcrt0.6_bge[[i]][[1]])
Freq_dcrt0.75_bge <- lapply(1:nsim, function(i) Dcrt0.75_bge[[i]][[1]])
Freq_dcrt0.9_bge <- lapply(1:nsim, function(i) Dcrt0.9_bge[[i]][[1]])


# Compare with true DAG
compDAGs_dcrt0.1_bge <- lapply(1:nsim, function(i) compareDAGs(matrix(Freq_dcrt0.1_bge[[i]]$StructureStr[[which.max(Freq_dcrt0.1_bge[[i]]$DAGscores)]],2,2),truemat))

compDAGs_dcrt0.25_bge <- lapply(1:nsim, function(i) compareDAGs(matrix(Freq_dcrt0.25_bge[[i]]$StructureStr[[which.max(Freq_dcrt0.25_bge[[i]]$DAGscores)]],2,2),truemat))

compDAGs_dcrt0.4_bge <- lapply(1:nsim, function(i) compareDAGs(matrix(Freq_dcrt0.4_bge[[i]]$StructureStr[[which.max(Freq_dcrt0.4_bge[[i]]$DAGscores)]],2,2),truemat))

compDAGs_dcrt0.6_bge <- lapply(1:nsim, function(i) compareDAGs(matrix(Freq_dcrt0.6_bge[[i]]$StructureStr[[which.max(Freq_dcrt0.6_bge[[i]]$DAGscores)]],2,2),truemat))

compDAGs_dcrt0.75_bge <- lapply(1:nsim, function(i) compareDAGs(matrix(Freq_dcrt0.75_bge[[i]]$StructureStr[[which.max(Freq_dcrt0.75_bge[[i]]$DAGscores)]],2,2),truemat))

compDAGs_dcrt0.9_bge <- lapply(1:nsim, function(i) compareDAGs(matrix(Freq_dcrt0.9_bge[[i]]$StructureStr[[which.max(Freq_dcrt0.9_bge[[i]]$DAGscores)]],2,2),truemat))

l <- length(compDAGs_dcrt0.1_bge)
# b = 0.1
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_dcrt0.1_bge[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Dcrt0.1_bge[[x]][[1]])))) # Avg no. of unique DAGs
# Number of correctly identified dags
length(which(sapply(1:nsim, function(i) identical(Dcrt0.1_bge[[i]][[2]],truemat)) == TRUE))


# b = 0.25
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_dcrt0.25_bge[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Dcrt0.25_bge[[x]][[1]])))) # Avg no. of unique DAGs
# Number of correctly identified dags
length(which(sapply(1:nsim, function(i) identical(Dcrt0.25_bge[[i]][[2]],truemat)) == TRUE))

# b = 0.4
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_dcrt0.4_bge[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Dcrt0.4_bge[[x]][[1]])))) # Avg no. of unique DAGs
# Number of correctly identified dags
length(which(sapply(1:nsim, function(i) identical(Dcrt0.4_bge[[i]][[2]],truemat)) == TRUE))

# b = 0.6
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_dcrt0.6_bge[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Dcrt0.6_bge[[x]][[1]])))) # Avg no. of unique DAGs
# Number of correctly identified dags
length(which(sapply(1:nsim, function(i) identical(Dcrt0.6_bge[[i]][[2]],truemat)) == TRUE))

# b = 0.75
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_dcrt0.75_bge[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Dcrt0.75_bge[[x]][[1]])))) # Avg no. of unique DAGs
# Number of correctly identified dags
length(which(sapply(1:nsim, function(i) identical(Dcrt0.75_bge[[i]][[2]],truemat)) == TRUE))

# b = 0.9
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_dcrt0.9_bge[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Dcrt0.9_bge[[x]][[1]])))) # Avg no. of unique DAGs
# Number of correctly identified dags
length(which(sapply(1:nsim, function(i) identical(Dcrt0.9_bge[[i]][[2]],truemat)) == TRUE))

# Find structure with no edge

getIdx <- function(n, freqtable, k){
  which(sapply(1:n, function(i) identical(freqtable[[k]]$StructureStr[[i]], rep(0,4))))
}

# Get the frequency for the incorrect dag
get_inc_freq <- function(freqtable, k, idx){
  freqtable[[k]]$Frequency[idx]
}

F_ratio <- function(freqtable){
  # Get the index of incorrect dags for each of the datasets
  inc_indx <- sapply(1:100, function(x) getIdx(nrow(freqtable[[x]]),freqtable, x))
  
  # Get the frequency of incorrect dags for each of the datasets
  inc_freq <- sapply(1:100, function(x) get_inc_freq(freqtable, x, inc_indx[[x]])) 
  
  inc_freq_unlist <- as.numeric(sapply(inc_freq, function(s) if (length(s) == 0) 0 else paste(s, collapse = " ")))
  
  sum(100000 - inc_freq_unlist)/sum(inc_freq_unlist)
}

# Get the index of incorrect dags for each of the datasets
inc_indx <- sapply(1:100, function(x) getIdx(nrow(Freq_dcrt0.25[[x]]),Freq_dcrt0.25, x))

inc_freq <- sapply(1:100, function(x) get_inc_freq(Freq_dcrt0.25, x, inc_indx[[x]])) 

fre_0.25_inc <- as.numeric(sapply(inc_freq, function(s) if (length(s) == 0) 0 else paste(s, collapse = " ")))

sum(100000 - fre_0.25_inc)/sum(fre_0.25_inc)

F_ratio(Freq_dcrt0.1_bge)
F_ratio(Freq_dcrt0.1)

