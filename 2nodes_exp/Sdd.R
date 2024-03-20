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
library(parallel)
library(doParallel)

source('2nodes_exp/Utils.R')

set.seed(2024)

nrep = 200
nvar = 2
nsim = 100
iter = 50000

# True DAG structure
model <- model2network("[A][B|A]")
graphviz.plot(model)
truemat <- as.matrix(get.adjacency(graph.edgelist(arcs(model))))

# Empty dataset
dat <- matrix(NA,nrow=nrep,ncol = nvar)
dat <- as_tibble(dat)
colnames(dat) <- LETTERS[1:nvar]

#####S_dd###########################################
#####DISC######
# List of datasets to be used for inference
Sdd_DISC_0.1_ds <- lapply(1:nsim, function(i) Sdd_Data(data=dat,beta=0.1,nrep))
Sdd_DISC_0.25_ds <- lapply(1:nsim, function(i) Sdd_Data(data=dat,beta=0.25,nrep))
Sdd_DISC_0.4_ds <- lapply(1:nsim, function(i) Sdd_Data(data=dat,beta=0.4,nrep))
Sdd_DISC_0.6_ds <- lapply(1:nsim, function(i) Sdd_Data(data=dat,beta=0.6,nrep))
Sdd_DISC_0.9_ds <- lapply(1:nsim, function(i) Sdd_Data(data=dat,beta=0.9,nrep))

# Running inference
numCores <- detectCores()
cl <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl) 

Sdd_DISC_0.1 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdd_DISC_0.1_ds[[i]])
  result$result <- runPartMCMC("bdecat", data, iter)
  return(result)
}

Sdd_DISC_0.25 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdd_DISC_0.25_ds[[i]])
  result$result <- runPartMCMC("bdecat", data, iter)
  return(result)
}

Sdd_DISC_0.4 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdd_DISC_0.4_ds[[i]])
  result$result <- runPartMCMC("bdecat", data, iter)
  return(result)
}

Sdd_DISC_0.6 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdd_DISC_0.6_ds[[i]])
  result$result <- runPartMCMC("bdecat", data, iter)
  return(result)
}

Sdd_DISC_0.9 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdd_DISC_0.9_ds[[i]])
  result$result <- runPartMCMC("bdecat", data, iter)
  return(result)
}

stopCluster(cl)

post_process(Sdd_DISC_0.1)
post_process(Sdd_DISC_0.25)
post_process(Sdd_DISC_0.4)
post_process(Sdd_DISC_0.6)
post_process(Sdd_DISC_0.9)

#####RAG######
# List of datasets to be used for inference
Sdd_RAG_0.1_ds <- lapply(1:nsim, function(i) Sdd_DISC_0.1_ds[[i]] %>% mutate_if(sapply(Sdd_DISC_0.1_ds[[i]], is.factor), as.numeric))
Sdd_RAG_0.25_ds <- lapply(1:nsim, function(i) Sdd_DISC_0.25_ds[[i]] %>% mutate_if(sapply(Sdd_DISC_0.25_ds[[i]], is.factor), as.numeric))
Sdd_RAG_0.4_ds <- lapply(1:nsim, function(i) Sdd_DISC_0.4_ds[[i]] %>% mutate_if(sapply(Sdd_DISC_0.4_ds[[i]], is.factor), as.numeric))
Sdd_RAG_0.6_ds <- lapply(1:nsim, function(i) Sdd_DISC_0.6_ds[[i]] %>% mutate_if(sapply(Sdd_DISC_0.6_ds[[i]], is.factor), as.numeric))
Sdd_RAG_0.9_ds <- lapply(1:nsim, function(i) Sdd_DISC_0.9_ds[[i]] %>% mutate_if(sapply(Sdd_DISC_0.9_ds[[i]], is.factor), as.numeric))

# Running inference
cl <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl) 

Sdd_RAG_0.1 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdd_RAG_0.1_ds[[i]])
  result$result <- runPartMCMC("bge", data, iter)
  return(result)
}

Sdd_RAG_0.25 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdd_RAG_0.25_ds[[i]])
  result$result <- runPartMCMC("bge", data, iter)
  return(result)
}

Sdd_RAG_0.4 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdd_RAG_0.4_ds[[i]])
  result$result <- runPartMCMC("bge", data, iter)
  return(result)
}

Sdd_RAG_0.6 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdd_RAG_0.6_ds[[i]])
  result$result <- runPartMCMC("bge", data, iter)
  return(result)
}

Sdd_RAG_0.9 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdd_RAG_0.9_ds[[i]])
  result$result <- runPartMCMC("bge", data, iter)
  return(result)
}

stopCluster(cl)

post_process(Sdd_RAG_0.1)
post_process(Sdd_RAG_0.25)
post_process(Sdd_RAG_0.4)
post_process(Sdd_RAG_0.6)
post_process(Sdd_RAG_0.9)

#save.image(file="2nodes_exp/Sdd.RData")