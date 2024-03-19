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

source('2nodes/Utils.R')

set.seed(3636)

nrep = 200
nvar = 2
nsim = 100
iter = 50000

#####S_cc###########################################

# True DAG structure
model <- model2network("[A][B|A]")
graphviz.plot(model)
truemat <- as.matrix(get.adjacency(graph.edgelist(arcs(model))))

# Empty dataset
dat <- matrix(NA,nrow=nrep,ncol = nvar)
dat <- as_tibble(dat)
colnames(dat) <- LETTERS[1:nvar]

#####RAG######
# List of datasets to be used for inference
# Continuous
Cont0.05_ds <- lapply(1:nsim, function(i) Scc_Data(data=dat,beta=0.05,nrep,nvar))
Cont0.5_ds <- lapply(1:nsim, function(i) Scc_Data(data=dat,beta=0.5,nrep,nvar))
Cont1_ds <- lapply(1:nsim, function(i) Scc_Data(data=dat,beta=1,nrep,nvar))
Cont1.5_ds <- lapply(1:nsim, function(i) Scc_Data(data=dat,beta=1.5,nrep,nvar))
Cont2_ds <- lapply(1:nsim, function(i) Scc_Data(data=dat,beta=2,nrep,nvar))

# Running inference
numCores <- detectCores()
cl <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl) 

Cont0.05 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Cont0.05_ds[[i]])
  result$result <- runPartMCMC("bge", data, iter)
  return(result)
}

Cont0.5 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Cont0.5_ds[[i]])
  result$result <- runPartMCMC("bge", data, iter)
  return(result)
}

Cont1 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Cont1_ds[[i]])
  result$result <- runPartMCMC("bge", data, iter)
  return(result)
}

Cont1.5 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Cont1.5_ds[[i]])
  result$result <- runPartMCMC("bge", data, iter)
  return(result)
}

Cont2 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Cont2_ds[[i]])
  result$result <- runPartMCMC("bge", data, iter)
  return(result)
}

stopCluster(cl)


post_process(Cont0.05)
post_process(Cont0.5)
post_process(Cont1)
post_process(Cont1.5)
post_process(Cont2)

#####DISC-2######
# List of datasets to be used for inference
DISC0.05_ds <- lapply(1:nsim, function(i) as.data.frame(Scc_Discretised(Cont0.05_ds[[i]])))
DISC0.5_ds <- lapply(1:nsim, function(i) as.data.frame(Scc_Discretised(Cont0.5_ds[[i]])))
DISC1_ds <- lapply(1:nsim, function(i) as.data.frame(Scc_Discretised(Cont1_ds[[i]])))
DISC1.5_ds <- lapply(1:nsim, function(i) as.data.frame(Scc_Discretised(Cont1.5_ds[[i]])))
DISC2_ds <- lapply(1:nsim, function(i) as.data.frame(Scc_Discretised(Cont2_ds[[i]])))

# Running inference
cl <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl) 

DISC0.05 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(DISC0.05_ds[[i]])
  result$result <- runPartMCMC("bdecat", data, iter)
  return(result)
}

DISC0.5 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(DISC0.5_ds[[i]])
  result$result <- runPartMCMC("bdecat", data, iter)
  return(result)
}

DISC1 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(DISC1_ds[[i]])
  result$result <- runPartMCMC("bdecat", data, iter)
  return(result)
}

DISC1.5 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(DISC1.5_ds[[i]])
  result$result <- runPartMCMC("bdecat", data, iter)
  return(result)
}

DISC2 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(DISC2_ds[[i]])
  result$result <- runPartMCMC("bdecat", data, iter)
  return(result)
}

stopCluster(cl)

post_process(DISC0.05)
post_process(DISC0.5)
post_process(DISC1)
post_process(DISC1.5)
post_process(DISC2)

save.image(file="2nodes/Scc.RData")

