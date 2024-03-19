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

# Load the necessary functions
source('StructureMCMC/combinations.R')
source('StructureMCMC/scoretables.R')
source('StructureMCMC/structurefns.R')
source('StructureMCMC/samplefns.R')
source('StructureMCMC/param_utils.R')
source('StructureMCMC/structureMCMC.R')
source('2nodes/Utils.R')

set.seed(3636)

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

#####S_dc###########################################
#####RAG######
# List of datasets to be used for inference
Sdc_0.05_ds <- lapply(1:nsim, function(i) Sdc_Data(data=dat,beta=0.05,nrep))
Sdc_0.5_ds <- lapply(1:nsim, function(i) Sdc_Data(data=dat,beta=0.5,nrep))
Sdc_1_ds <- lapply(1:nsim, function(i) Sdc_Data(data=dat,beta=1,nrep))
Sdc_1.5_ds <- lapply(1:nsim, function(i) Sdc_Data(data=dat,beta=1.5,nrep))
Sdc_2_ds <- lapply(1:nsim, function(i) Sdc_Data(data=dat,beta=2,nrep))

# Running inference
numCores <- detectCores()
cl <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl) 

Sdc_RAG_0.05 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdc_0.05_ds[[i]])
  result$result <- runPartMCMC("bge", data, iter)
  return(result)
}

Sdc_RAG_0.5 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdc_0.5_ds[[i]])
  result$result <- runPartMCMC("bge", data, iter)
  return(result)
}

Sdc_RAG_1 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdc_1_ds[[i]])
  result$result <- runPartMCMC("bge", data, iter)
  return(result)
}

Sdc_RAG_1.5 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdc_1.5_ds[[i]])
  result$result <- runPartMCMC("bge", data, iter)
  return(result)
}

Sdc_RAG_2 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdc_2_ds[[i]])
  result$result <- runPartMCMC("bge", data, iter)
  return(result)
}

stopCluster(cl)

post_process(Sdc_RAG_0.05)
post_process(Sdc_RAG_0.5)
post_process(Sdc_RAG_1)
post_process(Sdc_RAG_1.5)
post_process(Sdc_RAG_2)

#####DISC-2######
# List of datasets to be used for inference
Sdc_0.05_dcrt_ds <- lapply(1:nsim, function(i) as.data.frame(Sdc_Discretised(Sdc_0.05_ds[[i]])))
Sdc_0.5_dcrt_ds <- lapply(1:nsim, function(i) as.data.frame(Sdc_Discretised(Sdc_0.5_ds[[i]])))
Sdc_1_dcrt_ds <- lapply(1:nsim, function(i) as.data.frame(Sdc_Discretised(Sdc_1_ds[[i]])))
Sdc_1.5_dcrt_ds <- lapply(1:nsim, function(i) as.data.frame(Sdc_Discretised(Sdc_1.5_ds[[i]])))
Sdc_2_dcrt_ds <- lapply(1:nsim, function(i) as.data.frame(Sdc_Discretised(Sdc_2_ds[[i]])))

# Running inference
cl <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl) 

Sdc_DISC_0.05 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdc_0.05_dcrt_ds[[i]])
  result$result <- runPartMCMC("bdecat", data, iter)
  return(result)
}

Sdc_DISC_0.5 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdc_0.5_dcrt_ds[[i]])
  result$result <- runPartMCMC("bdecat", data, iter)
  return(result)
}

Sdc_DISC_1 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdc_1_dcrt_ds[[i]])
  result$result <- runPartMCMC("bdecat", data, iter)
  return(result)
}

Sdc_DISC_1.5 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdc_1.5_dcrt_ds[[i]])
  result$result <- runPartMCMC("bdecat", data, iter)
  return(result)
}

Sdc_DISC_2 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdc_2_dcrt_ds[[i]])
  result$result <- runPartMCMC("bdecat", data, iter)
  return(result)
}

stopCluster(cl)

post_process(Sdc_DISC_0.05)
post_process(Sdc_DISC_0.5)
post_process(Sdc_DISC_1)
post_process(Sdc_DISC_1.5)
post_process(Sdc_DISC_2)

#####CLG######
# List of datasets to be used for inference
Sdc_0.05_CLG_ds <- lapply(1:nsim, function(i) as.data.frame(Sdc_Discretised(Sdc_0.05_ds[[i]], clg=TRUE)))
Sdc_0.5_CLG_ds <- lapply(1:nsim, function(i) as.data.frame(Sdc_Discretised(Sdc_0.5_ds[[i]], clg=TRUE)))
Sdc_1_CLG_ds <- lapply(1:nsim, function(i) as.data.frame(Sdc_Discretised(Sdc_1_ds[[i]], clg=TRUE)))
Sdc_1.5_CLG_ds <- lapply(1:nsim, function(i) as.data.frame(Sdc_Discretised(Sdc_1.5_ds[[i]], clg=TRUE)))
Sdc_2_CLG_ds <- lapply(1:nsim, function(i) as.data.frame(Sdc_Discretised(Sdc_2_ds[[i]], clg=TRUE)))

Sdc_blklist <- matrix(0,2,2)
Sdc_blklist[2,1] <- 1
colnames(Sdc_blklist) <- rownames(Sdc_blklist) <- colnames(truemat)

# Running inference
cl <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl) 

Sdc_CLG_0.05 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdc_0.05_CLG_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=Sdc_blklist,scoretype="bic-cg",sample_parameters=FALSE)
  return(result)
}

Sdc_CLG_0.5 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdc_0.5_CLG_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=Sdc_blklist,scoretype="bic-cg",sample_parameters=FALSE)
  return(result)
}

Sdc_CLG_1 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdc_1_CLG_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=Sdc_blklist,scoretype="bic-cg",sample_parameters=FALSE)
  return(result)
}

Sdc_CLG_1.5 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdc_1.5_CLG_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=Sdc_blklist,scoretype="bic-cg",sample_parameters=FALSE)
  return(result)
}

Sdc_CLG_2 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdc_2_CLG_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=Sdc_blklist,scoretype="bic-cg",sample_parameters=FALSE)
  return(result)
}

stopCluster(cl)

post_process(Sdc_CLG_0.05)
post_process(Sdc_CLG_0.5)
post_process(Sdc_CLG_1)
post_process(Sdc_CLG_1.5)
post_process(Sdc_CLG_2)

save.image(file="2nodes/Sdc.RData")