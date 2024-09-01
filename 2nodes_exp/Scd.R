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
source('structureMCMC/combinations.R')
source('structureMCMC/scoretables.R')
source('structureMCMC/structurefns.R')
source('structureMCMC/samplefns.R')
source('structureMCMC/param_utils.R')
source('structureMCMC/structureMCMC.R')
source('2nodes_exp/Utils.R')

set.seed(2024)

nrep = 200
nvar = 2
nsim = 100
iter = 50000

# True DAG structure
model <- model2network("[A][B|A]")
graphviz.plot(model)
truemat <- as.matrix(as_adjacency_matrix(graph_from_edgelist(arcs(model))))

# Empty dataset
dat <- matrix(NA,nrow=nrep,ncol = nvar)
dat <- as_tibble(dat, .name_repair = "unique")
colnames(dat) <- LETTERS[1:nvar]

#####S_cd###########################################
#####RAG######
# List of datasets to be used for inference
Scd_0.05_ds <- lapply(1:nsim, function(i) Scd_Data(data=dat,muA=-1,b=0.05,nrep) 
                      %>% mutate(A = as.numeric(A), B = as.numeric(B)))
Scd_0.5_ds <- lapply(1:nsim, function(i) Scd_Data(data=dat,muA=-1,b=0.5,nrep)
                     %>% mutate(A = as.numeric(A), B = as.numeric(B)))
Scd_1_ds <- lapply(1:nsim, function(i) Scd_Data(data=dat,muA=-1,b=1,nrep)
                   %>% mutate(A = as.numeric(A), B = as.numeric(B)))
Scd_1.5_ds <- lapply(1:nsim, function(i) Scd_Data(data=dat,muA=-1,b=1.5,nrep)
                     %>% mutate(A = as.numeric(A), B = as.numeric(B)))
Scd_2_ds <- lapply(1:nsim, function(i) Scd_Data(data=dat,muA=-1,b=2,nrep)
                   %>% mutate(A = as.numeric(A), B = as.numeric(B)))

# Running inference
numCores <- detectCores()
cl <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl) 

Scd_RAG_0.05 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Scd_0.05_ds[[i]])
  result$result <- runPartMCMC("bge", data, iter)
  return(result)
}

Scd_RAG_0.5 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Scd_0.5_ds[[i]])
  result$result <- runPartMCMC("bge", data, iter)
  return(result)
}

Scd_RAG_1 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Scd_1_ds[[i]])
  result$result <- runPartMCMC("bge", data, iter)
  return(result)
}

Scd_RAG_1.5 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Scd_1.5_ds[[i]])
  result$result <- runPartMCMC("bge", data, iter)
  return(result)
}

Scd_RAG_2 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Scd_2_ds[[i]])
  result$result <- runPartMCMC("bge", data, iter)
  return(result)
}

stopCluster(cl)

post_process(Scd_RAG_0.05)
post_process(Scd_RAG_0.5)
post_process(Scd_RAG_1)
post_process(Scd_RAG_1.5)
post_process(Scd_RAG_2)

#####DISC-2######
# List of datasets to be used for inference
Scd_0.05_dcrt_ds <- lapply(1:nsim, function(i) as.data.frame(Scd_Discretised(Scd_0.05_ds[[i]])))
Scd_0.5_dcrt_ds <- lapply(1:nsim, function(i) as.data.frame(Scd_Discretised(Scd_0.5_ds[[i]])))
Scd_1_dcrt_ds <- lapply(1:nsim, function(i) as.data.frame(Scd_Discretised(Scd_1_ds[[i]])))
Scd_1.5_dcrt_ds <- lapply(1:nsim, function(i) as.data.frame(Scd_Discretised(Scd_1.5_ds[[i]])))
Scd_2_dcrt_ds <- lapply(1:nsim, function(i) as.data.frame(Scd_Discretised(Scd_2_ds[[i]])))

cl <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl) 

Scd_DISC_0.05 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Scd_0.05_dcrt_ds[[i]])
  result$result <- runPartMCMC("bdecat", data, iter)
  return(result)
}

Scd_DISC_0.5 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Scd_0.5_dcrt_ds[[i]])
  result$result <- runPartMCMC("bdecat", data, iter)
  return(result)
}

Scd_DISC_1 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Scd_1_dcrt_ds[[i]])
  result$result <- runPartMCMC("bdecat", data, iter)
  return(result)
}

Scd_DISC_1.5 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Scd_1.5_dcrt_ds[[i]])
  result$result <- runPartMCMC("bdecat", data, iter)
  return(result)
}

Scd_DISC_2 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Scd_2_dcrt_ds[[i]])
  result$result <- runPartMCMC("bdecat", data, iter)
  return(result)
}

stopCluster(cl)

post_process(Scd_DISC_0.05)
post_process(Scd_DISC_0.5)
post_process(Scd_DISC_1)
post_process(Scd_DISC_1.5)
post_process(Scd_DISC_2)

#####CLG######
# List of datasets to be used for inference
Scd_0.05_CLG_ds <- lapply(1:nsim, function(i) as.data.frame(Scd_Discretised(Scd_0.05_ds[[i]], clg=TRUE)))
Scd_0.5_CLG_ds <- lapply(1:nsim, function(i) as.data.frame(Scd_Discretised(Scd_0.5_ds[[i]], clg=TRUE)))
Scd_1_CLG_ds <- lapply(1:nsim, function(i) as.data.frame(Scd_Discretised(Scd_1_ds[[i]], clg=TRUE)))
Scd_1.5_CLG_ds <- lapply(1:nsim, function(i) as.data.frame(Scd_Discretised(Scd_1.5_ds[[i]], clg=TRUE)))
Scd_2_CLG_ds <- lapply(1:nsim, function(i) as.data.frame(Scd_Discretised(Scd_2_ds[[i]], clg=TRUE)))

Scd_blklist <- matrix(0,2,2)
Scd_blklist[1,2] <- 1
colnames(Scd_blklist) <- rownames(Scd_blklist) <- colnames(truemat)

# Running inference
cl <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl) 

Scd_CLG_0.05 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Scd_0.05_CLG_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=Scd_blklist,scoretype="bic-cg",sample_parameters=FALSE)
  return(result)
}

Scd_CLG_0.5 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Scd_0.5_CLG_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=Scd_blklist,scoretype="bic-cg",sample_parameters=FALSE)
  return(result)
}

Scd_CLG_1 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Scd_1_CLG_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=Scd_blklist,scoretype="bic-cg",sample_parameters=FALSE)
  return(result)
}

Scd_CLG_1.5 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Scd_1.5_CLG_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=Scd_blklist,scoretype="bic-cg",sample_parameters=FALSE)
  return(result)
}

Scd_CLG_2 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Scd_2_CLG_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=Scd_blklist,scoretype="bic-cg",sample_parameters=FALSE)
  return(result)
}

stopCluster(cl)

post_process(Scd_CLG_0.05)
post_process(Scd_CLG_0.5)
post_process(Scd_CLG_1)
post_process(Scd_CLG_1.5)
post_process(Scd_CLG_2)

#####RAG-Likelihood#####
set.seed(2024)
# Running inference
cl <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl) 

Scd_RAG_0.05_like <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Scd_0.05_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype="bic-g",sample_parameters=FALSE)
  return(result)
}

Scd_RAG_0.5_like <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Scd_0.5_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype="bic-g",sample_parameters=FALSE)
  return(result)
}

Scd_RAG_1_like <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Scd_1_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype="bic-g",sample_parameters=FALSE)
  return(result)
}

Scd_RAG_1.5_like <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Scd_1.5_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype="bic-g",sample_parameters=FALSE)
  return(result)
}

Scd_RAG_2_like <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Scd_2_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype="bic-g",sample_parameters=FALSE)
  return(result)
}

stopCluster(cl)

post_process(Scd_RAG_0.05_like)
post_process(Scd_RAG_0.5_like)
post_process(Scd_RAG_1_like)
post_process(Scd_RAG_1.5_like)
post_process(Scd_RAG_2_like)

save.image(file="2nodes_exp/Scd.RData")