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
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype=NULL,sample_parameters=TRUE)
  return(result)
}

Cont0.5 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Cont0.5_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype=NULL,sample_parameters=TRUE)
  return(result)
}

Cont1 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Cont1_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype=NULL,sample_parameters=TRUE)
  return(result)
}

Cont1.5 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Cont1.5_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype=NULL,sample_parameters=TRUE)
  return(result)
}

Cont2 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Cont2_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype=NULL,sample_parameters=TRUE)
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
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype="bde",sample_parameters=FALSE)
  return(result)
}

DISC0.5 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(DISC0.5_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype="bde",sample_parameters=FALSE)
  return(result)
}

DISC1 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(DISC1_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype="bde",sample_parameters=FALSE)
  return(result)
}

DISC1.5 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(DISC1.5_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype="bde",sample_parameters=FALSE)
  return(result)
}

DISC2 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(DISC2_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype="bde",sample_parameters=FALSE)
  return(result)
}

stopCluster(cl)

post_process(DISC0.05)
post_process(DISC0.5)
post_process(DISC1)
post_process(DISC1.5)
post_process(DISC2)


#####S_cd###########################################
#####RAG######
# List of datasets to be used for inference
Scd_0.05_ds <- lapply(1:nsim, function(i) Scd_Data(data=dat,muA=-1,b=0.05,nrep))
Scd_0.5_ds <- lapply(1:nsim, function(i) Scd_Data(data=dat,muA=-1,b=0.5,nrep))
Scd_1_ds <- lapply(1:nsim, function(i) Scd_Data(data=dat,muA=-1,b=1,nrep))
Scd_1.5_ds <- lapply(1:nsim, function(i) Scd_Data(data=dat,muA=-1,b=1.5,nrep))
Scd_2_ds <- lapply(1:nsim, function(i) Scd_Data(data=dat,muA=-1,b=2,nrep))

# Running inference
cl <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl) 

Scd_RAG_0.05 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Scd_0.05_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype=NULL,sample_parameters=TRUE)
  return(result)
}

Scd_RAG_0.5 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Scd_0.5_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype=NULL,sample_parameters=TRUE)
  return(result)
}

Scd_RAG_1 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Scd_1_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype=NULL,sample_parameters=TRUE)
  return(result)
}

Scd_RAG_1.5 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Scd_1.5_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype=NULL,sample_parameters=TRUE)
  return(result)
}

Scd_RAG_2 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Scd_2_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype=NULL,sample_parameters=TRUE)
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
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype="bde",sample_parameters=FALSE)
  return(result)
}

Scd_DISC_0.5 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Scd_0.5_dcrt_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype="bde",sample_parameters=FALSE)
  return(result)
}

Scd_DISC_1 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Scd_1_dcrt_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype="bde",sample_parameters=FALSE)
  return(result)
}

Scd_DISC_1.5 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Scd_1.5_dcrt_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype="bde",sample_parameters=FALSE)
  return(result)
}

Scd_DISC_2 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Scd_2_dcrt_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype="bde",sample_parameters=FALSE)
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

#####S_dc###########################################
#####RAG######
# List of datasets to be used for inference
Sdc_0.05_ds <- lapply(1:nsim, function(i) Sdc_Data(data=dat,beta=0.05,nrep))
Sdc_0.5_ds <- lapply(1:nsim, function(i) Sdc_Data(data=dat,beta=0.5,nrep))
Sdc_1_ds <- lapply(1:nsim, function(i) Sdc_Data(data=dat,beta=1,nrep))
Sdc_1.5_ds <- lapply(1:nsim, function(i) Sdc_Data(data=dat,beta=1.5,nrep))
Sdc_2_ds <- lapply(1:nsim, function(i) Sdc_Data(data=dat,beta=2,nrep))

# Running inference
cl <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl) 

Sdc_RAG_0.05 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdc_0.05_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype=NULL,sample_parameters=TRUE)
  return(result)
}

Sdc_RAG_0.5 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdc_0.5_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype=NULL,sample_parameters=TRUE)
  return(result)
}

Sdc_RAG_1 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdc_1_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype=NULL,sample_parameters=TRUE)
  return(result)
}

Sdc_RAG_1.5 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdc_1.5_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype=NULL,sample_parameters=TRUE)
  return(result)
}

Sdc_RAG_2 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdc_2_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype=NULL,sample_parameters=TRUE)
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
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype="bde",sample_parameters=FALSE)
  return(result)
}

Sdc_DISC_0.5 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdc_0.5_dcrt_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype="bde",sample_parameters=FALSE)
  return(result)
}

Sdc_DISC_1 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdc_1_dcrt_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype="bde",sample_parameters=FALSE)
  return(result)
}

Sdc_DISC_1.5 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdc_1.5_dcrt_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype="bde",sample_parameters=FALSE)
  return(result)
}

Sdc_DISC_2 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdc_2_dcrt_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype="bde",sample_parameters=FALSE)
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

#####S_dd###########################################
#####DISC######
# List of datasets to be used for inference
Sdd_DISC_0.1_ds <- lapply(1:nsim, function(i) Sdd_Data(data=dat,beta=0.1,nrep))
Sdd_DISC_0.25_ds <- lapply(1:nsim, function(i) Sdd_Data(data=dat,beta=0.25,nrep))
Sdd_DISC_0.4_ds <- lapply(1:nsim, function(i) Sdd_Data(data=dat,beta=0.4,nrep))
Sdd_DISC_0.6_ds <- lapply(1:nsim, function(i) Sdd_Data(data=dat,beta=0.6,nrep))
Sdd_DISC_0.9_ds <- lapply(1:nsim, function(i) Sdd_Data(data=dat,beta=0.9,nrep))

# Running inference
cl <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl) 

Sdd_DISC_0.1 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdd_DISC_0.1_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype="bde",sample_parameters=FALSE)
  return(result)
}

Sdd_DISC_0.25 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdd_DISC_0.25_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype="bde",sample_parameters=FALSE)
  return(result)
}

Sdd_DISC_0.4 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdd_DISC_0.4_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype="bde",sample_parameters=FALSE)
  return(result)
}

Sdd_DISC_0.6 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdd_DISC_0.6_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype="bde",sample_parameters=FALSE)
  return(result)
}

Sdd_DISC_0.9 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdd_DISC_0.9_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype="bde",sample_parameters=FALSE)
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
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype=NULL,sample_parameters=TRUE)
  return(result)
}

Sdd_RAG_0.25 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdd_RAG_0.25_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype=NULL,sample_parameters=TRUE)
  return(result)
}

Sdd_RAG_0.4 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdd_RAG_0.4_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype=NULL,sample_parameters=TRUE)
  return(result)
}

Sdd_RAG_0.6 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdd_RAG_0.6_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype=NULL,sample_parameters=TRUE)
  return(result)
}

Sdd_RAG_0.9 <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(Sdd_RAG_0.9_ds[[i]])
  result$result <- runStructMCMC(data,iterations = iter,blklist=NULL,scoretype=NULL,sample_parameters=TRUE)
  return(result)
}

stopCluster(cl)

post_process(Sdd_RAG_0.1)
post_process(Sdd_RAG_0.25)
post_process(Sdd_RAG_0.4)
post_process(Sdd_RAG_0.6)
post_process(Sdd_RAG_0.9)


save.image(file="RAG/2nodes_exp/Lemma1.RData")

