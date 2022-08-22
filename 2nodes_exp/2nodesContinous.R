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
model <- model2network("[A][B|A]")
graphviz.plot(model)
datmat <- as.matrix(get.adjacency(graph.edgelist(arcs(model))))

# Empty dataset
dat <- matrix(NA,nrow=nrep,ncol = nvar)
dat <- as_tibble(dat)
colnames(dat) <- LETTERS[1:nvar]

# Function to generate Continuous data
ContData <- function(data, beta, nrep, nvar){
  data$A <- rnorm(nrep,mean=-1, sd=1) #  has no parents
  data$B <- data$A * beta + rnorm(nrep)
  return(data)
}

# Function to discretise continuous data
ContDiscretised <- function(data){
  data <- apply(data, 2, function(x) cut(x, breaks = c(min(x),(min(x)+max(x))/2,max(x)),
                                        labels = 0:1,
                                        include.lowest = TRUE))
  return(data)
} 

## Beta = 0.05 ##################################################################
# List of datasets to be used for inference
# Continuous
Cont0.05_ds <- lapply(1:nsim, function(i) ContData(data=dat,beta=0.05,nrep,nvar))
# Discretise
Dcrt0.05_ds <- lapply(1:nsim, function(i) as.data.frame(ContDiscretised(Cont0.05_ds[[i]])))
Dcrt0.05_ds <- lapply(1:nsim, function(x) apply(Dcrt0.05_ds[[x]], 2, as.numeric))

# Partition MCMC
Cont0.05_bge <- lapply(1:nsim, function(i) Inference("bge", Cont0.05_ds[[i]]))
Dcrt0.05_bde <- lapply(1:nsim, function(i) Inference("bde", Dcrt0.05_ds[[i]]))

# Extracting frequency table for unique structures
Freq_cont0.05 <- lapply(1:nsim, function(i) Cont0.05_bge[[i]][[1]])
Freq_dcrt0.05 <- lapply(1:nsim, function(i) Dcrt0.05_bde[[i]][[1]])

# Compare with true DAG
compDAGs_cont0.05.bge <- lapply(1:nsim, 
                        function(i) compareDAGs(matrix(Freq_cont0.05[[i]]$StructureStr[[which.max(Freq_cont0.05[[i]]$DAGscores)]],2,2),datmat))
compDAGs_dcrt0.05.bde <- lapply(1:nsim, 
                        function(i) compareDAGs(matrix(Freq_dcrt0.05[[i]]$StructureStr[[which.max(Freq_dcrt0.05[[i]]$DAGscores)]],2,2),datmat))

l <- length(compDAGs_cont0.05.bge)
# Cont
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_cont0.05.bge[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Cont0.05_bge[[x]][[1]])))) # Avg no. of unique DAGs
# Number of correctly identified dags
length(which(sapply(1:nsim, function(i) identical(Cont0.05_bge[[i]][[2]],datmat)) == TRUE))
# Discrete
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_dcrt0.05.bde[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Dcrt0.05_bde[[x]][[1]])))) # Avg no. of unique DAGs
length(which(sapply(1:nsim, function(i) identical(Dcrt0.05_bde[[i]][[2]],datmat)) == TRUE))

## Beta = 0.1 ##################################################################
# List of datasets to be used for inference
# Continuous
Cont0.1_ds <- lapply(1:nsim, function(i) ContData(data=dat,beta=0.1,nrep,nvar))
# Discretise
Dcrt0.1_ds <- lapply(1:nsim, function(i) as.data.frame(ContDiscretised(Cont0.1_ds[[i]])))
Dcrt0.1_ds <- lapply(1:nsim, function(x) apply(Dcrt0.1_ds[[x]], 2, as.numeric))

# Partition MCMC
Cont0.1_bge <- lapply(1:nsim, function(i) Inference("bge", Cont0.1_ds[[i]]))
Dcrt0.1_bde <- lapply(1:nsim, function(i) Inference("bde", Dcrt0.1_ds[[i]]))

# Extracting frequency table for unique structures
Freq_cont0.1 <- lapply(1:nsim, function(i) Cont0.1_bge[[i]][[1]])
Freq_dcrt0.1 <- lapply(1:nsim, function(i) Dcrt0.1_bde[[i]][[1]])

# Compare with true DAG
compDAGs_cont0.1.bge <- lapply(1:nsim, 
                        function(i) compareDAGs(matrix(Freq_cont0.1[[i]]$StructureStr[[which.max(Freq_cont0.1[[i]]$DAGscores)]],2,2),datmat))
compDAGs_dcrt0.1.bde <- lapply(1:nsim, 
                        function(i) compareDAGs(matrix(Freq_dcrt0.1[[i]]$StructureStr[[which.max(Freq_dcrt0.1[[i]]$DAGscores)]],2,2),datmat))

# To plot
#sapply(1:nsim, 
#       function(i) plotdiffs(matrix(Cont_bge[[i]]$StructureStr[[which.max(Cont_bge[[i]]$DAGscores)]],2,2),datmat))
plotdiffs(Cont0.1_bge[[1]][[2]], datmat)

l <- length(compDAGs_cont0.1.bge)
# Cont
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_cont0.1.bge[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Cont0.1_bge[[x]][[1]])))) # Avg no. of unique DAGs
# Number of correctly identified dags
length(which(sapply(1:nsim, function(i) identical(Cont0.1_bge[[i]][[2]],datmat)) == TRUE))
# Discrete
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_dcrt0.1.bde[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Dcrt0.1_bde[[x]][[1]])))) # Avg no. of unique DAGs
length(which(sapply(1:nsim, function(i) identical(Dcrt0.1_bde[[i]][[2]],datmat)) == TRUE))

## Beta = 0.4 ##################################################################
# List of datasets to be used for inference
# Continuous
Cont0.4_ds <- lapply(1:nsim, function(i) ContData(data=dat,beta=0.4,nrep,nvar))
# Discretise
Dcrt0.4_ds <- lapply(1:nsim, function(i) as.data.frame(ContDiscretised(Cont0.4_ds[[i]])))
Dcrt0.4_ds <- lapply(1:nsim, function(x) apply(Dcrt0.4_ds[[x]], 2, as.numeric))

# Partition MCMC
Cont0.4_bge <- lapply(1:nsim, function(i) Inference("bge", Cont0.4_ds[[i]]))
Dcrt0.4_bde <- lapply(1:nsim, function(i) Inference("bde", Dcrt0.4_ds[[i]]))

# Extracting frequency table for unique structures
Freq_cont0.4 <- lapply(1:nsim, function(i) Cont0.4_bge[[i]][[1]])
Freq_dcrt0.4 <- lapply(1:nsim, function(i) Dcrt0.4_bde[[i]][[1]])

# Compare with true DAG
compDAGs_cont0.4.bge <- lapply(1:nsim, 
                               function(i) compareDAGs(matrix(Freq_cont0.4[[i]]$StructureStr[[which.max(Freq_cont0.4[[i]]$DAGscores)]],2,2),datmat))
compDAGs_dcrt0.4.bde <- lapply(1:nsim, 
                               function(i) compareDAGs(matrix(Freq_dcrt0.4[[i]]$StructureStr[[which.max(Freq_dcrt0.4[[i]]$DAGscores)]],2,2),datmat))

l <- length(compDAGs_cont0.4.bge)
# Cont
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_cont0.4.bge[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Cont0.4_bge[[x]][[1]])))) # Avg no. of unique DAGs
length(which(sapply(1:nsim, function(i) identical(Cont0.4_bge[[i]][[2]],datmat)) == TRUE))
# Discrete
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_dcrt0.4.bde[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Dcrt0.4_bde[[x]][[1]])))) # Avg no. of unique DAGs
length(which(sapply(1:nsim, function(i) identical(Dcrt0.4_bde[[i]][[2]],datmat)) == TRUE))

## Beta = 0.5 ##################################################################
# List of datasets to be used for inference
# Continuous
Cont0.5_ds <- lapply(1:nsim, function(i) ContData(data=dat,beta=0.5,nrep,nvar))
# Discretise
Dcrt0.5_ds <- lapply(1:nsim, function(i) as.data.frame(ContDiscretised(Cont0.5_ds[[i]])))
Dcrt0.5_ds <- lapply(1:nsim, function(x) apply(Dcrt0.5_ds[[x]], 2, as.numeric))

# Partition MCMC
Cont0.5_bge <- lapply(1:nsim, function(i) Inference("bge", Cont0.5_ds[[i]]))
Dcrt0.5_bde <- lapply(1:nsim, function(i) Inference("bde", Dcrt0.5_ds[[i]]))

# Extracting frequency table for unique structures
Freq_cont0.5 <- lapply(1:nsim, function(i) Cont0.5_bge[[i]][[1]])
Freq_dcrt0.5 <- lapply(1:nsim, function(i) Dcrt0.5_bde[[i]][[1]])

# Compare with true DAG
compDAGs_cont0.5.bge <- lapply(1:nsim, 
                               function(i) compareDAGs(matrix(Freq_cont0.5[[i]]$StructureStr[[which.max(Freq_cont0.5[[i]]$DAGscores)]],2,2),datmat))
compDAGs_dcrt0.5.bde <- lapply(1:nsim, 
                               function(i) compareDAGs(matrix(Freq_dcrt0.5[[i]]$StructureStr[[which.max(Freq_dcrt0.5[[i]]$DAGscores)]],2,2),datmat))

l <- length(compDAGs_cont0.5.bge)
# Cont
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_cont0.5.bge[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Cont0.5_bge[[x]][[1]])))) # Avg no. of unique DAGs
length(which(sapply(1:nsim, function(i) identical(Cont0.5_bge[[i]][[2]],datmat)) == TRUE))
# Discrete
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_dcrt0.5.bde[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Dcrt0.5_bde[[x]][[1]])))) # Avg no. of unique DAGs
length(which(sapply(1:nsim, function(i) identical(Dcrt0.5_bde[[i]][[2]],datmat)) == TRUE))

## Beta = 0.8 ##################################################################
# List of datasets to be used for inference
# Continuous
Cont0.8_ds <- lapply(1:nsim, function(i) ContData(data=dat,beta=0.8,nrep,nvar))
# Discretise
Dcrt0.8_ds <- lapply(1:nsim, function(i) as.data.frame(ContDiscretised(Cont0.8_ds[[i]])))
Dcrt0.8_ds <- lapply(1:nsim, function(x) apply(Dcrt0.8_ds[[x]], 2, as.numeric))

# Partition MCMC
Cont0.8_bge <- lapply(1:nsim, function(i) Inference("bge", Cont0.8_ds[[i]]))
Dcrt0.8_bde <- lapply(1:nsim, function(i) Inference("bde", Dcrt0.8_ds[[i]]))

# Extracting frequency table for unique structures
Freq_cont0.8 <- lapply(1:nsim, function(i) Cont0.8_bge[[i]][[1]])
Freq_dcrt0.8 <- lapply(1:nsim, function(i) Dcrt0.8_bde[[i]][[1]])

# Compare with true DAG
compDAGs_cont0.8.bge <- lapply(1:nsim, 
                               function(i) compareDAGs(matrix(Freq_cont0.8[[i]]$StructureStr[[which.max(Freq_cont0.8[[i]]$DAGscores)]],2,2),datmat))
compDAGs_dcrt0.8.bde <- lapply(1:nsim, 
                               function(i) compareDAGs(matrix(Freq_dcrt0.8[[i]]$StructureStr[[which.max(Freq_dcrt0.8[[i]]$DAGscores)]],2,2),datmat))

l <- length(compDAGs_cont0.8.bge)
# Cont
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_cont0.8.bge[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Cont0.8_bge[[x]][[1]])))) # Avg no. of unique DAGs
length(which(sapply(1:nsim, function(i) identical(Cont0.8_bge[[i]][[2]],datmat)) == TRUE))

# Discrete
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_dcrt0.8.bde[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Dcrt0.8_bde[[x]][[1]])))) # Avg no. of unique DAGs
length(which(sapply(1:nsim, function(i) identical(Dcrt0.8_bde[[i]][[2]],datmat)) == TRUE))


## Beta = 1 ##################################################################
# List of datasets to be used for inference
# Continuous
Cont1_ds <- lapply(1:nsim, function(i) ContData(data=dat,beta=1,nrep,nvar))
# Discretise
Dcrt1_ds <- lapply(1:nsim, function(i) as.data.frame(ContDiscretised(Cont1_ds[[i]])))
Dcrt1_ds <- lapply(1:nsim, function(x) apply(Dcrt1_ds[[x]], 2, as.numeric))

# Partition MCMC
Cont1_bge <- lapply(1:nsim, function(i) Inference("bge", Cont1_ds[[i]]))
Dcrt1_bde <- lapply(1:nsim, function(i) Inference("bde", Dcrt1_ds[[i]]))

# Extracting frequency table for unique structures
Freq_cont1 <- lapply(1:nsim, function(i) Cont1_bge[[i]][[1]])
Freq_dcrt1 <- lapply(1:nsim, function(i) Dcrt1_bde[[i]][[1]])

# Compare with true DAG
compDAGs_cont1.bge <- lapply(1:nsim, 
                      function(i) compareDAGs(matrix(Freq_cont1[[i]]$StructureStr[[which.max(Freq_cont1[[i]]$DAGscores)]],2,2),datmat))
compDAGs_dcrt1.bde <- lapply(1:nsim, 
                      function(i) compareDAGs(matrix(Freq_dcrt1[[i]]$StructureStr[[which.max(Freq_dcrt1[[i]]$DAGscores)]],2,2),datmat))

l <- length(compDAGs_cont1.bge)
# Cont
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_cont1.bge[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Cont1_bge[[x]][[1]])))) # Avg no. of unique DAGs
length(which(sapply(1:nsim, function(i) identical(Cont1_bge[[i]][[2]],datmat)) == TRUE))

# Discrete
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_dcrt1.bde[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Dcrt1_bde[[x]][[1]])))) # Avg no. of unique DAGs
length(which(sapply(1:nsim, function(i) identical(Dcrt1_bde[[i]][[2]],datmat)) == TRUE))


## Beta = 1.5 ##################################################################
# List of datasets to be used for inference
# Continuous
Cont1.5_ds <- lapply(1:nsim, function(i) ContData(data=dat,beta=1.5,nrep,nvar))
# Discretise
Dcrt1.5_ds <- lapply(1:nsim, function(i) as.data.frame(ContDiscretised(Cont1.5_ds[[i]])))
Dcrt1.5_ds <- lapply(1:nsim, function(x) apply(Dcrt1.5_ds[[x]], 2, as.numeric))

# Partition MCMC
Cont1.5_bge <- lapply(1:nsim, function(i) Inference("bge", Cont1.5_ds[[i]]))
Dcrt1.5_bde <- lapply(1:nsim, function(i) Inference("bde", Dcrt1.5_ds[[i]]))

# Extracting frequency table for unique structures
Freq_cont1.5 <- lapply(1:nsim, function(i) Cont1.5_bge[[i]][[1]])
Freq_dcrt1.5 <- lapply(1:nsim, function(i) Dcrt1.5_bde[[i]][[1]])

# Compare with true DAG
compDAGs_cont1.5.bge <- lapply(1:nsim, 
                             function(i) compareDAGs(matrix(Freq_cont1.5[[i]]$StructureStr[[which.max(Freq_cont1.5[[i]]$DAGscores)]],2,2),datmat))
compDAGs_dcrt1.5.bde <- lapply(1:nsim, 
                             function(i) compareDAGs(matrix(Freq_dcrt1.5[[i]]$StructureStr[[which.max(Freq_dcrt1.5[[i]]$DAGscores)]],2,2),datmat))

# Cont
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_cont1.5.bge[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Cont1.5_bge[[x]][[1]])))) # Avg no. of unique DAGs
length(which(sapply(1:nsim, function(i) identical(Cont1.5_bge[[i]][[2]],datmat)) == TRUE))

# Discrete
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_dcrt1.5.bde[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Dcrt1.5_bde[[x]][[1]])))) # Avg no. of unique DAGs
length(which(sapply(1:nsim, function(i) identical(Dcrt1.5_bde[[i]][[2]],datmat)) == TRUE))


## Beta = 2 ##################################################################
# List of datasets to be used for inference
# Continuous
Cont2_ds <- lapply(1:nsim, function(i) ContData(data=dat,beta=2,nrep,nvar))
# Discretise
Dcrt2_ds <- lapply(1:nsim, function(i) as.data.frame(ContDiscretised(Cont2_ds[[i]])))
Dcrt2_ds <- lapply(1:nsim, function(x) apply(Dcrt2_ds[[x]], 2, as.numeric))

# Partition MCMC
Cont2_bge <- lapply(1:nsim, function(i) Inference("bge", Cont2_ds[[i]]))
Dcrt2_bde <- lapply(1:nsim, function(i) Inference("bde", Dcrt2_ds[[i]]))

# Extracting frequency table for unique structures
Freq_cont2 <- lapply(1:nsim, function(i) Cont2_bge[[i]][[1]])
Freq_dcrt2 <- lapply(1:nsim, function(i) Dcrt2_bde[[i]][[1]])

# Compare with true DAG
compDAGs_cont2.bge <- lapply(1:nsim, 
                               function(i) compareDAGs(matrix(Freq_cont2[[i]]$StructureStr[[which.max(Freq_cont2[[i]]$DAGscores)]],2,2),datmat))
compDAGs_dcrt2.bde <- lapply(1:nsim, 
                               function(i) compareDAGs(matrix(Freq_dcrt2[[i]]$StructureStr[[which.max(Freq_dcrt2[[i]]$DAGscores)]],2,2),datmat))

# Cont
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_cont2.bge[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Cont2_bge[[x]][[1]])))) # Avg no. of unique DAGs
length(which(sapply(1:nsim, function(i) identical(Cont2_bge[[i]][[2]],datmat)) == TRUE))

# Discrete
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_dcrt2.bde[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Dcrt2_bde[[x]][[1]])))) # Avg no. of unique DAGs
length(which(sapply(1:nsim, function(i) identical(Dcrt2_bde[[i]][[2]],datmat)) == TRUE))

## Beta = 2.5 ##################################################################
# List of datasets to be used for inference
# Continuous
Cont2.5_ds <- lapply(1:nsim, function(i) ContData(data=dat,beta=2.5,nrep,nvar))
# Discretise
Dcrt2.5_ds <- lapply(1:nsim, function(i) as.data.frame(ContDiscretised(Cont2.5_ds[[i]])))
Dcrt2.5_ds <- lapply(1:nsim, function(x) apply(Dcrt2.5_ds[[x]], 2, as.numeric))

# Partition MCMC
Cont2.5_bge <- lapply(1:nsim, function(i) Inference("bge", Cont2.5_ds[[i]]))
Dcrt2.5_bde <- lapply(1:nsim, function(i) Inference("bde", Dcrt2.5_ds[[i]]))

# Extracting frequency table for unique structures
Freq_cont2.5 <- lapply(1:nsim, function(i) Cont2.5_bge[[i]][[1]])
Freq_dcrt2.5 <- lapply(1:nsim, function(i) Dcrt2.5_bde[[i]][[1]])

# Compare with true DAG
compDAGs_cont2.5.bge <- lapply(1:nsim, 
                      function(i) compareDAGs(matrix(Freq_cont2.5[[i]]$StructureStr[[which.max(Freq_cont2.5[[i]]$DAGscores)]],2,2),datmat))
compDAGs_dcrt2.5.bde <- lapply(1:nsim, 
                      function(i) compareDAGs(matrix(Freq_dcrt2.5[[i]]$StructureStr[[which.max(Freq_dcrt2.5[[i]]$DAGscores)]],2,2),datmat))

# Cont
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_cont2.5.bge[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Cont2.5_bge[[x]][[1]])))) # Avg no. of unique DAGs
length(which(sapply(1:nsim, function(i) identical(Cont2.5_bge[[i]][[2]],datmat)) == TRUE))

# Discrete
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_dcrt2.5.bde[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Dcrt2.5_bde[[x]][[1]])))) # Avg no. of unique DAGs
length(which(sapply(1:nsim, function(i) identical(Dcrt2.5_bde[[i]][[2]],datmat)) == TRUE))

## Beta = 3 ##################################################################
# List of datasets to be used for inference
# Continuous
Cont3_ds <- lapply(1:nsim, function(i) ContData(data=dat,beta=3,nrep,nvar))
# Discretise
Dcrt3_ds <- lapply(1:nsim, function(i) as.data.frame(ContDiscretised(Cont3_ds[[i]])))
Dcrt3_ds <- lapply(1:nsim, function(x) apply(Dcrt3_ds[[x]], 2, as.numeric))

# Partition MCMC
Cont3_bge <- lapply(1:nsim, function(i) Inference("bge", Cont3_ds[[i]]))
Dcrt3_bde <- lapply(1:nsim, function(i) Inference("bde", Dcrt3_ds[[i]]))

# Extracting frequency table for unique structures
Freq_cont3 <- lapply(1:nsim, function(i) Cont3_bge[[i]][[1]])
Freq_dcrt3 <- lapply(1:nsim, function(i) Dcrt3_bde[[i]][[1]])

# Compare with true DAG
compDAGs_cont3.bge <- lapply(1:nsim, 
                      function(i) compareDAGs(matrix(Freq_cont3[[i]]$StructureStr[[which.max(Freq_cont3[[i]]$DAGscores)]],2,2),datmat))
compDAGs_dcrt3.bde <- lapply(1:nsim, 
                      function(i) compareDAGs(matrix(Freq_dcrt3[[i]]$StructureStr[[which.max(Freq_dcrt3[[i]]$DAGscores)]],2,2),datmat))

# Cont
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_cont3.bge[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Cont3_bge[[x]][[1]])))) # Avg no. of unique DAGs
length(which(sapply(1:nsim, function(i) identical(Cont3_bge[[i]][[2]],datmat)) == TRUE))

# Discrete
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_dcrt3.bde[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Dcrt3_bde[[x]][[1]])))) # Avg no. of unique DAGs
length(which(sapply(1:nsim, function(i) identical(Dcrt3_bde[[i]][[2]],datmat)) == TRUE))



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


F_ratio(Freq_cont0.1)
F_ratio(Freq_dcrt0.1)