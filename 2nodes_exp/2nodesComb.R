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

CombData <- function(data, muA, b, nrep){
  data$A <- rnorm(nrep,mean=muA, sd=1) #  has no parents
  # Bern probability
  p <- function(a, muA) exp(b*(a - muA))/(1 + exp(b*(a - muA)))
  for (i in 1:nrep){
    data$B[i] <- rbern(1, p(data$A[i], muA))
  }
  return(data)
}

# List of datasets to be used for inference
Comb0.05_ds <- lapply(1:nsim, function(i) CombData(data=dat,muA=-1,b=0.05,nrep))
Comb0.1_ds <- lapply(1:nsim, function(i) CombData(data=dat,muA=-1,b=0.1,nrep))
Comb0.5_ds <- lapply(1:nsim, function(i) CombData(data=dat,muA=-1,b=0.5,nrep))
Comb1_ds <- lapply(1:nsim, function(i) CombData(data=dat,muA=-1,b=1,nrep))
Comb1.5_ds <- lapply(1:nsim, function(i) CombData(data=dat,muA=-1,b=1.5,nrep))
Comb2_ds <- lapply(1:nsim, function(i) CombData(data=dat,muA=-1,b=2,nrep))
Comb2.5_ds <- lapply(1:nsim, function(i) CombData(data=dat,muA=-1,b=2.5,nrep))
Comb3_ds <- lapply(1:nsim, function(i) CombData(data=dat,muA=-1,b=3,nrep))
Comb5_ds <- lapply(1:nsim, function(i) CombData(data=dat,muA=-1,b=5,nrep))
Comb10_ds <- lapply(1:nsim, function(i) CombData(data=dat,muA=-1,b=10,nrep))

# Partition MCMC
Comb0.05_bge <- lapply(1:nsim, function(i) Inference("bge", Comb0.05_ds[[i]]))
Comb0.1_bge <- lapply(1:nsim, function(i) Inference("bge", Comb0.1_ds[[i]]))
Comb0.5_bge <- lapply(1:nsim, function(i) Inference("bge", Comb0.5_ds[[i]]))
Comb1_bge <- lapply(1:nsim, function(i) Inference("bge", Comb1_ds[[i]]))
Comb1.5_bge <- lapply(1:nsim, function(i) Inference("bge", Comb1.5_ds[[i]]))
Comb2_bge <- lapply(1:nsim, function(i) Inference("bge", Comb2_ds[[i]]))
Comb2.5_bge <- lapply(1:nsim, function(i) Inference("bge", Comb2.5_ds[[i]]))
Comb3_bge <- lapply(1:nsim, function(i) Inference("bge", Comb3_ds[[i]]))
Comb5_bge <- lapply(1:nsim, function(i) Inference("bge", Comb5_ds[[i]]))
Comb10_bge <- lapply(1:nsim, function(i) Inference("bge", Comb10_ds[[i]]))

# Extracting frequency table for unique structures
Freq_comb0.05.bge <- lapply(1:nsim, function(i) Comb0.05_bge[[i]][[1]])
Freq_comb0.1.bge <- lapply(1:nsim, function(i) Comb0.1_bge[[i]][[1]])
Freq_comb0.5.bge <- lapply(1:nsim, function(i) Comb0.5_bge[[i]][[1]])
Freq_comb1.bge <- lapply(1:nsim, function(i) Comb1_bge[[i]][[1]])
Freq_comb1.5.bge <- lapply(1:nsim, function(i) Comb1.5_bge[[i]][[1]])
Freq_comb2.bge <- lapply(1:nsim, function(i) Comb2_bge[[i]][[1]])
Freq_comb2.5.bge <- lapply(1:nsim, function(i) Comb2.5_bge[[i]][[1]])
Freq_comb3.bge <- lapply(1:nsim, function(i) Comb3_bge[[i]][[1]])
Freq_comb5.bge <- lapply(1:nsim, function(i) Comb5_bge[[i]][[1]])
Freq_comb10.bge <- lapply(1:nsim, function(i) Comb10_bge[[i]][[1]])

# Compare with true DAG
compDAGs_comb0.05.bge <- lapply(1:nsim, 
                         function(i) compareDAGs(matrix(Freq_comb0.05.bge[[i]]$StructureStr[[which.max(Freq_comb0.05.bge[[i]]$DAGscores)]],2,2),truemat))

compDAGs_comb0.1.bge <- lapply(1:nsim, 
                        function(i) compareDAGs(matrix(Freq_comb0.1.bge[[i]]$StructureStr[[which.max(Freq_comb0.1.bge[[i]]$DAGscores)]],2,2),truemat))

compDAGs_comb0.5.bge <- lapply(1:nsim, 
                             function(i) compareDAGs(matrix(Freq_comb0.5.bge[[i]]$StructureStr[[which.max(Freq_comb0.5.bge[[i]]$DAGscores)]],2,2),truemat))

compDAGs_comb1.bge <- lapply(1:nsim, 
                             function(i) compareDAGs(matrix(Freq_comb1.bge[[i]]$StructureStr[[which.max(Freq_comb1.bge[[i]]$DAGscores)]],2,2),truemat))

compDAGs_comb1.5.bge <- lapply(1:nsim, 
                               function(i) compareDAGs(matrix(Freq_comb1.5.bge[[i]]$StructureStr[[which.max(Freq_comb1.5.bge[[i]]$DAGscores)]],2,2),truemat))

compDAGs_comb2.bge <- lapply(1:nsim, 
                             function(i) compareDAGs(matrix(Freq_comb2.bge[[i]]$StructureStr[[which.max(Freq_comb2.bge[[i]]$DAGscores)]],2,2),truemat))

compDAGs_comb2.5.bge <- lapply(1:nsim, 
                             function(i) compareDAGs(matrix(Freq_comb2.5.bge[[i]]$StructureStr[[which.max(Freq_comb2.5.bge[[i]]$DAGscores)]],2,2),truemat))

compDAGs_comb3.bge <- lapply(1:nsim, 
                             function(i) compareDAGs(matrix(Freq_comb3.bge[[i]]$StructureStr[[which.max(Freq_comb3.bge[[i]]$DAGscores)]],2,2),truemat))

compDAGs_comb5.bge <- lapply(1:nsim, 
                             function(i) compareDAGs(matrix(Freq_comb5.bge[[i]]$StructureStr[[which.max(Freq_comb5.bge[[i]]$DAGscores)]],2,2),truemat))

compDAGs_comb10.bge <- lapply(1:nsim, 
                              function(i) compareDAGs(matrix(Freq_comb10.bge[[i]]$StructureStr[[which.max(Freq_comb10.bge[[i]]$DAGscores)]],2,2),truemat))

l <- length(compDAGs_comb0.05.bge)
# b = 0.05
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_comb0.05.bge[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Comb0.05_bge[[x]][[1]])))) # Avg no. of unique DAGs
# Number of correctly identified dags
length(which(sapply(1:nsim, function(i) identical(Comb0.05_bge[[i]][[2]],truemat)) == TRUE))

# b = 0.1
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_comb0.1.bge[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Comb0.1_bge[[x]][[1]])))) # Avg no. of unique DAGs
# Number of correctly identified dags
length(which(sapply(1:nsim, function(i) identical(Comb0.1_bge[[i]][[2]],truemat)) == TRUE))

# b = 0.5
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_comb0.5.bge[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Comb0.5_bge[[x]][[1]])))) # Avg no. of unique DAGs
# Number of correctly identified dags
length(which(sapply(1:nsim, function(i) identical(Comb0.5_bge[[i]][[2]],truemat)) == TRUE))

# b = 1
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_comb1.bge[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Comb1_bge[[x]][[1]])))) # Avg no. of unique DAGs
# Number of correctly identified dags
length(which(sapply(1:nsim, function(i) identical(Comb1_bge[[i]][[2]],truemat)) == TRUE))

# b=1.5
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_comb1.5.bge[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Comb1.5_bge[[x]][[1]])))) # Avg no. of unique DAGs
# Number of correctly identified dags
length(which(sapply(1:nsim, function(i) identical(Comb1.5_bge[[i]][[2]],truemat)) == TRUE))

# b=2
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_comb2.bge[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Comb2_bge[[x]][[1]])))) # Avg no. of unique DAGs
# Number of correctly identified dags
length(which(sapply(1:nsim, function(i) identical(Comb2_bge[[i]][[2]],truemat)) == TRUE))

# b=2.5
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_comb2.5.bge[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Comb2.5_bge[[x]][[1]])))) # Avg no. of unique DAGs
# Number of correctly identified dags
length(which(sapply(1:nsim, function(i) identical(Comb2.5_bge[[i]][[2]],truemat)) == TRUE))

# b=3
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_comb3.bge[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Comb3_bge[[x]][[1]])))) # Avg no. of unique DAGs
# Number of correctly identified dags
length(which(sapply(1:nsim, function(i) identical(Comb3_bge[[i]][[2]],truemat)) == TRUE))

# b=5
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_comb5.bge[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Comb5_bge[[x]][[1]])))) # Avg no. of unique DAGs
# Number of correctly identified dags
length(which(sapply(1:nsim, function(i) identical(Comb5_bge[[i]][[2]],truemat)) == TRUE))

# b=10
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_comb10.bge[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Comb10_bge[[x]][[1]])))) # Avg no. of unique DAGs
# Number of correctly identified dags
length(which(sapply(1:nsim, function(i) identical(Comb10_bge[[i]][[2]],truemat)) == TRUE))

# Discretised
CombDiscretised <- function(data){
  Discretised <- function(x) cut(x, breaks = c(min(x),(min(x)+max(x))/2,max(x)),
                                 labels = 0:1,
                                 include.lowest = TRUE)
  
  data[,1] <- Discretised((data %>% pull("A")))
  return(data)
} 

# List of datasets to be used for inference
# b=0.05
Comb0.05_dcrt_ds <- lapply(1:nsim, function(i) as.data.frame(CombDiscretised(Comb0.05_ds[[i]])))
Comb0.05_dcrt_ds <- lapply(1:nsim, function(x) apply(Comb0.05_dcrt_ds[[x]], 2, as.numeric))
# b=0.1
Comb0.1_dcrt_ds <- lapply(1:nsim, function(i) as.data.frame(CombDiscretised(Comb0.1_ds[[i]])))
Comb0.1_dcrt_ds <- lapply(1:nsim, function(x) apply(Comb0.1_dcrt_ds[[x]], 2, as.numeric))
# b=0.5
Comb0.5_dcrt_ds <- lapply(1:nsim, function(i) as.data.frame(CombDiscretised(Comb0.5_ds[[i]])))
Comb0.5_dcrt_ds <- lapply(1:nsim, function(x) apply(Comb0.5_dcrt_ds[[x]], 2, as.numeric))
# b=1
Comb1_dcrt_ds <- lapply(1:nsim, function(i) as.data.frame(CombDiscretised(Comb1_ds[[i]])))
Comb1_dcrt_ds <- lapply(1:nsim, function(x) apply(Comb1_dcrt_ds[[x]], 2, as.numeric))
# b=1.5
Comb1.5_dcrt_ds <- lapply(1:nsim, function(i) as.data.frame(CombDiscretised(Comb1.5_ds[[i]])))
Comb1.5_dcrt_ds <- lapply(1:nsim, function(x) apply(Comb1.5_dcrt_ds[[x]], 2, as.numeric))
# b=2
Comb2_dcrt_ds <- lapply(1:nsim, function(i) as.data.frame(CombDiscretised(Comb2_ds[[i]])))
Comb2_dcrt_ds <- lapply(1:nsim, function(x) apply(Comb2_dcrt_ds[[x]], 2, as.numeric))
# b=2.5
Comb2.5_dcrt_ds <- lapply(1:nsim, function(i) as.data.frame(CombDiscretised(Comb2.5_ds[[i]])))
Comb2.5_dcrt_ds <- lapply(1:nsim, function(x) apply(Comb2.5_dcrt_ds[[x]], 2, as.numeric))
# b=3
Comb3_dcrt_ds <- lapply(1:nsim, function(i) as.data.frame(CombDiscretised(Comb3_ds[[i]])))
Comb3_dcrt_ds <- lapply(1:nsim, function(x) apply(Comb3_dcrt_ds[[x]], 2, as.numeric))
# b=5
Comb5_dcrt_ds <- lapply(1:nsim, function(i) as.data.frame(CombDiscretised(Comb5_ds[[i]])))
Comb5_dcrt_ds <- lapply(1:nsim, function(x) apply(Comb5_dcrt_ds[[x]], 2, as.numeric))
# b=10
Comb10_dcrt_ds <- lapply(1:nsim, function(i) as.data.frame(CombDiscretised(Comb10_ds[[i]])))
Comb10_dcrt_ds <- lapply(1:nsim, function(x) apply(Comb10_dcrt_ds[[x]], 2, as.numeric))


# Partition MCMC
Comb0.05_bde <- lapply(1:nsim, function(i) Inference("bde", Comb0.05_dcrt_ds[[i]]))
Comb0.1_bde <- lapply(1:nsim, function(i) Inference("bde", Comb0.1_dcrt_ds[[i]]))
Comb0.5_bde <- lapply(1:nsim, function(i) Inference("bde", Comb0.5_dcrt_ds[[i]]))
Comb1_bde <- lapply(1:nsim, function(i) Inference("bde", Comb1_dcrt_ds[[i]]))
Comb1.5_bde <- lapply(1:nsim, function(i) Inference("bde", Comb1.5_dcrt_ds[[i]]))
Comb2_bde <- lapply(1:nsim, function(i) Inference("bde", Comb2_dcrt_ds[[i]]))
Comb2.5_bde <- lapply(1:nsim, function(i) Inference("bde", Comb2.5_dcrt_ds[[i]]))
Comb3_bde <- lapply(1:nsim, function(i) Inference("bde", Comb3_dcrt_ds[[i]]))
Comb5_bde <- lapply(1:nsim, function(i) Inference("bde", Comb5_dcrt_ds[[i]]))
Comb10_bde <- lapply(1:nsim, function(i) Inference("bde", Comb10_dcrt_ds[[i]]))


# Extracting frequency table for unique structures
Freq_comb0.01.bde <- lapply(1:nsim, function(i) Comb0.01_bde[[i]][[1]])
Freq_comb0.05.bde <- lapply(1:nsim, function(i) Comb0.05_bde[[i]][[1]])
Freq_comb0.1.bde <- lapply(1:nsim, function(i) Comb0.1_bde[[i]][[1]])
Freq_comb0.5.bde <- lapply(1:nsim, function(i) Comb0.5_bde[[i]][[1]])
Freq_comb1.bde <- lapply(1:nsim, function(i) Comb1_bde[[i]][[1]])
Freq_comb1.5.bde <- lapply(1:nsim, function(i) Comb1.5_bde[[i]][[1]])
Freq_comb2.bde <- lapply(1:nsim, function(i) Comb2_bde[[i]][[1]])
Freq_comb2.5.bde <- lapply(1:nsim, function(i) Comb2.5_bde[[i]][[1]])
Freq_comb3.bde <- lapply(1:nsim, function(i) Comb3_bde[[i]][[1]])
Freq_comb5.bde <- lapply(1:nsim, function(i) Comb5_bde[[i]][[1]])
Freq_comb10.bde <- lapply(1:nsim, function(i) Comb10_bde[[i]][[1]])

# Compare with true DAG
compDAGs_comb0.05.bde <- lapply(1:nsim, 
                         function(i) compareDAGs(matrix(Freq_comb0.05.bde[[i]]$StructureStr[[which.max(Freq_comb0.05.bde[[i]]$DAGscores)]],2,2),truemat))

compDAGs_comb0.1.bde <- lapply(1:nsim, 
                        function(i) compareDAGs(matrix(Freq_comb0.1.bde[[i]]$StructureStr[[which.max(Freq_comb0.1.bde[[i]]$DAGscores)]],2,2),truemat))

compDAGs_comb0.5.bde <- lapply(1:nsim, 
                        function(i) compareDAGs(matrix(Freq_comb0.5.bde[[i]]$StructureStr[[which.max(Freq_comb0.5.bde[[i]]$DAGscores)]],2,2),truemat))

compDAGs_comb1.bde <- lapply(1:nsim, 
                      function(i) compareDAGs(matrix(Freq_comb1.bde[[i]]$StructureStr[[which.max(Freq_comb1.bde[[i]]$DAGscores)]],2,2),truemat))

compDAGs_comb1.5.bde <- lapply(1:nsim, 
                        function(i) compareDAGs(matrix(Freq_comb1.5.bde[[i]]$StructureStr[[which.max(Freq_comb1.5.bde[[i]]$DAGscores)]],2,2),truemat))

compDAGs_comb2.bde <- lapply(1:nsim, 
                      function(i) compareDAGs(matrix(Freq_comb2.bde[[i]]$StructureStr[[which.max(Freq_comb2.bde[[i]]$DAGscores)]],2,2),truemat))

compDAGs_comb2.5.bde <- lapply(1:nsim, 
                        function(i) compareDAGs(matrix(Freq_comb2.5.bde[[i]]$StructureStr[[which.max(Freq_comb2.5.bde[[i]]$DAGscores)]],2,2),truemat))

compDAGs_comb3.bde <- lapply(1:nsim, 
                      function(i) compareDAGs(matrix(Freq_comb3.bde[[i]]$StructureStr[[which.max(Freq_comb3.bde[[i]]$DAGscores)]],2,2),truemat))

compDAGs_comb5.bde <- lapply(1:nsim, 
                             function(i) compareDAGs(matrix(Freq_comb5.bde[[i]]$StructureStr[[which.max(Freq_comb5.bde[[i]]$DAGscores)]],2,2),truemat))

compDAGs_comb10.bde <- lapply(1:nsim, 
                              function(i) compareDAGs(matrix(Freq_comb10.bde[[i]]$StructureStr[[which.max(Freq_comb10.bde[[i]]$DAGscores)]],2,2),truemat))

l <- length(compDAGs_comb2.bde)

# b=0.05
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_comb0.05.bde[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Comb0.05_bde[[x]][[1]])))) # Avg no. of unique DAGs
# Number of correctly identified dags
length(which(sapply(1:nsim, function(i) identical(Comb0.05_bde[[i]][[2]],truemat)) == TRUE))

# b=0.1
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_comb0.1.bde[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Comb0.1_bde[[x]][[1]])))) # Avg no. of unique DAGs
# Number of correctly identified dags
length(which(sapply(1:nsim, function(i) identical(Comb0.1_bde[[i]][[2]],truemat)) == TRUE))

# b=0.5
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_comb0.5.bde[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Comb0.5_bde[[x]][[1]])))) # Avg no. of unique DAGs
# Number of correctly identified dags
length(which(sapply(1:nsim, function(i) identical(Comb0.5_bde[[i]][[2]],truemat)) == TRUE))

# b=1
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_comb1.bde[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Comb1_bde[[x]][[1]])))) # Avg no. of unique DAGs
# Number of correctly identified dags
length(which(sapply(1:nsim, function(i) identical(Comb1_bde[[i]][[2]],truemat)) == TRUE))

# b=1.5
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_comb1.5.bde[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Comb1.5_bde[[x]][[1]])))) # Avg no. of unique DAGs
# Number of correctly identified dags
length(which(sapply(1:nsim, function(i) identical(Comb1.5_bde[[i]][[2]],truemat)) == TRUE))

# b=2
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_comb2.bde[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Comb2_bde[[x]][[1]])))) # Avg no. of unique DAGs
# Number of correctly identified dags
length(which(sapply(1:nsim, function(i) identical(Comb2_bde[[i]][[2]],truemat)) == TRUE))

# b=2.5
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_comb2.5.bde[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Comb2.5_bde[[x]][[1]])))) # Avg no. of unique DAGs
# Number of correctly identified dags
length(which(sapply(1:nsim, function(i) identical(Comb2.5_bde[[i]][[2]],truemat)) == TRUE))

# b=3
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_comb3.bde[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Comb3_bde[[x]][[1]])))) # Avg no. of unique DAGs
# Number of correctly identified dags
length(which(sapply(1:nsim, function(i) identical(Comb3_bde[[i]][[2]],truemat)) == TRUE))

# b=5
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_comb5.bde[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Comb5_bde[[x]][[1]])))) # Avg no. of unique DAGs
# Number of correctly identified dags
length(which(sapply(1:nsim, function(i) identical(Comb5_bde[[i]][[2]],truemat)) == TRUE))


# b=10
sapply(1:8, function (x) mean(as.numeric(sapply(1:l, function(i) compDAGs_comb10.bde[[i]][x])))) # Compare to True DAGs
mean(as.numeric(lapply(1:l, function (x) nrow(Comb10_bde[[x]][[1]])))) # Avg no. of unique DAGs
# Number of correctly identified dags
length(which(sapply(1:nsim, function(i) identical(Comb10_bde[[i]][[2]],truemat)) == TRUE))



# Freq ratio

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


F_ratio(Freq_comb0.1.bge)
F_ratio(Freq_comb0.1.bde)