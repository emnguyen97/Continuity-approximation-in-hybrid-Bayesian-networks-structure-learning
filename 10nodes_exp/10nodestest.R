# Libraries
# Packages such as 'graph', 'Rgraphviz' and 'RBGL' that are required for BiDAG and some others are no longer available from CRAN. 
# You must use BiocManager to install them. This needs to be done manually. eg: BiocManager::install("Rgraphviz")
# Libraries

library(BiDAG)
library(bnlearn)
library(tibble)
library(ggdag)
library(dagitty)
library(bnlearn)
library(reshape2)
library(igraph)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork) # For arranging multiple plots with equal panel sizes
library(RColorBrewer)

# Load the necessary functions
source('10nodes/MethodsUtil.R')
source('10nodes/DataGenUtil.R')
source('structureMCMC/combinations.R')
source('structureMCMC/scoretables.R')
source('structureMCMC/structurefns.R')
source('structureMCMC/samplefns.R')
source('structureMCMC/param_utils.R')
source('structureMCMC/structureMCMC.R')


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

numCores <- detectCores()

nrep <- 200 #No of observations
nvar <- 10 #No of variables
nsim <- 100 #No of data replicates

##### S_cc #####################################################################
set.seed(1097)
# True adj matrix
datmat <- rDAG(nvar,0.2)
#datmat[10,2] <- datmat[10,3] <- datmat[10,8] <- 0
datmat[10,3] <- 0
datmat[5,2] <- datmat[5,6] <- datmat[2,7] <- 1
plotadj(datmat)

# List of datasets to be used for inference
W_mat <- datmat*matrix(runif(nvar*nvar, 1, 5), nvar, nvar)*sample(c(-1, 1), size = nvar*nvar, replace = TRUE)
D <- diag(runif(nvar, 0, 1))
colnames(D) <- rownames(D) <- colnames(datmat)
set.seed(2023)
datasets <- lapply(1:nsim, function(i) Scc_Data(datmat, W_mat, D, nrep, nvar))

### RAG - bge
#1. TABU
tabur <- lapply(1:nsim, function(x) runtabu(as.data.frame(datasets[[x]]), datmat, "bge"))
get_result(tabur,datmat)

#2. HC
hcr <- lapply(1:nsim, function(x) runhc(as.data.frame(datasets[[x]]), datmat, "bge"))
get_result(hcr,datmat)

#3. Partition MCMC

partMCMC <- list()
cl <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl) 

partMCMC <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(datasets[[i]])
  result$result <- runpartMCMC(data, 'bge')
  return(result)
}
stopCluster(cl)

options(digits = 4)
get_result(partMCMC,datmat)

### DISC
contnodes_Scc <- colnames(datasets[[1]])

runDiscretise <- function(datasets, contnodes, runtabu, runhc, get_result, runpartMCMC, multiResultClass){
  
  # Discretise datasets
  datasets_interval_2b <- lapply(1:nsim, function(x) Discretised(datasets[[x]], contnodes, method="interval", ncategory=2))
  datasets_interval_4b <- lapply(1:nsim, function(x) Discretised(datasets[[x]], contnodes, method="interval", ncategory=4))
  datasets_clust_2b <- lapply(1:nsim, function(x) Discretised(datasets[[x]], contnodes, method="cluster", ncategory=2))
  datasets_clust_4b <- lapply(1:nsim, function(x) Discretised(datasets[[x]], contnodes, method="cluster", ncategory=4))
  
  ## Strategy: equal interval width
  # TABU
  tabur_interval_2b <- lapply(1:nsim, function(x) runtabu(as.data.frame(datasets_interval_2b[[x]]), datmat, "bde"))
  print(get_result(tabur_interval_2b,datmat))
  
  tabur_interval_4b <- lapply(1:nsim, function(x) runtabu(as.data.frame(datasets_interval_4b[[x]]), datmat, "bde"))
  print(get_result(tabur_interval_4b,datmat))
  
  # HC
  hcr_interval_2b <- lapply(1:nsim, function(x) runhc(as.data.frame(datasets_interval_2b[[x]]), datmat, "bde"))
  print(get_result(hcr_interval_2b,datmat))
  
  hcr_interval_4b <- lapply(1:nsim, function(x) runhc(as.data.frame(datasets_interval_4b[[x]]), datmat, "bde"))
  print(get_result(hcr_interval_4b,datmat))
  
  # PMCMC
  partMCMC_interval_2b <- list()
  partMCMC_interval_4b <- list()
  cl <- parallel::makeCluster(numCores)
  doParallel::registerDoParallel(cl) 
  
  partMCMC_interval_2b <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
    result <- multiResultClass()
    data <- as.data.frame(datasets_interval_2b[[i]])
    result$result <- runpartMCMC(data, 'bdecat')
    return(result)
  }
  
  partMCMC_interval_4b <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
    result <- multiResultClass()
    data <- as.data.frame(datasets_interval_4b[[i]])
    result$result <- runpartMCMC(data, 'bdecat')
    return(result)
  }
  stopCluster(cl)
  
  print(get_result(partMCMC_interval_2b,datmat))
  print(get_result(partMCMC_interval_4b,datmat))
  
  ## Strategy: k-means clustering discretisation
  # TABU
  tabur_clust_2b <- lapply(1:nsim, function(x) runtabu(as.data.frame(datasets_clust_2b[[x]]), datmat, "bde"))
  print(get_result(tabur_clust_2b,datmat))
  
  tabur_clust_4b <- lapply(1:nsim, function(x) runtabu(as.data.frame(datasets_clust_4b[[x]]), datmat, "bde"))
  print(get_result(tabur_clust_4b,datmat))
  
  # HC
  hcr_clust_2b <- lapply(1:nsim, function(x) runhc(as.data.frame(datasets_clust_2b[[x]]), datmat, "bde"))
  print(get_result(hcr_clust_2b,datmat))
  
  hcr_clust_4b <- lapply(1:nsim, function(x) runhc(as.data.frame(datasets_clust_4b[[x]]), datmat, "bde"))
  print(get_result(hcr_clust_4b,datmat))
  
  # PMCMC
  partMCMC_clust_2b <- list()
  partMCMC_clust_4b <- list()
  cl <- parallel::makeCluster(numCores)
  doParallel::registerDoParallel(cl) 
  
  partMCMC_clust_2b <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
    result <- multiResultClass()
    data <- as.data.frame(datasets_clust_2b[[i]])
    result$result <- runpartMCMC(data, 'bdecat')
    return(result)
  }
  
  partMCMC_clust_4b <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
    result <- multiResultClass()
    data <- as.data.frame(datasets_clust_4b[[i]])
    result$result <- runpartMCMC(data, 'bdecat')
    return(result)
  }
  stopCluster(cl)
  
  print(get_result(partMCMC_clust_2b,datmat))
  print(get_result(partMCMC_clust_4b,datmat))
  
  return(list(tabur_interval_2b=tabur_interval_2b, tabur_interval_4b=tabur_interval_4b, hcr_interval_2b=hcr_interval_2b, hcr_interval_4b=hcr_interval_4b, partMCMC_interval_2b=partMCMC_interval_2b, partMCMC_interval_4b=partMCMC_interval_4b, tabur_clust_2b=tabur_clust_2b, tabur_clust_4b=tabur_clust_4b, hcr_clust_2b=hcr_clust_2b, hcr_clust_4b= hcr_clust_4b, partMCMC_clust_2b=partMCMC_clust_2b, partMCMC_clust_4b=partMCMC_clust_4b))
}

DCRT_Scc <- runDiscretise(datasets, contnodes_Scc, runtabu, runhc, get_result, runpartMCMC, multiResultClass)

tabur_Scc_interval_2b <- DCRT_Scc$tabur_interval_2b
tabur_Scc_interval_4b <- DCRT_Scc$tabur_interval_4b
hcr_Scc_interval_2b <- DCRT_Scc$hcr_interval_2b
hcr_Scc_interval_4b <- DCRT_Scc$hcr_interval_4b
partMCMC_Scc_interval_2b <- DCRT_Scc$partMCMC_interval_2b
partMCMC_Scc_interval_4b <- DCRT_Scc$partMCMC_interval_4b

tabur_Scc_clust_2b <- DCRT_Scc$tabur_clust_2b
tabur_Scc_clust_4b <- DCRT_Scc$tabur_clust_4b
hcr_Scc_clust_2b <- DCRT_Scc$hcr_clust_2b
hcr_Scc_clust_4b <- DCRT_Scc$hcr_clust_4b
partMCMC_Scc_clust_2b <- DCRT_Scc$partMCMC_clust_2b
partMCMC_Scc_clust_4b <- DCRT_Scc$partMCMC_clust_4b


### RAG - Likelihood
# TABU
tabur_like <- lapply(1:nsim, function(x) runtabu(as.data.frame(datasets[[x]]), datmat, "bic-g"))
get_result(tabur_like,datmat)

# HC
hcr_like <- lapply(1:nsim, function(x) runhc(as.data.frame(datasets[[x]]), datmat, "bic-g"))
get_result(hcr_like,datmat)

# Structure MCMC
structMCMC <- list()
cl <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl) 

structMCMC <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  tryCatch({
    result <- multiResultClass()
    data <- as.data.frame(datasets[[i]])
    result$result <- runStructMCMC(data,iterations = 100000,blklist=NULL,scoretype="bic-g",sample_parameters=FALSE)
    return(result)
  }, error = function(e) {
    message(sprintf("Error in iteration %d: %s", i, e$message))
    return(NULL)
  })
}

stopCluster(cl)

options(digits = 4)
get_result(structMCMC,datmat)

### Tabulate results

S_cc_SHD <- metric_Scc("SHD", datmat, tabur, hcr, partMCMC, tabur_like, hcr_like, structMCMC, tabur_Scc_interval_2b, hcr_Scc_interval_2b, partMCMC_Scc_interval_2b, tabur_Scc_interval_4b, hcr_Scc_interval_4b, partMCMC_Scc_interval_4b, tabur_Scc_clust_2b, hcr_Scc_clust_2b, partMCMC_Scc_clust_2b, tabur_Scc_clust_4b, hcr_Scc_clust_4b, partMCMC_Scc_clust_4b)

S_cc_TPR <- metric_Scc("TPR", datmat, tabur, hcr, partMCMC, tabur_like, hcr_like, structMCMC, tabur_Scc_interval_2b, hcr_Scc_interval_2b, partMCMC_Scc_interval_2b, tabur_Scc_interval_4b, hcr_Scc_interval_4b, partMCMC_Scc_interval_4b, tabur_Scc_clust_2b, hcr_Scc_clust_2b, partMCMC_Scc_clust_2b, tabur_Scc_clust_4b, hcr_Scc_clust_4b, partMCMC_Scc_clust_4b)

S_cc_TP <- metric_Scc("TP", datmat, tabur, hcr, partMCMC, tabur_like, hcr_like, structMCMC, tabur_Scc_interval_2b, hcr_Scc_interval_2b, partMCMC_Scc_interval_2b, tabur_Scc_interval_4b, hcr_Scc_interval_4b, partMCMC_Scc_interval_4b, tabur_Scc_clust_2b, hcr_Scc_clust_2b, partMCMC_Scc_clust_2b, tabur_Scc_clust_4b, hcr_Scc_clust_4b, partMCMC_Scc_clust_4b)

S_cc_FP <- metric_Scc("FP", datmat, tabur, hcr, partMCMC, tabur_like, hcr_like, structMCMC, tabur_Scc_interval_2b, hcr_Scc_interval_2b, partMCMC_Scc_interval_2b, tabur_Scc_interval_4b, hcr_Scc_interval_4b, partMCMC_Scc_interval_4b, tabur_Scc_clust_2b, hcr_Scc_clust_2b, partMCMC_Scc_clust_2b, tabur_Scc_clust_4b, hcr_Scc_clust_4b, partMCMC_Scc_clust_4b)

S_cc_FN <- metric_Scc("FN", datmat, tabur, hcr, partMCMC, tabur_like, hcr_like, structMCMC, tabur_Scc_interval_2b, hcr_Scc_interval_2b, partMCMC_Scc_interval_2b, tabur_Scc_interval_4b, hcr_Scc_interval_4b, partMCMC_Scc_interval_4b, tabur_Scc_clust_2b, hcr_Scc_clust_2b, partMCMC_Scc_clust_2b, tabur_Scc_clust_4b, hcr_Scc_clust_4b, partMCMC_Scc_clust_4b)

S_cc_SHD$F1 <- F1score(S_cc_TP$metric, S_cc_FP$metric, S_cc_FN$metric)

plot_metric(S_cc_SHD, "SHD") 
plot_metric(S_cc_TPR, "TPR") 

##### S_cd #####################################################################
# Only children nodes are discrete

# List of datasets to be used for inference
set.seed(23)
# Discrete nodes are stored as factors
datasets_Scd <- lapply(1:nsim, function(i) Scd_Data(datmat, W_mat, D, nrep, nvar))
# All nodes are stored as numeric
datasets_Scd_num <- lapply(1:nsim, function(x) datasets_Scd[[x]] %>% mutate_if(sapply(datasets_Scd[[x]], is.factor), as.numeric))

### RAG
# TABU
tabur_Scd <- lapply(1:nsim, function(x) runtabu(as.data.frame(datasets_Scd_num[[x]]), datmat, "bge"))
get_result(tabur_Scd,datmat)

# HC
hcr_Scd <- lapply(1:nsim, function(x) runhc(as.data.frame(datasets_Scd_num[[x]]), datmat, "bge"))
get_result(hcr_Scd,datmat)

# PMCMC
partMCMC_Scd <- list()
cl <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl) 

partMCMC_Scd <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(datasets_Scd_num[[i]])
  result$result <- runpartMCMC(data, 'bge')
  return(result)
}
stopCluster(cl)

options(digits = 4)
get_result(partMCMC_Scd,datmat) 

### CLG
tabur_Scd_cg <- lapply(1:nsim, function(x) runtabu(as.data.frame(datasets_Scd[[x]]), datmat, "bic-cg"))
get_result(tabur_Scd_cg,datmat)

hcr_Scd_cg <- lapply(1:nsim, function(x) runhc(as.data.frame(datasets_Scd[[x]]), datmat, "bic-cg"))
get_result(hcr_Scd_cg,datmat)

dcrtnodes_Scd <- c("A","C","G","H")
contnodes_Scd <- setdiff(colnames(datmat),dcrtnodes_Scd)
Scd_blklist <- matrix(0, nvar, nvar)
colnames(Scd_blklist) <- rownames(Scd_blklist) <- colnames(datmat)
diag(Scd_blklist) <- 1
Scd_blklist[contnodes_Scd,dcrtnodes_Scd] <- 1

structMCMC_Scd_clg <- list()
cl <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl) 

structMCMC_Scd_clg <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  tryCatch({
    result <- multiResultClass()
    data <- as.data.frame(datasets_Scd[[i]])
    result$result <- runStructMCMC(data, iterations = 100000, blklist = Scd_blklist, scoretype = "bic-cg", sample_parameters = FALSE)
    return(result)
  }, error = function(e) {
    message(sprintf("Error in iteration %d: %s", i, e$message))
    return(NULL)
  })
}

stopCluster(cl)

options(digits = 4)
get_result(structMCMC_Scd_clg, datmat)

### Discretise

DCRT_Scd <- runDiscretise(datasets_Scd, contnodes_Scd, runtabu, runhc, get_result, runpartMCMC, multiResultClass)

tabur_Scd_interval_2b <- DCRT_Scd$tabur_interval_2b
tabur_Scd_interval_4b <- DCRT_Scd$tabur_interval_4b
hcr_Scd_interval_2b <- DCRT_Scd$hcr_interval_2b
hcr_Scd_interval_4b <- DCRT_Scd$hcr_interval_4b
partMCMC_Scd_interval_2b <- DCRT_Scd$partMCMC_interval_2b
partMCMC_Scd_interval_4b <- DCRT_Scd$partMCMC_interval_4b

tabur_Scd_clust_2b <- DCRT_Scd$tabur_clust_2b
tabur_Scd_clust_4b <- DCRT_Scd$tabur_clust_4b
hcr_Scd_clust_2b <- DCRT_Scd$hcr_clust_2b
hcr_Scd_clust_4b <- DCRT_Scd$hcr_clust_4b
partMCMC_Scd_clust_2b <- DCRT_Scd$partMCMC_clust_2b
partMCMC_Scd_clust_4b <- DCRT_Scd$partMCMC_clust_4b

### RAG - Likelihood
# TABU
tabur_Scd_like <- lapply(1:nsim, function(x) runtabu(as.data.frame(datasets_Scd_num[[x]]), datmat, "bic-g"))
get_result(tabur_Scd_like,datmat)

# HC
hcr_Scd_like <- lapply(1:nsim, function(x) runhc(as.data.frame(datasets_Scd_num[[x]]), datmat, "bic-g"))
get_result(hcr_Scd_like,datmat)

# Structure MCMC 
structMCMC_Scd <- list()
cl <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl) 

structMCMC_Scd <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(datasets_Scd_num[[i]])
  result$result <- runStructMCMC(data,iterations = 100000,blklist=NULL,scoretype="bic-g",sample_parameters=FALSE)
  return(result)
}

stopCluster(cl)

options(digits = 4)
get_result(structMCMC_Scd,datmat)


### Tabulate results

S_cd_SHD <- metric_comb("SHD", datmat, tabur_Scd, hcr_Scd, partMCMC_Scd, tabur_Scd_like, hcr_Scd_like, structMCMC_Scd, tabur_Scd_cg, hcr_Scd_cg, structMCMC_Scd_clg, tabur_Scd_interval_2b, hcr_Scd_interval_2b, partMCMC_Scd_interval_2b, tabur_Scd_interval_4b, hcr_Scd_interval_4b, partMCMC_Scd_interval_4b, tabur_Scd_clust_2b, hcr_Scd_clust_2b, partMCMC_Scd_clust_2b, tabur_Scd_clust_4b, hcr_Scd_clust_4b, partMCMC_Scd_clust_4b)

S_cd_TP <- metric_comb("TP", datmat, tabur_Scd, hcr_Scd, partMCMC_Scd, tabur_Scd_like, hcr_Scd_like, structMCMC_Scd, tabur_Scd_cg, hcr_Scd_cg, structMCMC_Scd_clg, tabur_Scd_interval_2b, hcr_Scd_interval_2b, partMCMC_Scd_interval_2b, tabur_Scd_interval_4b, hcr_Scd_interval_4b, partMCMC_Scd_interval_4b, tabur_Scd_clust_2b, hcr_Scd_clust_2b, partMCMC_Scd_clust_2b, tabur_Scd_clust_4b, hcr_Scd_clust_4b, partMCMC_Scd_clust_4b)

S_cd_FP <- metric_comb("FP", datmat, tabur_Scd, hcr_Scd, partMCMC_Scd, tabur_Scd_like, hcr_Scd_like, structMCMC_Scd, tabur_Scd_cg, hcr_Scd_cg, structMCMC_Scd_clg, tabur_Scd_interval_2b, hcr_Scd_interval_2b, partMCMC_Scd_interval_2b, tabur_Scd_interval_4b, hcr_Scd_interval_4b, partMCMC_Scd_interval_4b, tabur_Scd_clust_2b, hcr_Scd_clust_2b, partMCMC_Scd_clust_2b, tabur_Scd_clust_4b, hcr_Scd_clust_4b, partMCMC_Scd_clust_4b)

S_cd_FN <- metric_comb("FN", datmat, tabur_Scd, hcr_Scd, partMCMC_Scd, tabur_Scd_like, hcr_Scd_like, structMCMC_Scd, tabur_Scd_cg, hcr_Scd_cg, structMCMC_Scd_clg, tabur_Scd_interval_2b, hcr_Scd_interval_2b, partMCMC_Scd_interval_2b, tabur_Scd_interval_4b, hcr_Scd_interval_4b, partMCMC_Scd_interval_4b, tabur_Scd_clust_2b, hcr_Scd_clust_2b, partMCMC_Scd_clust_2b, tabur_Scd_clust_4b, hcr_Scd_clust_4b, partMCMC_Scd_clust_4b)

S_cd_TPR <- metric_comb("TPR", datmat, tabur_Scd, hcr_Scd, partMCMC_Scd, tabur_Scd_like, hcr_Scd_like, structMCMC_Scd, tabur_Scd_cg, hcr_Scd_cg, structMCMC_Scd_clg, tabur_Scd_interval_2b, hcr_Scd_interval_2b, partMCMC_Scd_interval_2b, tabur_Scd_interval_4b, hcr_Scd_interval_4b, partMCMC_Scd_interval_4b, tabur_Scd_clust_2b, hcr_Scd_clust_2b, partMCMC_Scd_clust_2b, tabur_Scd_clust_4b, hcr_Scd_clust_4b, partMCMC_Scd_clust_4b)

S_cd_SHD$F1 <- F1score(S_cd_TP$metric, S_cd_FP$metric, S_cd_FN$metric)

plot_metric(S_cd_SHD, "SHD") 
plot_metric(S_cd_TPR, "TPR") 
# plot_metric(S_cd_F1, "F1") 

##### S_dc #####################################################################
# Only children nodes are continuous

# List of datasets to be used for inference
set.seed(11)
# All nodes are stored as numeric
datasets_Sdc_num <- lapply(1:nsim, function(i) Sdc_Data(datmat, W_mat, D, nrep, nvar))
datasets_Sdc_num <- lapply(1:nsim, function(x) datasets_Sdc_num[[x]] %>% mutate_if(sapply(datasets_Sdc_num[[x]], is.integer), as.numeric))
# Discrete nodes are stored as factors
contnodes_dc <- c("A","C","D","G","H")
dcrtnodes_dc <- setdiff(colnames(datmat),contnodes_dc)
datasets_Sdc <- lapply(1:nsim, function(x) datasets_Sdc_num[[x]] %>% mutate_at(dcrtnodes_dc, as.factor))


### RAG - bge
# TABU
tabur_Sdc <- lapply(1:nsim, function(x) runtabu(as.data.frame(datasets_Sdc_num[[x]]), datmat, "bge"))
get_result(tabur_Sdc,datmat)

# HC
hcr_Sdc <- lapply(1:nsim, function(x) runhc(as.data.frame(datasets_Sdc_num[[x]]), datmat, "bge"))
get_result(hcr_Sdc,datmat)


# PMCMC
partMCMC_Sdc <- list()
cl <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl) 

partMCMC_Sdc <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(datasets_Sdc_num[[i]])
  result$result <- runpartMCMC(data, 'bge')
  return(result)
}
stopCluster(cl)

options(digits = 4)
get_result(partMCMC_Sdc,datmat)

### CLG
tabur_Sdc_cg <- lapply(1:nsim, function(x) runtabu(as.data.frame(datasets_Sdc[[x]]), datmat, "bic-cg"))
get_result(tabur_Sdc_cg,datmat)

hcr_Sdc_cg <- lapply(1:nsim, function(x) runhc(as.data.frame(datasets_Sdc[[x]]), datmat, "bic-cg"))
get_result(hcr_Sdc_cg,datmat)

Sdc_blklist <- matrix(0, nvar, nvar)
colnames(Sdc_blklist) <- rownames(Sdc_blklist) <- colnames(datmat)
diag(Sdc_blklist) <- 1
Sdc_blklist[contnodes_dc,dcrtnodes_dc] <- 1

structMCMC_Sdc_clg <- list()
cl <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl) 

structMCMC_Sdc_clg <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  tryCatch({
    result <- multiResultClass()
    data <- as.data.frame(datasets_Sdc[[i]])
    result$result <- runStructMCMC(data,iterations = 100000,blklist=Sdc_blklist,scoretype="bic-cg",sample_parameters=FALSE)
    return(result)
  }, error = function(e) {
    message(sprintf("Error in iteration %d: %s", i, e$message))
    return(NULL)
  })
}

stopCluster(cl)

options(digits = 4)
get_result(structMCMC_Sdc_clg,datmat)


### Discretise
DCRT_Sdc <- runDiscretise(datasets_Sdc, contnodes_dc, runtabu, runhc, get_result, runpartMCMC, multiResultClass)

tabur_Sdc_interval_2b <- DCRT_Sdc$tabur_interval_2b
tabur_Sdc_interval_4b <- DCRT_Sdc$tabur_interval_4b
hcr_Sdc_interval_2b <- DCRT_Sdc$hcr_interval_2b
hcr_Sdc_interval_4b <- DCRT_Sdc$hcr_interval_4b
partMCMC_Sdc_interval_2b <- DCRT_Sdc$partMCMC_interval_2b
partMCMC_Sdc_interval_4b <- DCRT_Sdc$partMCMC_interval_4b

tabur_Sdc_clust_2b <- DCRT_Sdc$tabur_clust_2b
tabur_Sdc_clust_4b <- DCRT_Sdc$tabur_clust_4b
hcr_Sdc_clust_2b <- DCRT_Sdc$hcr_clust_2b
hcr_Sdc_clust_4b <- DCRT_Sdc$hcr_clust_4b
partMCMC_Sdc_clust_2b <- DCRT_Sdc$partMCMC_clust_2b
partMCMC_Sdc_clust_4b <- DCRT_Sdc$partMCMC_clust_4b

## RAG - Likelihood
# TABU
tabur_Sdc_like <- lapply(1:nsim, function(x) runtabu(as.data.frame(datasets_Sdc_num[[x]]), datmat, "bic-g"))
get_result(tabur_Sdc_like,datmat)

# HC
hcr_Sdc_like <- lapply(1:nsim, function(x) runhc(as.data.frame(datasets_Sdc_num[[x]]), datmat, "bic-g"))
get_result(hcr_Sdc_like,datmat)

# Structure MCMC 
structMCMC_Sdc <- list()
cl <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl) 

structMCMC_Sdc <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(datasets_Sdc_num[[i]])
  result$result <- runStructMCMC(data,iterations = 100000,blklist=NULL,scoretype="bic-g",sample_parameters=FALSE)
  return(result)
}

stopCluster(cl)

options(digits = 4)
get_result(structMCMC_Sdc,datmat)


### Tabulate results
S_dc_SHD <- metric_comb("SHD", datmat, tabur_Sdc, hcr_Sdc, partMCMC_Sdc, tabur_Sdc_like, hcr_Sdc_like, structMCMC_Sdc, tabur_Sdc_cg, hcr_Sdc_cg, structMCMC_Sdc_clg, tabur_Sdc_interval_2b, hcr_Sdc_interval_2b, partMCMC_Sdc_interval_2b, tabur_Sdc_interval_4b, hcr_Sdc_interval_4b, partMCMC_Sdc_interval_4b, tabur_Sdc_clust_2b, hcr_Sdc_clust_2b, partMCMC_Sdc_clust_2b, tabur_Sdc_clust_4b, hcr_Sdc_clust_4b, partMCMC_Sdc_clust_4b)

S_dc_TPR <- metric_comb("TPR", datmat, tabur_Sdc, hcr_Sdc, partMCMC_Sdc, tabur_Sdc_like, hcr_Sdc_like, structMCMC_Sdc, tabur_Sdc_cg, hcr_Sdc_cg, structMCMC_Sdc_clg, tabur_Sdc_interval_2b, hcr_Sdc_interval_2b, partMCMC_Sdc_interval_2b, tabur_Sdc_interval_4b, hcr_Sdc_interval_4b, partMCMC_Sdc_interval_4b, tabur_Sdc_clust_2b, hcr_Sdc_clust_2b, partMCMC_Sdc_clust_2b, tabur_Sdc_clust_4b, hcr_Sdc_clust_4b, partMCMC_Sdc_clust_4b)

S_dc_TP <- metric_comb("TP", datmat, tabur_Sdc, hcr_Sdc, partMCMC_Sdc, tabur_Sdc_like, hcr_Sdc_like, structMCMC_Sdc, tabur_Sdc_cg, hcr_Sdc_cg, structMCMC_Sdc_clg, tabur_Sdc_interval_2b, hcr_Sdc_interval_2b, partMCMC_Sdc_interval_2b, tabur_Sdc_interval_4b, hcr_Sdc_interval_4b, partMCMC_Sdc_interval_4b, tabur_Sdc_clust_2b, hcr_Sdc_clust_2b, partMCMC_Sdc_clust_2b, tabur_Sdc_clust_4b, hcr_Sdc_clust_4b, partMCMC_Sdc_clust_4b)

S_dc_FP <- metric_comb("FP", datmat, tabur_Sdc, hcr_Sdc, partMCMC_Sdc, tabur_Sdc_like, hcr_Sdc_like, structMCMC_Sdc, tabur_Sdc_cg, hcr_Sdc_cg, structMCMC_Sdc_clg, tabur_Sdc_interval_2b, hcr_Sdc_interval_2b, partMCMC_Sdc_interval_2b, tabur_Sdc_interval_4b, hcr_Sdc_interval_4b, partMCMC_Sdc_interval_4b, tabur_Sdc_clust_2b, hcr_Sdc_clust_2b, partMCMC_Sdc_clust_2b, tabur_Sdc_clust_4b, hcr_Sdc_clust_4b, partMCMC_Sdc_clust_4b)

S_dc_FN <- metric_comb("FN", datmat, tabur_Sdc, hcr_Sdc, partMCMC_Sdc, tabur_Sdc_like, hcr_Sdc_like, structMCMC_Sdc, tabur_Sdc_cg, hcr_Sdc_cg, structMCMC_Sdc_clg, tabur_Sdc_interval_2b, hcr_Sdc_interval_2b, partMCMC_Sdc_interval_2b, tabur_Sdc_interval_4b, hcr_Sdc_interval_4b, partMCMC_Sdc_interval_4b, tabur_Sdc_clust_2b, hcr_Sdc_clust_2b, partMCMC_Sdc_clust_2b, tabur_Sdc_clust_4b, hcr_Sdc_clust_4b, partMCMC_Sdc_clust_4b)


S_dc_SHD$F1 <- F1score(S_dc_TP$metric, S_dc_FP$metric, S_dc_FN$metric)


plot_metric(S_dc_SHD, "SHD")
plot_metric(S_dc_TPR, "TPR")

##### S_dd #####################################################################  

# List of datasets to be used for inference
set.seed(2024)
# All nodes are stored as numeric
datasets_Sdd_num <- lapply(1:nsim, function(i) Sdd_Data(datmat, W_mat, nrep, nvar))
datasets_Sdd_num <- lapply(1:nsim, function(x) datasets_Sdd_num[[x]] %>% mutate_if(sapply(datasets_Sdd_num[[x]], is.integer), as.numeric))
# All nodes are stored as factors
datasets_Sdd <- lapply(1:nsim, function(x) datasets_Sdd_num[[x]] %>% mutate_if(sapply(datasets_Sdd_num[[x]], is.numeric), as.factor))

## RAG

# TABU
tabur_Sdd <- lapply(1:nsim, function(x) runtabu(as.data.frame(datasets_Sdd_num[[x]]), datmat, "bge"))
get_result(tabur_Sdd,datmat)

# HC
hcr_Sdd <- lapply(1:nsim, function(x) runhc(as.data.frame(datasets_Sdd_num[[x]]), datmat, "bge"))
get_result(hcr_Sdd,datmat)

# PMCMC
partMCMC_Sdd <- list()
cl <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl) 

partMCMC_Sdd <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(datasets_Sdd_num[[i]])
  result$result <- runpartMCMC(data, 'bge')
  return(result)
}
stopCluster(cl)

options(digits = 4)
get_result(partMCMC_Sdd,datmat)

## Discretised
# TABU
tabur_Sdd_dcrt <- lapply(1:nsim, function(x) runtabu(as.data.frame(datasets_Sdd[[x]]), datmat, "bde"))
get_result(tabur_Sdd_dcrt,datmat)

# HC
hcr_Sdd_dcrt <- lapply(1:nsim, function(x) runhc(as.data.frame(datasets_Sdd[[x]]), datmat, "bde"))
get_result(hcr_Sdd_dcrt,datmat)

# PMCMC
partMCMC_Sdd_dcrt <- list()
cl <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl) 

partMCMC_Sdd_dcrt <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(datasets_Sdd[[i]])
  result$result <- runpartMCMC(data, 'bdecat')
  return(result)
}
stopCluster(cl)

options(digits = 4)
get_result(partMCMC_Sdd_dcrt,datmat)

# RAG - Likelihood
# TABU
tabur_Sdd_like <- lapply(1:nsim, function(x) runtabu(as.data.frame(datasets_Sdd_num[[x]]), datmat, "bic-g"))
get_result(tabur_Sdd_like,datmat)

# HC
hcr_Sdd_like <- lapply(1:nsim, function(x) runhc(as.data.frame(datasets_Sdd_num[[x]]), datmat, "bic-g"))
get_result(hcr_Sdd_like,datmat)

# Structure MCMC 
structMCMC_Sdd <- list()
cl <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl) 

structMCMC_Sdd <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  tryCatch({
    result <- multiResultClass()
    data <- as.data.frame(datasets_Sdd_num[[i]])
    result$result <- runStructMCMC(data,iterations = 100000,blklist=NULL,scoretype="bic-g",sample_parameters=FALSE)
    return(result)
  }, error = function(e) {
    message(sprintf("Error in iteration %d: %s", i, e$message))
    return(NULL)
  })
}

stopCluster(cl)

options(digits = 4)
get_result(structMCMC_Sdd,datmat)


### Tabulate results

S_dd_SHD <- metric_Sdd("SHD", datmat, tabur_Sdd, hcr_Sdd, partMCMC_Sdd, tabur_Sdd_like, hcr_Sdd_like, structMCMC_Sdd, tabur_Sdd_dcrt, hcr_Sdd_dcrt, partMCMC_Sdd_dcrt)
S_dd_TPR <- metric_Sdd("TPR", datmat, tabur_Sdd, hcr_Sdd, partMCMC_Sdd, tabur_Sdd_like, hcr_Sdd_like, structMCMC_Sdd, tabur_Sdd_dcrt, hcr_Sdd_dcrt, partMCMC_Sdd_dcrt)

S_dd_TP <- metric_Sdd("TP", datmat, tabur_Sdd, hcr_Sdd, partMCMC_Sdd, tabur_Sdd_like, hcr_Sdd_like, structMCMC_Sdd, tabur_Sdd_dcrt, hcr_Sdd_dcrt, partMCMC_Sdd_dcrt)
S_dd_FP <- metric_Sdd("FP", datmat, tabur_Sdd, hcr_Sdd, partMCMC_Sdd, tabur_Sdd_like, hcr_Sdd_like, structMCMC_Sdd, tabur_Sdd_dcrt, hcr_Sdd_dcrt, partMCMC_Sdd_dcrt)
S_dd_FN <- metric_Sdd("FN", datmat, tabur_Sdd, hcr_Sdd, partMCMC_Sdd, tabur_Sdd_like, hcr_Sdd_like, structMCMC_Sdd, tabur_Sdd_dcrt, hcr_Sdd_dcrt, partMCMC_Sdd_dcrt)

S_dd_SHD$F1 <- F1score(S_dd_TP$metric, S_dd_FP$metric, S_dd_FN$metric)


plot_metric(S_dd_SHD, "SHD") 

#################################################################################
S_cc_TPR$Scenario <- S_cc_SHD$Scenario <- "S_cc"
S_cd_TPR$Scenario <- S_cd_SHD$Scenario <- "S_cd"
S_dc_TPR$Scenario <- S_dc_SHD$Scenario <- "S_dc"
S_dd_TPR$Scenario <- S_dd_SHD$Scenario <- "S_dd"

rtable_SHD <- rbind(S_cc_SHD, S_cd_SHD, S_dc_SHD, S_dd_SHD)
rtable_TPR <- rbind(S_cc_TPR, S_cd_TPR, S_dc_TPR, S_dd_TPR)

rtable_SHD = rtable_SHD %>% mutate( Scenario_new = factor(Scenario,
                                                          levels = c("S_cc", "S_cd", "S_dc", "S_dd"),
                                                          labels = c(expression(S[cc]),expression(S[cd]),expression(S[dc]), expression(S[dd]))))

rtable_TPR = rtable_TPR %>% mutate( Scenario_new = factor(Scenario,
                                                          levels = c("S_cc", "S_cd", "S_dc", "S_dd"),
                                                          labels = c(expression(S[cc]),expression(S[cd]),expression(S[dc]), expression(S[dd]))))

# ggplot(rtable_SHD, aes(x=Strategy, y=metric)) +  
#   geom_boxplot(aes(colour = Method)) + 
#   facet_grid( Scenario_new ~ Method , labeller = labeller(Scenario_new = label_parsed)) +
#   stat_summary(fun = mean, color = "darkblue", position = position_dodge(0.75),
#                geom = "point", shape = 18, size = 2,
#                show.legend = FALSE) +
#   theme_light() + labs(y= "SHD") + 
#   theme(text = element_text(size=9), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position ='None')
# 
# ggplot(rtable_TPR, aes(x=Strategy, y=metric)) +  
#   geom_boxplot(aes(colour = Method)) + 
#   facet_grid( Scenario_new ~ Method , labeller = labeller(Scenario_new = label_parsed)) +
#   stat_summary(fun = mean, color = "darkblue", position = position_dodge(0.75),
#                geom = "point", shape = 18, size = 2,
#                show.legend = FALSE) +
#   theme_light() + labs(y= "TPR") + 
#   theme(text = element_text(size=9), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position ='None')

# Pivot the data to long format
rtable_long <- rtable_SHD %>% 
  filter(Method == "MCMC") %>% 
  gather(key = "variable", value = "value", metric, F1)

# Add Color column to rtable_long based on Strategy
rtable_long <- rtable_long %>%
  mutate(Color = case_when(
    grepl("DISC", Strategy) ~ "DISC",
    grepl("RAG", Strategy) ~ "RAG",
    TRUE ~ "CLG"  # Assumes any non-DISC or RAG strategy is CLG
  ))

# Define custom colors with a new palette
custom_colors <- c("DISC" = brewer.pal(3, "Dark2")[2],  # Orange
                   "RAG" = brewer.pal(3, "Dark2")[1],   # Green
                   "CLG" = brewer.pal(3, "Dark2")[3])  # Purple

# Assuming there are four scenarios to fit into one row, adjust ncol as needed
ncol_scenario <- 4  # Adjust based on the number of unique scenarios you have

# Plot for F1 (with x-axis labels removed)
p_f1 <- ggplot(rtable_long %>% filter(variable == "F1"), 
  aes(x = Strategy, y = value, color = Color)) +  
  geom_boxplot() + 
  facet_wrap(~ Scenario_new, scales = "free_x", nrow = 1, ncol = ncol_scenario, labeller = labeller(Scenario_new = label_parsed)) +
  stat_summary(fun = mean, color = "#666666", position = position_dodge(0.75),
               geom = "point", shape = 18, size = 2,
               show.legend = FALSE) +
  scale_color_manual(values = custom_colors) +
  scale_y_continuous(limits = c(0.3, 1)) +
  theme_light() + 
  labs(y = "F1", x = NULL) + 
  theme(text = element_text(size=11), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        legend.position = 'None',
        strip.text = element_text(size = 14))  # Increase facet title text size


# Plot for SHD without y-axis limits and with no panel titles
p_shd <- ggplot(rtable_long %>% filter(variable == "metric"), 
  aes(x = Strategy, y = value, color = Color)) +  
  geom_boxplot() + 
  facet_wrap(~ Scenario_new, scales = "free_x", nrow = 1, ncol = ncol_scenario, labeller = labeller(Scenario_new = label_parsed)) +
  stat_summary(fun = mean, color = "#666666", position = position_dodge(0.75),
               geom = "point", shape = 18, size = 2,
               show.legend = FALSE) +
  scale_color_manual(values = custom_colors) +
  theme_light() + 
  labs(y = "SHD", x = "Strategy") + 
  theme(text = element_text(size=11), 
        axis.text.x = element_text(angle = 45, hjust=1), 
        legend.position = 'None',
        strip.text = element_blank())  # This removes the facet titles

# Combine the two plots with patchwork
combined_plot <- p_f1 / p_shd + plot_layout(heights = c(1, 1))

# Print the combined plot
print(combined_plot)

ggsave("MCMC_10_nodes_Plot.pdf", combined_plot, device="pdf",width = 11, height = 7.5)
# save.image(file = "10nodes/Results.Rdata")
