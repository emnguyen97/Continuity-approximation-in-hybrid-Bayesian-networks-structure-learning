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
library(gRain)
library(dplyr)
library(ggplot2)

# Load the necessary functions
source('10nodes/MethodsUtil.R')
source('10nodes/DataGenUtil.R')

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

### RAG
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

partMCMC <- foreach(i = seq_len(100), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(datasets[[i]])
  result$result <- runpartMCMC(data, 'bge')
  return(result)
}
stopCluster(cl)

options(digits = 4)
get_result(partMCMC,datmat)

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
  
  partMCMC_interval_2b <- foreach(i = seq_len(100), .combine = 'c') %dopar% {
    result <- multiResultClass()
    data <- as.data.frame(datasets_interval_2b[[i]])
    result$result <- runpartMCMC(data, 'bdecat')
    return(result)
  }
  
  partMCMC_interval_4b <- foreach(i = seq_len(100), .combine = 'c') %dopar% {
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
  
  partMCMC_clust_2b <- foreach(i = seq_len(100), .combine = 'c') %dopar% {
    result <- multiResultClass()
    data <- as.data.frame(datasets_clust_2b[[i]])
    result$result <- runpartMCMC(data, 'bdecat')
    return(result)
  }
  
  partMCMC_clust_4b <- foreach(i = seq_len(100), .combine = 'c') %dopar% {
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

### Tabulate results

metric_Scc <- function(metric, datmat, tabur, hcr, partMCMC, tabur_interval_2b, hcr_interval_2b, partMCMC_interval_2b, tabur_interval_4b, hcr_interval_4b, partMCMC_interval_4b, tabur_clust_2b, hcr_clust_2b, partMCMC_clust_2b, tabur_clust_4b, hcr_clust_4b, partMCMC_clust_4b){
  RAG <- cbind(TABU = get_metric(tabur,datmat,metric), HC = get_metric(hcr,datmat,metric), MCMC = get_metric(partMCMC,datmat,metric))
  res_long <- tidyr::gather(as.data.frame(RAG), Method, metric, TABU:MCMC, factor_key=TRUE)
  res_long$Strategy <- "RAG (BGe)"
  
  DCRT_inv_2b <- cbind(TABU = get_metric(tabur_interval_2b,datmat,metric), HC = get_metric(hcr_interval_2b,datmat,metric), MCMC = get_metric(partMCMC_interval_2b,datmat,metric))
  DCRT_inv_4b <- cbind(TABU = get_metric(tabur_interval_4b,datmat,metric), HC = get_metric(hcr_interval_4b,datmat,metric), MCMC = get_metric(partMCMC_interval_4b,datmat,metric))
  DCRT_clust_2b <- cbind(TABU = get_metric(tabur_clust_2b,datmat,metric), HC = get_metric(hcr_clust_2b,datmat,metric), MCMC = get_metric(partMCMC_clust_2b,datmat,metric))
  DCRT_clust_4b <- cbind(TABU = get_metric(tabur_clust_4b,datmat,metric), HC = get_metric(hcr_clust_4b,datmat,metric), MCMC = get_metric(partMCMC_clust_4b,datmat,metric))
  DCRT_inv_long_2b <- as.data.frame(tidyr::gather(as.data.frame(DCRT_inv_2b), Method, metric, TABU:MCMC, factor_key=TRUE))
  DCRT_inv_long_4b <- as.data.frame(tidyr::gather(as.data.frame(DCRT_inv_4b), Method, metric, TABU:MCMC, factor_key=TRUE))
  DCRT_clust_long_2b <- as.data.frame(tidyr::gather(as.data.frame(DCRT_clust_2b), Method, metric, TABU:MCMC, factor_key=TRUE))
  DCRT_clust_long_4b <- as.data.frame(tidyr::gather(as.data.frame(DCRT_clust_4b), Method, metric, TABU:MCMC, factor_key=TRUE))
  DCRT_inv_long_2b$Strategy <- "DISC2 - I (BDe)"
  DCRT_clust_long_2b$Strategy <- "DISC2 - C (BDe)"
  DCRT_inv_long_4b$Strategy <- "DISC4 - I (BDe)"
  DCRT_clust_long_4b$Strategy <- "DISC4 - C (BDe)"
  res_long <- rbind(res_long, DCRT_inv_long_2b, DCRT_inv_long_4b, DCRT_clust_long_2b, DCRT_clust_long_4b)
  
  return(res_long)
}

plot_metric <- function(res_long, metric){
  
  P <- ggplot(res_long, aes(x=Strategy, y=metric)) +  
    geom_boxplot(aes(colour = Method)) + 
    facet_wrap(~Method) +
    stat_summary(fun = mean, color = "darkblue", position = position_dodge(0.75),
                 geom = "point", shape = 18, size = 3,
                 show.legend = FALSE) +
    theme_light() + labs(y= "SHD") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position ='None')
  
  return(P)
}

F1score <- function(TP, FP, FN) {
  return( TP / (TP + (1/2) * (FP + FN) ) )
}

S_cc_SHD <- metric_Scc("SHD", datmat, tabur, hcr, partMCMC, tabur_Scc_interval_2b, hcr_Scc_interval_2b, partMCMC_Scc_interval_2b, tabur_Scc_interval_4b, hcr_Scc_interval_4b, partMCMC_Scc_interval_4b, tabur_Scc_clust_2b, hcr_Scc_clust_2b, partMCMC_Scc_clust_2b, tabur_Scc_clust_4b, hcr_Scc_clust_4b, partMCMC_Scc_clust_4b)

S_cc_TPR <- metric_Scc("TPR", datmat, tabur, hcr, partMCMC, tabur_Scc_interval_2b, hcr_Scc_interval_2b, partMCMC_Scc_interval_2b, tabur_Scc_interval_4b, hcr_Scc_interval_4b, partMCMC_Scc_interval_4b, tabur_Scc_clust_2b, hcr_Scc_clust_2b, partMCMC_Scc_clust_2b, tabur_Scc_clust_4b, hcr_Scc_clust_4b, partMCMC_Scc_clust_4b)

S_cc_TP <- metric_Scc("TP", datmat, tabur, hcr, partMCMC, tabur_Scc_interval_2b, hcr_Scc_interval_2b, partMCMC_Scc_interval_2b, tabur_Scc_interval_4b, hcr_Scc_interval_4b, partMCMC_Scc_interval_4b, tabur_Scc_clust_2b, hcr_Scc_clust_2b, partMCMC_Scc_clust_2b, tabur_Scc_clust_4b, hcr_Scc_clust_4b, partMCMC_Scc_clust_4b)

S_cc_FP <- metric_Scc("FP", datmat, tabur, hcr, partMCMC, tabur_Scc_interval_2b, hcr_Scc_interval_2b, partMCMC_Scc_interval_2b, tabur_Scc_interval_4b, hcr_Scc_interval_4b, partMCMC_Scc_interval_4b, tabur_Scc_clust_2b, hcr_Scc_clust_2b, partMCMC_Scc_clust_2b, tabur_Scc_clust_4b, hcr_Scc_clust_4b, partMCMC_Scc_clust_4b)

S_cc_FN <- metric_Scc("FN", datmat, tabur, hcr, partMCMC, tabur_Scc_interval_2b, hcr_Scc_interval_2b, partMCMC_Scc_interval_2b, tabur_Scc_interval_4b, hcr_Scc_interval_4b, partMCMC_Scc_interval_4b, tabur_Scc_clust_2b, hcr_Scc_clust_2b, partMCMC_Scc_clust_2b, tabur_Scc_clust_4b, hcr_Scc_clust_4b, partMCMC_Scc_clust_4b)

S_cc_SHD$F1 <- F1score(S_cc_TP$metric, S_cc_FP$metric, S_cc_FN$metric)

plot_metric(S_cc_SHD, "SHD") 
plot_metric(S_cc_TPR, "TPR") 

# plot_metric_Scc <- function(metric, datmat, tabur, hcr, partMCMC, tabur_interval_2b, hcr_interval_2b, partMCMC_interval_2b, tabur_interval_4b, hcr_interval_4b, partMCMC_interval_4b, tabur_clust_2b, hcr_clust_2b, partMCMC_clust_2b, tabur_clust_4b, hcr_clust_4b, partMCMC_clust_4b){
#   RAG <- cbind(TABU = get_metric(tabur,datmat,metric), HC = get_metric(hcr,datmat,metric), MCMC = get_metric(partMCMC,datmat,metric))
#   res_long <- tidyr::gather(as.data.frame(RAG), Method, metric, TABU:MCMC, factor_key=TRUE)
#   res_long$Strategy <- "RAG (BGe)"
#   
#   DCRT_inv_2b <- cbind(TABU = get_metric(tabur_interval_2b,datmat,metric), HC = get_metric(hcr_interval_2b,datmat,metric), MCMC = get_metric(partMCMC_interval_2b,datmat,metric))
#   DCRT_inv_4b <- cbind(TABU = get_metric(tabur_interval_4b,datmat,metric), HC = get_metric(hcr_interval_4b,datmat,metric), MCMC = get_metric(partMCMC_interval_4b,datmat,metric))
#   DCRT_clust_2b <- cbind(TABU = get_metric(tabur_clust_2b,datmat,metric), HC = get_metric(hcr_clust_2b,datmat,metric), MCMC = get_metric(partMCMC_clust_2b,datmat,metric))
#   DCRT_clust_4b <- cbind(TABU = get_metric(tabur_clust_4b,datmat,metric), HC = get_metric(hcr_clust_4b,datmat,metric), MCMC = get_metric(partMCMC_clust_4b,datmat,metric))
#   DCRT_inv_long_2b <- as.data.frame(tidyr::gather(as.data.frame(DCRT_inv_2b), Method, metric, TABU:MCMC, factor_key=TRUE))
#   DCRT_inv_long_4b <- as.data.frame(tidyr::gather(as.data.frame(DCRT_inv_4b), Method, metric, TABU:MCMC, factor_key=TRUE))
#   DCRT_clust_long_2b <- as.data.frame(tidyr::gather(as.data.frame(DCRT_clust_2b), Method, metric, TABU:MCMC, factor_key=TRUE))
#   DCRT_clust_long_4b <- as.data.frame(tidyr::gather(as.data.frame(DCRT_clust_4b), Method, metric, TABU:MCMC, factor_key=TRUE))
#   DCRT_inv_long_2b$Strategy <- "DISC2 - INTV (BDe)"
#   DCRT_clust_long_2b$Strategy <- "DISC2 - CLST (BDe)"
#   DCRT_inv_long_4b$Strategy <- "DISC4 - INTV (BDe)"
#   DCRT_clust_long_4b$Strategy <- "DISC4 - CLST (BDe)"
#   res_long <- rbind(res_long, DCRT_inv_long_2b, DCRT_inv_long_4b, DCRT_clust_long_2b, DCRT_clust_long_4b)
#   
#   P <- ggplot(res_long, aes(x=Strategy, y=metric)) +  
#     geom_boxplot(aes(colour = Method)) + 
#     facet_wrap(~Method) +
#     stat_summary(fun.y = mean, color = "darkblue", position = position_dodge(0.75),
#                  geom = "point", shape = 18, size = 3,
#                  show.legend = FALSE) +
#     theme_light() 
#   return(P)
# }
# 
# # SHD
# p <- plot_metric_Scc("SHD", datmat, tabur, hcr, partMCMC, tabur_Scc_interval_2b, hcr_Scc_interval_2b, partMCMC_Scc_interval_2b, tabur_Scc_interval_4b, hcr_Scc_interval_4b, partMCMC_Scc_interval_4b, tabur_Scc_clust_2b, hcr_Scc_clust_2b, partMCMC_Scc_clust_2b, tabur_Scc_clust_4b, hcr_Scc_clust_4b, partMCMC_Scc_clust_4b) 
# #p + labs(y = "SHD") + theme(legend.position = c(0.92, 0.86))
# p + labs(y= "SHD") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position ='None')

#ggsave('10nodes/Scc_SHD.pdf', plot = p, width = 6, height = 4)

# # TPR
# p2 <- plot_metric_Scc("TPR", datmat, tabur, hcr, partMCMC, tabur_Scc_interval_2b, hcr_Scc_interval_2b, partMCMC_Scc_interval_2b, tabur_Scc_interval_4b, hcr_Scc_interval_4b, partMCMC_Scc_interval_4b, tabur_Scc_clust_2b, hcr_Scc_clust_2b, partMCMC_Scc_clust_2b, tabur_Scc_clust_4b, hcr_Scc_clust_4b, partMCMC_Scc_clust_4b)
# #p2 + labs(y= "TPR") + theme(legend.position = c(0.92, 0.14))
# p2 + labs(y= "TPR") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position ='None')

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

partMCMC_Scd <- foreach(i = seq_len(100), .combine = 'c') %dopar% {
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

### Discretise
dcrtnodes_Scd <- c("A","C","G","H")
contnodes_Scd <- setdiff(colnames(datmat),dcrtnodes_Scd)

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


### Tabulate results

metric_comb <- function(metric, datmat, tabur, hcr, partMCMC, tabur_cg, hcr_cg, tabur_interval_2b, hcr_interval_2b, partMCMC_interval_2b, tabur_interval_4b, hcr_interval_4b, partMCMC_interval_4b, tabur_clust_2b, hcr_clust_2b, partMCMC_clust_2b, tabur_clust_4b, hcr_clust_4b, partMCMC_clust_4b){
  
  RAG <- cbind(TABU = get_metric(tabur,datmat,metric), HC = get_metric(hcr,datmat,metric), MCMC = get_metric(partMCMC,datmat,metric))
  res_long <- tidyr::gather(as.data.frame(RAG), Method, metric, TABU:MCMC, factor_key=TRUE)
  res_long$Strategy <- "RAG (BGe)"
  
  CLG <- cbind(TABU = get_metric(tabur_cg,datmat,metric), HC = get_metric(hcr_cg,datmat,metric))
  CLG_long <- tidyr::gather(as.data.frame(CLG), Method, metric, TABU:HC, factor_key=TRUE)
  CLG_long$Strategy <- "CLG"
  res_long <- rbind(res_long, CLG_long)
  
  DCRT_inv_2b <- cbind(TABU = get_metric(tabur_interval_2b,datmat,metric), HC = get_metric(hcr_interval_2b,datmat,metric), MCMC = get_metric(partMCMC_interval_2b,datmat,metric))
  DCRT_inv_4b <- cbind(TABU = get_metric(tabur_interval_4b,datmat,metric), HC = get_metric(hcr_interval_4b,datmat,metric), MCMC = get_metric(partMCMC_interval_4b,datmat,metric))
  DCRT_clust_2b <- cbind(TABU = get_metric(tabur_clust_2b,datmat,metric), HC = get_metric(hcr_clust_2b,datmat,metric), MCMC = get_metric(partMCMC_clust_2b,datmat,metric))
  DCRT_clust_4b <- cbind(TABU = get_metric(tabur_clust_4b,datmat,metric), HC = get_metric(hcr_clust_4b,datmat,metric), MCMC = get_metric(partMCMC_clust_4b,datmat,metric))
  DCRT_inv_long_2b <- as.data.frame(tidyr::gather(as.data.frame(DCRT_inv_2b), Method, metric, TABU:MCMC, factor_key=TRUE))
  DCRT_inv_long_4b <- as.data.frame(tidyr::gather(as.data.frame(DCRT_inv_4b), Method, metric, TABU:MCMC, factor_key=TRUE))
  DCRT_clust_long_2b <- as.data.frame(tidyr::gather(as.data.frame(DCRT_clust_2b), Method, metric, TABU:MCMC, factor_key=TRUE))
  DCRT_clust_long_4b <- as.data.frame(tidyr::gather(as.data.frame(DCRT_clust_4b), Method, metric, TABU:MCMC, factor_key=TRUE))
  DCRT_inv_long_2b$Strategy <- "DISC2 - I (BDe)"
  DCRT_clust_long_2b$Strategy <- "DISC2 - C (BDe)"
  DCRT_inv_long_4b$Strategy <- "DISC4 - I (BDe)"
  DCRT_clust_long_4b$Strategy <- "DISC4 - C (BDe)"
  res_long <- rbind(res_long, DCRT_inv_long_2b, DCRT_inv_long_4b, DCRT_clust_long_2b, DCRT_clust_long_4b)
  
  return(res_long)
}

S_cd_SHD <- metric_comb("SHD", datmat, tabur_Scd, hcr_Scd, partMCMC_Scd, tabur_Scd_cg, hcr_Scd_cg, tabur_Scd_interval_2b, hcr_Scd_interval_2b, partMCMC_Scd_interval_2b, tabur_Scd_interval_4b, hcr_Scd_interval_4b, partMCMC_Scd_interval_4b, tabur_Scd_clust_2b, hcr_Scd_clust_2b, partMCMC_Scd_clust_2b, tabur_Scd_clust_4b, hcr_Scd_clust_4b, partMCMC_Scd_clust_4b)

S_cd_TP <- metric_comb("TP", datmat, tabur_Scd, hcr_Scd, partMCMC_Scd, tabur_Scd_cg, hcr_Scd_cg, tabur_Scd_interval_2b, hcr_Scd_interval_2b, partMCMC_Scd_interval_2b, tabur_Scd_interval_4b, hcr_Scd_interval_4b, partMCMC_Scd_interval_4b, tabur_Scd_clust_2b, hcr_Scd_clust_2b, partMCMC_Scd_clust_2b, tabur_Scd_clust_4b, hcr_Scd_clust_4b, partMCMC_Scd_clust_4b)

S_cd_FP <- metric_comb("FP", datmat, tabur_Scd, hcr_Scd, partMCMC_Scd, tabur_Scd_cg, hcr_Scd_cg, tabur_Scd_interval_2b, hcr_Scd_interval_2b, partMCMC_Scd_interval_2b, tabur_Scd_interval_4b, hcr_Scd_interval_4b, partMCMC_Scd_interval_4b, tabur_Scd_clust_2b, hcr_Scd_clust_2b, partMCMC_Scd_clust_2b, tabur_Scd_clust_4b, hcr_Scd_clust_4b, partMCMC_Scd_clust_4b)

S_cd_FN <- metric_comb("FN", datmat, tabur_Scd, hcr_Scd, partMCMC_Scd, tabur_Scd_cg, hcr_Scd_cg, tabur_Scd_interval_2b, hcr_Scd_interval_2b, partMCMC_Scd_interval_2b, tabur_Scd_interval_4b, hcr_Scd_interval_4b, partMCMC_Scd_interval_4b, tabur_Scd_clust_2b, hcr_Scd_clust_2b, partMCMC_Scd_clust_2b, tabur_Scd_clust_4b, hcr_Scd_clust_4b, partMCMC_Scd_clust_4b)

S_cd_TPR <- metric_comb("TPR", datmat, tabur_Scd, hcr_Scd, partMCMC_Scd, tabur_Scd_cg, hcr_Scd_cg, tabur_Scd_interval_2b, hcr_Scd_interval_2b, partMCMC_Scd_interval_2b, tabur_Scd_interval_4b, hcr_Scd_interval_4b, partMCMC_Scd_interval_4b, tabur_Scd_clust_2b, hcr_Scd_clust_2b, partMCMC_Scd_clust_2b, tabur_Scd_clust_4b, hcr_Scd_clust_4b, partMCMC_Scd_clust_4b)

S_cd_SHD$F1 <- F1score(S_cd_TP$metric, S_cd_FP$metric, S_cd_FN$metric)

plot_metric(S_cd_SHD, "SHD") 
plot_metric(S_cd_TPR, "TPR") 


# plot_metric_comb <- function(metric, datmat, tabur, hcr, partMCMC, tabur_cg, hcr_cg, tabur_interval_2b, hcr_interval_2b, partMCMC_interval_2b, tabur_interval_4b, hcr_interval_4b, partMCMC_interval_4b, tabur_clust_2b, hcr_clust_2b, partMCMC_clust_2b, tabur_clust_4b, hcr_clust_4b, partMCMC_clust_4b){
#   
#   RAG <- cbind(TABU = get_metric(tabur,datmat,metric), HC = get_metric(hcr,datmat,metric), MCMC = get_metric(partMCMC,datmat,metric))
#   res_long <- tidyr::gather(as.data.frame(RAG), Method, metric, TABU:MCMC, factor_key=TRUE)
#   res_long$Strategy <- "RAG (BGe)"
#   
#   CLG <- cbind(TABU = get_metric(tabur_cg,datmat,metric), HC = get_metric(hcr_cg,datmat,metric))
#   CLG_long <- tidyr::gather(as.data.frame(CLG), Method, metric, TABU:HC, factor_key=TRUE)
#   CLG_long$Strategy <- "CLG"
#   res_long <- rbind(res_long, CLG_long)
#   
#   DCRT_inv_2b <- cbind(TABU = get_metric(tabur_interval_2b,datmat,metric), HC = get_metric(hcr_interval_2b,datmat,metric), MCMC = get_metric(partMCMC_interval_2b,datmat,metric))
#   DCRT_inv_4b <- cbind(TABU = get_metric(tabur_interval_4b,datmat,metric), HC = get_metric(hcr_interval_4b,datmat,metric), MCMC = get_metric(partMCMC_interval_4b,datmat,metric))
#   DCRT_clust_2b <- cbind(TABU = get_metric(tabur_clust_2b,datmat,metric), HC = get_metric(hcr_clust_2b,datmat,metric), MCMC = get_metric(partMCMC_clust_2b,datmat,metric))
#   DCRT_clust_4b <- cbind(TABU = get_metric(tabur_clust_4b,datmat,metric), HC = get_metric(hcr_clust_4b,datmat,metric), MCMC = get_metric(partMCMC_clust_4b,datmat,metric))
#   DCRT_inv_long_2b <- as.data.frame(tidyr::gather(as.data.frame(DCRT_inv_2b), Method, metric, TABU:MCMC, factor_key=TRUE))
#   DCRT_inv_long_4b <- as.data.frame(tidyr::gather(as.data.frame(DCRT_inv_4b), Method, metric, TABU:MCMC, factor_key=TRUE))
#   DCRT_clust_long_2b <- as.data.frame(tidyr::gather(as.data.frame(DCRT_clust_2b), Method, metric, TABU:MCMC, factor_key=TRUE))
#   DCRT_clust_long_4b <- as.data.frame(tidyr::gather(as.data.frame(DCRT_clust_4b), Method, metric, TABU:MCMC, factor_key=TRUE))
#   DCRT_inv_long_2b$Strategy <- "DISC2 - INTV (BDe)"
#   DCRT_clust_long_2b$Strategy <- "DISC2 - CLST (BDe)"
#   DCRT_inv_long_4b$Strategy <- "DISC4 - INTV (BDe)"
#   DCRT_clust_long_4b$Strategy <- "DISC4 - CLST (BDe)"
#   res_long <- rbind(res_long, DCRT_inv_long_2b, DCRT_inv_long_4b, DCRT_clust_long_2b, DCRT_clust_long_4b)
#   
#   P <- ggplot(res_long, aes(x=Strategy, y=metric)) +  
#     geom_boxplot(aes(colour = Method)) + 
#     facet_wrap(~Method) +
#     stat_summary(fun.y = mean, color = "darkblue", position = position_dodge(0.75),
#                  geom = "point", shape = 18, size = 3,
#                  show.legend = FALSE) +
#     theme_light() 
#   return(P)
# }
# 
# # SHD
# p <- plot_metric_comb("SHD", datmat, tabur_Scd, hcr_Scd, partMCMC_Scd, tabur_Scd_cg, hcr_Scd_cg, tabur_Scd_interval_2b, hcr_Scd_interval_2b, partMCMC_Scd_interval_2b, tabur_Scd_interval_4b, hcr_Scd_interval_4b, partMCMC_Scd_interval_4b, tabur_Scd_clust_2b, hcr_Scd_clust_2b, partMCMC_Scd_clust_2b, tabur_Scd_clust_4b, hcr_Scd_clust_4b, partMCMC_Scd_clust_4b)
# #p + labs(y= "SHD") + theme(legend.position = c(0.92, 0.86))
# p + labs(y= "SHD") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position ='None')
# 
# # TPR
# p2 <- plot_metric_comb("TPR", datmat, tabur_Scd, hcr_Scd, partMCMC_Scd, tabur_Scd_cg, hcr_Scd_cg, tabur_Scd_interval_2b, hcr_Scd_interval_2b, partMCMC_Scd_interval_2b, tabur_Scd_interval_4b, hcr_Scd_interval_4b, partMCMC_Scd_interval_4b, tabur_Scd_clust_2b, hcr_Scd_clust_2b, partMCMC_Scd_clust_2b, tabur_Scd_clust_4b, hcr_Scd_clust_4b, partMCMC_Scd_clust_4b)
# #p2 + labs(y= "TPR") + theme(legend.position = c(0.92, 0.14))
# p2 + labs(y= "TPR") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position ='None')

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


### RAG
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

partMCMC_Sdc <- foreach(i = seq_len(100), .combine = 'c') %dopar% {
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

### Tabulate results
S_dc_SHD <- metric_comb("SHD", datmat, tabur_Sdc, hcr_Sdc, partMCMC_Sdc, tabur_Sdc_cg, hcr_Sdc_cg, tabur_Sdc_interval_2b, hcr_Sdc_interval_2b, partMCMC_Sdc_interval_2b, tabur_Sdc_interval_4b, hcr_Sdc_interval_4b, partMCMC_Sdc_interval_4b, tabur_Sdc_clust_2b, hcr_Sdc_clust_2b, partMCMC_Sdc_clust_2b, tabur_Sdc_clust_4b, hcr_Sdc_clust_4b, partMCMC_Sdc_clust_4b)

S_dc_TPR <- metric_comb("TPR", datmat, tabur_Sdc, hcr_Sdc, partMCMC_Sdc, tabur_Sdc_cg, hcr_Sdc_cg, tabur_Sdc_interval_2b, hcr_Sdc_interval_2b, partMCMC_Sdc_interval_2b, tabur_Sdc_interval_4b, hcr_Sdc_interval_4b, partMCMC_Sdc_interval_4b, tabur_Sdc_clust_2b, hcr_Sdc_clust_2b, partMCMC_Sdc_clust_2b, tabur_Sdc_clust_4b, hcr_Sdc_clust_4b, partMCMC_Sdc_clust_4b)

S_dc_TP <- metric_comb("TP", datmat, tabur_Sdc, hcr_Sdc, partMCMC_Sdc, tabur_Sdc_cg, hcr_Sdc_cg, tabur_Sdc_interval_2b, hcr_Sdc_interval_2b, partMCMC_Sdc_interval_2b, tabur_Sdc_interval_4b, hcr_Sdc_interval_4b, partMCMC_Sdc_interval_4b, tabur_Sdc_clust_2b, hcr_Sdc_clust_2b, partMCMC_Sdc_clust_2b, tabur_Sdc_clust_4b, hcr_Sdc_clust_4b, partMCMC_Sdc_clust_4b)

S_dc_FP <- metric_comb("FP", datmat, tabur_Sdc, hcr_Sdc, partMCMC_Sdc, tabur_Sdc_cg, hcr_Sdc_cg, tabur_Sdc_interval_2b, hcr_Sdc_interval_2b, partMCMC_Sdc_interval_2b, tabur_Sdc_interval_4b, hcr_Sdc_interval_4b, partMCMC_Sdc_interval_4b, tabur_Sdc_clust_2b, hcr_Sdc_clust_2b, partMCMC_Sdc_clust_2b, tabur_Sdc_clust_4b, hcr_Sdc_clust_4b, partMCMC_Sdc_clust_4b)

S_dc_FN <- metric_comb("FN", datmat, tabur_Sdc, hcr_Sdc, partMCMC_Sdc, tabur_Sdc_cg, hcr_Sdc_cg, tabur_Sdc_interval_2b, hcr_Sdc_interval_2b, partMCMC_Sdc_interval_2b, tabur_Sdc_interval_4b, hcr_Sdc_interval_4b, partMCMC_Sdc_interval_4b, tabur_Sdc_clust_2b, hcr_Sdc_clust_2b, partMCMC_Sdc_clust_2b, tabur_Sdc_clust_4b, hcr_Sdc_clust_4b, partMCMC_Sdc_clust_4b)


S_dc_SHD$F1 <- F1score(S_dc_TP$metric, S_dc_FP$metric, S_dc_FN$metric)


plot_metric(S_dc_SHD, "SHD")
plot_metric(S_dc_TPR, "TPR")

# # SHD
# p <- plot_metric_comb("SHD", datmat, tabur_Sdc, hcr_Sdc, partMCMC_Sdc, tabur_Sdc_cg, hcr_Sdc_cg, tabur_Sdc_interval_2b, hcr_Sdc_interval_2b, partMCMC_Sdc_interval_2b, tabur_Sdc_interval_4b, hcr_Sdc_interval_4b, partMCMC_Sdc_interval_4b, tabur_Sdc_clust_2b, hcr_Sdc_clust_2b, partMCMC_Sdc_clust_2b, tabur_Sdc_clust_4b, hcr_Sdc_clust_4b, partMCMC_Sdc_clust_4b)
# #p + labs(y= "SHD") + theme(legend.position = c(0.06, 0.86))
# p + labs(y= "SHD") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position ='None')
# 
# # TPR
# p2 <- plot_metric_comb("TPR", datmat, tabur_Sdc, hcr_Sdc, partMCMC_Sdc, tabur_Sdc_cg, hcr_Sdc_cg, tabur_Sdc_interval_2b, hcr_Sdc_interval_2b, partMCMC_Sdc_interval_2b, tabur_Sdc_interval_4b, hcr_Sdc_interval_4b, partMCMC_Sdc_interval_4b, tabur_Sdc_clust_2b, hcr_Sdc_clust_2b, partMCMC_Sdc_clust_2b, tabur_Sdc_clust_4b, hcr_Sdc_clust_4b, partMCMC_Sdc_clust_4b)
# #p2 + labs(y= "TPR") + theme(legend.position = c(0.06, 0.86))
# p2 + labs(y= "TPR") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position ='None')

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

partMCMC_Sdd <- foreach(i = seq_len(100), .combine = 'c') %dopar% {
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

partMCMC_Sdd_dcrt <- foreach(i = seq_len(100), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(datasets_Sdd[[i]])
  result$result <- runpartMCMC(data, 'bdecat')
  return(result)
}
stopCluster(cl)

options(digits = 4)
get_result(partMCMC_Sdd_dcrt,datmat)

### Tabulate results

metric_Sdd <- function(metric, datmat, tabur, hcr, partMCMC, tabur_dcrt, hcr_dcrt, partMCMC_dcrt){
  RAG <- cbind(TABU = get_metric(tabur,datmat,metric), HC = get_metric(hcr,datmat,metric), MCMC = get_metric(partMCMC,datmat,metric))
  res_long <- tidyr::gather(as.data.frame(RAG), Method, metric, TABU:MCMC, factor_key=TRUE)
  res_long$Strategy <- "RAG (BGe)"
  
  DISC <- cbind(TABU = get_metric(tabur_dcrt,datmat,metric), HC = get_metric(hcr_dcrt,datmat,metric), MCMC = get_metric(partMCMC_dcrt,datmat,metric))
  DISC_long <- tidyr::gather(as.data.frame(DISC), Method, metric, TABU:MCMC, factor_key=TRUE)
  DISC_long$Strategy <- "DISC (BDe)"
  res_long <- rbind(res_long, DISC_long)

  return(res_long)
}

S_dd_SHD <- metric_Sdd("SHD", datmat, tabur_Sdd, hcr_Sdd, partMCMC_Sdd, tabur_Sdd_dcrt, hcr_Sdd_dcrt, partMCMC_Sdd_dcrt)
S_dd_TPR <- metric_Sdd("TPR", datmat, tabur_Sdd, hcr_Sdd, partMCMC_Sdd, tabur_Sdd_dcrt, hcr_Sdd_dcrt, partMCMC_Sdd_dcrt)

S_dd_TP <- metric_Sdd("TP", datmat, tabur_Sdd, hcr_Sdd, partMCMC_Sdd, tabur_Sdd_dcrt, hcr_Sdd_dcrt, partMCMC_Sdd_dcrt)
S_dd_FP <- metric_Sdd("FP", datmat, tabur_Sdd, hcr_Sdd, partMCMC_Sdd, tabur_Sdd_dcrt, hcr_Sdd_dcrt, partMCMC_Sdd_dcrt)
S_dd_FN <- metric_Sdd("FN", datmat, tabur_Sdd, hcr_Sdd, partMCMC_Sdd, tabur_Sdd_dcrt, hcr_Sdd_dcrt, partMCMC_Sdd_dcrt)

S_dd_SHD$F1 <- F1score(S_dd_TP$metric, S_dd_FP$metric, S_dd_FN$metric)


plot_metric(S_dd_SHD, "SHD") 

# 
# 
# plot_metric_Sdd <- function(metric, datmat, tabur, hcr, partMCMC, tabur_dcrt, hcr_dcrt, partMCMC_dcrt){
#   RAG <- cbind(TABU = get_metric(tabur,datmat,metric), HC = get_metric(hcr,datmat,metric), MCMC = get_metric(partMCMC,datmat,metric))
#   res_long <- tidyr::gather(as.data.frame(RAG), Method, metric, TABU:MCMC, factor_key=TRUE)
#   res_long$Strategy <- "RAG (BGe)"
#   
#   DISC <- cbind(TABU = get_metric(tabur_dcrt,datmat,metric), HC = get_metric(hcr_dcrt,datmat,metric), MCMC = get_metric(partMCMC_dcrt,datmat,metric))
#   DISC_long <- tidyr::gather(as.data.frame(DISC), Method, metric, TABU:MCMC, factor_key=TRUE)
#   DISC_long$Strategy <- "DISC (BDe)"
#   res_long <- rbind(res_long, DISC_long)
#   
#   P <- ggplot(res_long, aes(x=Strategy, y=metric)) +  
#     geom_boxplot(aes(colour = Method)) + 
#     facet_wrap(~Method) +
#     stat_summary(fun.y = mean, color = "darkblue", position = position_dodge(0.75),
#                  geom = "point", shape = 18, size = 3,
#                  show.legend = FALSE) +
#     theme_light() 
#   return(P)
# }
# 
# # SHD
# p <- plot_metric_Sdd("SHD", datmat, tabur_Sdd, hcr_Sdd, partMCMC_Sdd, tabur_Sdd_dcrt, hcr_Sdd_dcrt, partMCMC_Sdd_dcrt)
# #p + labs(y= "SHD") + theme(legend.position = c(0.92, 0.86))
# p + labs(y= "SHD") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position ='None')
# 
# # TPR
# p2 <- plot_metric_Sdd("TPR", datmat, tabur_Sdd, hcr_Sdd, partMCMC_Sdd, tabur_Sdd_dcrt, hcr_Sdd_dcrt, partMCMC_Sdd_dcrt)
# #p2 + labs(y= "TPR") + theme(legend.position = c(0.92, 0.14))
# p2 + labs(y= "TPR") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position ='None')

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

ggplot(rtable_SHD, aes(x=Strategy, y=metric)) +  
  geom_boxplot(aes(colour = Method)) + 
  facet_grid( Scenario_new ~ Method , labeller = labeller(Scenario_new = label_parsed)) +
  stat_summary(fun = mean, color = "darkblue", position = position_dodge(0.75),
               geom = "point", shape = 18, size = 2,
               show.legend = FALSE) +
  theme_light() + labs(y= "SHD") + 
  theme(text = element_text(size=9), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position ='None')

ggplot(rtable_TPR, aes(x=Strategy, y=metric)) +  
  geom_boxplot(aes(colour = Method)) + 
  facet_grid( Scenario_new ~ Method , labeller = labeller(Scenario_new = label_parsed)) +
  stat_summary(fun = mean, color = "darkblue", position = position_dodge(0.75),
               geom = "point", shape = 18, size = 2,
               show.legend = FALSE) +
  theme_light() + labs(y= "TPR") + 
  theme(text = element_text(size=9), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position ='None')

save.image(file = "10nodes/13-01Results.Rdata")
