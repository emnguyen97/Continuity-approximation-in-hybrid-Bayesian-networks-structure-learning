

# Implement Tabu search (TABU) greedy search from bnlearn and return adjacency matrix
runtabu <- function(data, datmat, score){
  rtabu <- bnlearn::tabu(data, score=score)
  tabumat <- as.matrix(as_adjacency_matrix(graph_from_edgelist(rtabu$arcs)))
  tabumat <- as.data.frame(tabumat)[,order(colnames(as.data.frame(tabumat)))]
  tabumat <- tabumat[order(rownames(tabumat)),]
  tabumat <- as.matrix(tabumat)
  
  # If the dimension of the estimated matrix does not match the number of nodes, fill in these empty links
  nvar <- ncol(data)
  num <- ncol(tabumat)
  if (num < nvar){
    tbadded <- setdiff(colnames(datmat),colnames(tabumat))
    n <- length(tbadded)
    fillc <- matrix(0, nvar-n, n)
    fillr <- matrix(0, n, nvar)
    tabumat <- cbind(tabumat,fillc)
    tabumat <- rbind(tabumat,fillr)
    colnames(tabumat)[(num+1):(num+n)] <- rownames(tabumat)[(num+1):(num+n)] <- tbadded
    tabumat <- tabumat[LETTERS[1:nvar],LETTERS[1:nvar]]
  }
  
  return(tabumat)
}

# Implement PC algorithm from bnlearn and return adjacency matrix
runhc <- function(data, datmat, score){
  rhc <- bnlearn::hc(as.data.frame(data), score=score)
  hcmat <- as.matrix(as_adjacency_matrix(graph_from_edgelist(rhc$arcs)))
  hcmat <- as.data.frame(hcmat)[,order(colnames(as.data.frame(hcmat)))]
  hcmat <- hcmat[order(rownames(hcmat)),]
  hcmat <- as.matrix(hcmat)
  
  # If the dimension of the estimated matrix does not match the number of nodes, fill in these empty links
  nvar <- ncol(data)
  num <- ncol(hcmat)
  if (num < nvar){
    tbadded <- setdiff(colnames(datmat),colnames(hcmat))
    n <- length(tbadded)
    fillc <- matrix(0, nvar-n, n)
    fillr <- matrix(0, n, nvar)
    hcmat <- cbind(hcmat,fillc)
    hcmat <- rbind(hcmat,fillr)
    colnames(hcmat)[(num+1):(num+n)] <- rownames(hcmat)[(num+1):(num+n)] <- tbadded
    hcmat <- hcmat[LETTERS[1:nvar],LETTERS[1:nvar]]
  }
  
  return(hcmat)
}

# Implement PMCMC from BiDAG and return adjacency matrix

modelpcore<-function(MCMCchain, p, pdag=FALSE, burnin=0.2, DBN=FALSE, nsmall=0, dyn=0, b=0) {
  
  varlabels<-colnames(MCMCchain[[1]])
  n<-nrow(MCMCchain[[1]])
  incidence<-matrix(rep(0, n*n), nrow=n, ncol=n)
  endstep<-length(MCMCchain)
  startstep<-max(as.integer(burnin*endstep),1)
  if (pdag) {
    cpdags<-lapply(MCMCchain[startstep:endstep],dagadj2cpadj)
    incidence[which(as.matrix(Reduce('+', cpdags)/(endstep-startstep+1))>p)]<-1
  } else {
    incidence[which(as.matrix(Reduce('+', MCMCchain[startstep:endstep])/(endstep-startstep+1))>p)]<-1
  }
  colnames(incidence)<-varlabels
  rownames(incidence)<-varlabels
  if(DBN) {
    incidence<-DBNcut(incidence,dyn=dyn,b=b)
    incidence.init<-DBNinit(incidence,dyn=dyn,b=b)
    incidence[1:(dyn+b),1:(dyn+b)]<-incidence.init
  }
  return(incidence)
}

runpartMCMC <- function(data,score){
  
  modelpcore<-function(MCMCchain, p, pdag=FALSE, burnin=0.2, DBN=FALSE, nsmall=0, dyn=0, b=0) {
    
    varlabels<-colnames(MCMCchain[[1]])
    n<-nrow(MCMCchain[[1]])
    incidence<-matrix(rep(0, n*n), nrow=n, ncol=n)
    endstep<-length(MCMCchain)
    startstep<-max(as.integer(burnin*endstep),1)
    if (pdag) {
      cpdags<-lapply(MCMCchain[startstep:endstep],dagadj2cpadj)
      incidence[which(as.matrix(Reduce('+', cpdags)/(endstep-startstep+1))>p)]<-1
    } else {
      incidence[which(as.matrix(Reduce('+', MCMCchain[startstep:endstep])/(endstep-startstep+1))>p)]<-1
    }
    colnames(incidence)<-varlabels
    rownames(incidence)<-varlabels
    if(DBN) {
      incidence<-DBNcut(incidence,dyn=dyn,b=b)
      incidence.init<-DBNinit(incidence,dyn=dyn,b=b)
      incidence[1:(dyn+b),1:(dyn+b)]<-incidence.init
    }
    return(incidence)
  }
  
  myScore<-BiDAG::scoreparameters(score, data)
  
  partfit<-BiDAG::partitionMCMC(myScore, stepsave = 1)
  
  result <- modelpcore(partfit$traceadd$incidence, 0.5)
  
  return(as.matrix(result))
}

# Running Structure MCMC experiment
runStructMCMC <- function(data,iterations,blklist=NULL,scoretype,sample_parameters){
  
  n <- ncol(data)
  
  # Choose maximum number of parents
  maxparents<-n-1 # Maximum number of parents to allow
  
  # Fill up a matrix with possible parents
  parenttable<-listpossibleparents(maxparents,c(1:n))
  
  moveprobs<-1 # having length 1 disallows the new edge reversal move
  stepsave<-1 #stepsave<-iterations/1000 #how often to save the result
  
  # Initialize the starting DAG
  if (grepl("cg", scoretype)){
    startDAG <- matrix(0, nrow = n, ncol = n) # starting DAG is empty
  } else {
    startDAG<-rDAG(n, 0.5)
  }
  
  if (!is.null(blklist) && all(!is.na(startDAG - blklist))) {
    if (sum(abs(startDAG - blklist), na.rm = TRUE) == 0) {
      startDAG <- matrix(0, nrow = n, ncol = n) # starting DAG is empty
    } else {
      na_mask <- !is.na(startDAG) & !is.na(blklist)
      startDAG[na_mask & (startDAG == blklist) & (startDAG != 0)] <- 0
    }
  } else {
    startDAG <- matrix(0, nrow = n, ncol = n) # starting DAG is empty
  }
  
  colnames(startDAG) <- rownames(startDAG) <- colnames(data)
  
  revallowed<-1 # allow standard edge reversals
  
  example<-structureMCMC(n,data,startDAG,iterations,stepsave,maxparents,parenttable,scoretable=NULL,revallowed,moveprobs,blklist,scoretype,sample_parameters)
  
  burnIn <- floor(0.05 * iterations)
  
  # Extract and process results after burn-in
  incidence_list <- example$incidence[-c(1:burnIn)]
  # score_list <- example$DAGlogscore[-c(1:burnIn)]
  # uniqstructure <- unique(incidence_list)
  # 
  # freq <- tabulate(match(incidence_list, uniqstructure))
  # 
  # bgefreq <- data.frame(
  #   Indx = match(uniqstructure, incidence_list),
  #   StructureStr = I(lapply(uniqstructure, function(x) as.vector(x))),
  #   DAGscores = sapply(match(uniqstructure, incidence_list), function(x) score_list[[x]]),
  #   Frequency = freq
  # )
  
  DAG <- modelpcore(incidence_list, 0.5)
  colnames(DAG) <- rownames(DAG) <- colnames(data)
  
  return(as.matrix(DAG))
}

# Summary of results given a list of experimental results
get_result <- function(results,truemat){
  com <- lapply(1:length(results), function(i) compareDAGs(results[[i]],truemat,cpdag = TRUE))
  results <- sapply(1:8, function (x) mean(as.numeric(lapply(1:length(results), function(i) com[[i]][x])),na.rm = TRUE)) # Compare to True DAGs
  names(results) <- names(com[[1]])
  return(results)
}

get_metric <- function(results,truemat,metric){
  com <- lapply(1:length(results), function(i) compareDAGs(results[[i]],truemat,cpdag = TRUE))
  results <- as.numeric(lapply(1:length(results), function(i) com[[i]][metric])) # Compare to True DAGs
  return(results)
}

MCMCres <- function(MCMC,truemat){
  compDags <- lapply(1:nsim, function(i) compareDAGs(MCMC[[i]][[2]],truemat,cpdag = TRUE))
  results <- sapply(1:8, function (x) mean(as.numeric(lapply(1:nsim, function(i) compDags[[i]][x])),na.rm = TRUE)) # Compare to True DAGs
  names(results) <- names(compDags[[1]])
  return(results)
}

get_metric_MC <- function(MCMC,truemat,metric){
  com <- lapply(1:nsim, function(i) compareDAGs(MCMC[[i]][[2]],truemat,cpdag = TRUE))
  results <- as.numeric(lapply(1:nsim, function(i) com[[i]][metric])) # Compare to True DAGs
  return(results)
}

# Postprocessing functions

### Tabulate results

metric_Scc <- function(metric, datmat, tabur, hcr, partMCMC, tabur_like, hcr_like, structMCMC, tabur_interval_2b, hcr_interval_2b, partMCMC_interval_2b, tabur_interval_4b, hcr_interval_4b, partMCMC_interval_4b, tabur_clust_2b, hcr_clust_2b, partMCMC_clust_2b, tabur_clust_4b, hcr_clust_4b, partMCMC_clust_4b){
  RAG <- cbind(TABU = get_metric(tabur,datmat,metric), HC = get_metric(hcr,datmat,metric), MCMC = get_metric(partMCMC,datmat,metric))
  res_long <- tidyr::gather(as.data.frame(RAG), Method, metric, TABU:MCMC, factor_key=TRUE)
  res_long$Strategy <- "RAG (BGe)"
  
  RAG_like <- cbind(TABU = get_metric(tabur_like,datmat,metric), HC = get_metric(hcr_like,datmat,metric), MCMC = get_metric(structMCMC,datmat,metric))
  RAG_like_long <- tidyr::gather(as.data.frame(RAG_like), Method, metric, TABU:MCMC, factor_key=TRUE)
  RAG_like_long$Strategy <- "RAG (L)"
  
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
  res_long <- rbind(res_long, RAG_like_long, DCRT_inv_long_2b, DCRT_inv_long_4b, DCRT_clust_long_2b, DCRT_clust_long_4b)
  
  return(res_long)
}

metric_comb <- function(metric, datmat, tabur, hcr, partMCMC, tabur_like, hcr_like, structMCMC, tabur_cg, hcr_cg, MCMC_cg, tabur_interval_2b, hcr_interval_2b, partMCMC_interval_2b, tabur_interval_4b, hcr_interval_4b, partMCMC_interval_4b, tabur_clust_2b, hcr_clust_2b, partMCMC_clust_2b, tabur_clust_4b, hcr_clust_4b, partMCMC_clust_4b){
  
  RAG <- cbind(TABU = get_metric(tabur,datmat,metric), HC = get_metric(hcr,datmat,metric), MCMC = get_metric(partMCMC,datmat,metric))
  res_long <- tidyr::gather(as.data.frame(RAG), Method, metric, TABU:MCMC, factor_key=TRUE)
  res_long$Strategy <- "RAG (BGe)"
  
  RAG_like <- cbind(TABU = get_metric(tabur_like,datmat,metric), HC = get_metric(hcr_like,datmat,metric), MCMC = get_metric(structMCMC,datmat,metric))
  RAG_like_long <- tidyr::gather(as.data.frame(RAG_like), Method, metric, TABU:MCMC, factor_key=TRUE)
  RAG_like_long$Strategy <- "RAG (L)"
  
  CLG <- cbind(TABU = get_metric(tabur_cg,datmat,metric), HC = get_metric(hcr_cg,datmat,metric), MCMC = get_metric(MCMC_cg, datmat, metric))
  CLG_long <- tidyr::gather(as.data.frame(CLG), Method, metric, TABU:MCMC, factor_key=TRUE)
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
  res_long <- rbind(res_long, RAG_like_long, DCRT_inv_long_2b, DCRT_inv_long_4b, DCRT_clust_long_2b, DCRT_clust_long_4b)
  
  return(res_long)
}

metric_Sdd <- function(metric, datmat, tabur, hcr, partMCMC, tabur_like, hcr_like, structMCMC, tabur_dcrt, hcr_dcrt, partMCMC_dcrt){
  RAG <- cbind(TABU = get_metric(tabur,datmat,metric), HC = get_metric(hcr,datmat,metric), MCMC = get_metric(partMCMC,datmat,metric))
  res_long <- tidyr::gather(as.data.frame(RAG), Method, metric, TABU:MCMC, factor_key=TRUE)
  res_long$Strategy <- "RAG (BGe)"
  
  RAG_like <- cbind(TABU = get_metric(tabur_like,datmat,metric), HC = get_metric(hcr_like,datmat,metric), MCMC = get_metric(structMCMC,datmat,metric))
  RAG_like_long <- tidyr::gather(as.data.frame(RAG_like), Method, metric, TABU:MCMC, factor_key=TRUE)
  RAG_like_long$Strategy <- "RAG (L)"
  
  DISC <- cbind(TABU = get_metric(tabur_dcrt,datmat,metric), HC = get_metric(hcr_dcrt,datmat,metric), MCMC = get_metric(partMCMC_dcrt,datmat,metric))
  DISC_long <- tidyr::gather(as.data.frame(DISC), Method, metric, TABU:MCMC, factor_key=TRUE)
  DISC_long$Strategy <- "DISC (BDe)"
  res_long <- rbind(res_long, RAG_like_long, DISC_long)
  
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
