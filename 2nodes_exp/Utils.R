## Helper functions
library(dplyr)

# Function to generate random DAG with probability p
rDAG <- function(n, p) {
  DAG = matrix(0, n, n)
  colnames(DAG) = rownames(DAG) = 1:n
  DAG[lower.tri(DAG)] = stats::rbinom(n = n * (n - 1)/2, size = 1, 
                                      prob = p)
  return(DAG)
}

# Function to generate Continuous data
Scc_Data <- function(data, beta, nrep, nvar){
  data$A <- rnorm(nrep,mean=-1, sd=1) #  has no parents
  data$B <- data$A * beta + rnorm(nrep)
  return(data)
}

# Function to discretise continuous data
Discretised <- function(x) {
  cut(x, breaks = c(min(x), (min(x) + max(x)) / 2, max(x)),
      labels = c("a", "b"), include.lowest = TRUE)
}

Scc_Discretised <- function(data){
  data <- data %>%
    mutate(
      A = Discretised(A),
      B = Discretised(B)
    )
  return(data)
}

# Function to generate combination data

# Scd_Data <- function(data, muA, b, nrep, char=FALSE){
#   data$A <- rnorm(nrep,mean=muA, sd=1) #  has no parents
#   # Bern probability
#   p <- function(a, muA) exp(b*(a - muA))/(1 + exp(b*(a - muA)))
#   for (i in 1:nrep){
#     data$B[i] <- rbern(1, p(data$A[i], muA), char)
#   }
#   return(data)
# }

rbern <- function(n, p, char=FALSE){
  if (char) {
    # Directly sample characters
    sims <- sample(c("a", "b"), size = n, replace = TRUE, prob = c(1-p, p))
  } else {
    # Directly sample 0s and 1s
    sims <- rbinom(n, size = 1, prob = p)
  }
  return(sims)
}

Scd_Data <- function(data, muA, b, nrep, char=FALSE){
  # Generate all values of A in one step
  data$A <- rnorm(nrep, mean = muA, sd = 1)  # A has no parents
  
  # Calculate probabilities for B based on A in one step
  probs <- exp(b * (data$A - muA)) / (1 + exp(b * (data$A - muA)))
  
  # Generate all values of B in one step using vectorized rbern
  data$B <- rbern(nrep, probs, char)
  
  return(data)
}


# Sdc_Data <- function(data, beta, nrep, char=FALSE){
#   
#   if (char == TRUE){
#     data$A <-as.factor(data$A)
#     levels(data$A) <- c("a","b")
#   } else {
#     data$A <- rbern(nrep, p=0.5, char) #  has no parents
#   }
#   data$B <- rnorm(nrep, mean = data$A*beta + 1, sd = 1)
#   
#   return(data)
# }

Sdc_Data <- function(data, beta, nrep, char=FALSE){
  
  if (char) {
    data$A <- factor(rbern(nrep, 0.5, char), levels = c("a", "b"))
  } else {
    data$A <- rbern(nrep, p = 0.5, char) #  has no parents
  }
  
  # Vectorized calculation for B
  data$B <- rnorm(nrep, mean = data$A*beta + 1, sd = 1)
  
  return(data)
}


# Discretised
Scd_Discretised <- function(data, clg = FALSE){
  if (!clg) {
    data <- data %>%
      mutate(A = Discretised(A))
  }
  
  data <- data %>%
    mutate(B = as.factor(B))
  
  return(data)
}

Sdc_Discretised <- function(data, clg = FALSE){
  if (clg) {
    data <- data %>%
      mutate(
        A = as.factor(A)
      )
  } else {
    data <- data %>%
      mutate(
        B = Discretised(B),
        A = as.factor(A)
      )
  }
  return(data)
}

# Sdd_Data <- function(data, beta, nrep){
#   data$A <- rbern(nrep, 0.2)
#   for (i in 1:nrep){
#     if(data$A[i] == 0) {
#       data$B[i] <- rbern(1, 1-beta/2-0.5)
#     } else if(data$A[i] == 1) {
#       data$B[i] <- rbern(1, beta/2 + 0.5)
#     } 
#   }
#   data <- data %>% mutate_if(sapply(data, is.numeric), as.factor)
#   return(data)
# }

Sdd_Data <- function(data, beta, nrep){
  data$A <- rbern(nrep, 0.2)
  
  # Vectorized calculation for B
  probs <- ifelse(data$A == 0, 1 - beta/2 - 0.5, beta/2 + 0.5)
  data$B <- rbern(nrep, probs)
  
  # Convert numeric columns to factors in one step
  data <- data %>% mutate(across(where(is.numeric), as.factor))
  
  return(data)
}


# Function to run partition MCMC from BiDAG and return results as a summary table
# Summary table of unique dags over iterations, with their respective index, adjacency matrix in strings, Dag score, and number of times they appear.

runPartMCMC <- function(score, data, iterations){
  data <- as.data.frame(data)
  
  myScore<-BiDAG::scoreparameters(score, data)
  
  partfit<-BiDAG::partitionMCMC(myScore, stepsave=1, iterations=iterations)
  
  DAG <- modelpcore(partfit$traceadd$incidence, 0.5)
  
  options(digits = 22)
  
  uniqstructure <- unique(partfit$traceadd$incidence)
  
  freq <- tabulate(match(partfit$traceadd$incidence, unique(partfit$traceadd$incidence)))
  
  bgefreq <- matrix(NA,nrow=length(uniqstructure),ncol = 4)
  colnames(bgefreq) <- c("Indx", "StructureStr", "DAGscores", "Frequency")
  bgefreq <- as.data.frame(bgefreq)
  
  bgefreq$Frequency <- freq
  bgefreq$Indx <- match(unique(partfit$traceadd$incidence),partfit$traceadd$incidence)
  bgefreq$DAGscores <- lapply(bgefreq$Indx, function(x) partfit$trace[[x]])
  bgefreq$StructureStr <- lapply(bgefreq$Indx, function(x) as.vector(partfit$traceadd$incidence[[x]]))
  return(list(summary=bgefreq, DAGlogscores=partfit$trace, DAG=as.matrix(DAG)))
}

modelpcore<-function(MCMCchain, p, pdag=FALSE, burnin=0.0, DBN=FALSE, nsmall=0, dyn=0, b=0) {
  
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

# Running Structure MCMC experiment
runStructMCMC <- function(data,iterations,blklist,scoretype,sample_parameters){
  
  n <- ncol(data)
  
  # Choose maximum number of parents
  
  maxparents<-n-1 # Maximum number of parents to allow
  
  # Fill up a matrix with possible parents
  
  parenttable<-listpossibleparents(maxparents,c(1:n))
  tablelength<-nrow(parenttable[[1]]) # size of the table
  moveprobs<-1 # having length 1 disallows the new edge reversal move
  #if(!(length(moveprobs)==1)){print('Vector of move probabilities has the wrong length!')}
  burnIn <- 0.05*iterations
  
  stepsave<-1 #stepsave<-iterations/1000 #how often to save the result
  
  startDAG<-rDAG(n, 0.5)
  if (sum(abs(startDAG-blklist)) == 0){
    startDAG<-matrix(numeric(n*n),nrow=n) # starting DAG is empty say
  }
  #startDAG<-matrix(numeric(n*n),nrow=n) # starting DAG is empty say
  #set.seed(2024)
  
  revallowed<-1 # allow standard edge reversals
  
  example<-structureMCMC(n,data,startDAG,iterations,stepsave,maxparents,parenttable,scoretable=NULL,revallowed,moveprobs,blklist,scoretype,sample_parameters)
  
  incidence_list <- example$incidence[-c(1:burnIn)]
  score_list <- example$DAGlogscore[-c(1:burnIn)]
  #incidence_list <- example$incidence
  #score_list <- example$DAGlogscore
  uniqstructure <- unique(incidence_list)
  
  freq <- tabulate(match(incidence_list, uniqstructure))
  
  bgefreq <- matrix(NA,nrow=length(uniqstructure),ncol = 4)
  colnames(bgefreq) <- c("Indx", "StructureStr", "DAGscores", "Frequency")
  bgefreq <- as.data.frame(bgefreq)
  bgefreq$Frequency <- freq
  bgefreq$Indx <- match(uniqstructure,incidence_list)
  bgefreq$DAGscores <- sapply(bgefreq$Indx, function(x) score_list[[x]])
  bgefreq$StructureStr <- lapply(bgefreq$Indx, function(x) as.vector(incidence_list[[x]]))
  
  DAGscores<-unlist(score_list)
  #maxDAG<-incidence_list[which(DAGscores==max(DAGscores))][[1]]
  DAG <- modelpcore(incidence_list, 0.5)
  colnames(DAG) <- rownames(DAG) <- colnames(data)
  
  return(list(summary=bgefreq, DAGlogscores=DAGscores, DAG=DAG))
}

multiResultClass <- function(result=NULL){
  me <- list(
    result = result
  )
  
  ## Set the name for the class
  class(me) <- append(class(me),"multiResultClass")
  return(me)
}

# Metric: Getting frequency ratios

getIdx <- function(n, freqtable, k){
  which(sapply(1:n, function(i) identical(freqtable[[k]]$StructureStr[[i]], rep(0,4))))
}

# Get the frequency for the correct dag
get_freq <- function(freqtable, k, idx){
  freqtable[[k]]$Frequency[idx]
}

F_ratio <- function(freqtable){
  
  iter <- sum(freqtable[[1]]$Frequency)-1
  
  # Get the index of incorrect dags for each of the datasets
  inc_indx <- sapply(1:length(freqtable), function(x) getIdx(nrow(freqtable[[x]]),freqtable, x))
  
  # Get the frequency of incorrect dags for each of the datasets
  inc_freq <- sapply(1:length(freqtable), function(x) get_freq(freqtable, x, inc_indx[[x]])) 
  
  inc_freq_unlist <- as.numeric(sapply(inc_freq, function(s) if (length(s) == 0) 0 else paste(s, collapse = " ")))
  
  sum(iter - inc_freq_unlist)/sum(inc_freq_unlist)
}

# Function to extract results
post_process <- function(res_object){
  # Extracting frequency table for unique structures
  Freqtab <- lapply(1:length(res_object), function(i) res_object[[i]][[1]])
  
  compDags <- lapply(1:length(res_object), function(i) BiDAG::compareDAGs(res_object[[i]][[3]],truemat, cpdag = TRUE))
  
  options(digits = 5)
  res <- sapply(1:8, function (x) mean(as.numeric(lapply(1:length(res_object), function(i) compDags[[i]][x])),na.rm = TRUE))
  res <- c(res,F_ratio(Freqtab))
  
  names(res) <- c(names(compDags[[1]]),"FR")
  
  return(res) # Compare to True DAGs
}
