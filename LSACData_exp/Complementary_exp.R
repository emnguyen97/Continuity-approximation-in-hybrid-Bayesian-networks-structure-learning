# LSAC data
load("LSACData_exp/imputed.dat.all.RData")
# Load the necessary functions
source('StructureMCMC/combinations.R')
source('StructureMCMC/scoretables.R')
source('StructureMCMC/structurefns.R')
source('StructureMCMC/samplefns.R')
source('StructureMCMC/param_utils.R')
source('StructureMCMC/structureMCMC.R')
source('LSACData_exp/blackpartition.R')

library(BiDAG)
library(bnlearn)
library(tidyverse)
library(tibble)
library(ggdag)
library(dagitty)
library(bnlearn)
library(reshape2)
library(igraph)
library(dplyr)
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

nsim <- 100

# Function to run pmcmc on LSAC data
runPmcmc <- function(score, blacklist, data){
  data <- as.data.frame(data)
  
  partition_score = BiDAG::scoreparameters(scoretype = score,
                                    data = data, 
                                    bgepar = list(edgepf = floor(2*log(ncol(data)))))
  
  results_partition = BiDAG::partitionMCMC(scorepar = partition_score,
                                    iterations = 100000, 
                                    blacklist=blklist)
  
  # pip_partition = as.matrix(edgep(results_partition))
  
  return(results_partition)
}

runStrcmcmc <- function(blacklist, data){
  
  # Choose maximum number of parents
  n <- ncol(data)
  maxparents<-n-1 # Maximum number of parents to allow
  
  # Fill up a matrix with possible parents
  
  parenttable<-listpossibleparents(maxparents,c(1:n))
  
  iterations<-100000 #number of iterations in the chain
  moveprobs<-1 # having length 1 disallows the new edge reversal move
  #if(!(length(moveprobs)==1)){print('Vector of move probabilities has the wrong length!')}
  
  stepsave<-1 #stepsave<-iterations/1000 #how often to save the result
  
  set.seed(97)
  
  startDAG<-matrix(numeric(n*n),nrow=n) # starting DAG is empty say
  
  revallowed<-1 # allow standard edge reversals
  
  example<-structureMCMC(n,data,startDAG,iterations,stepsave,maxparents,parenttable,scoretable=NULL,revallowed,moveprobs,blacklist = blacklist,scoretype = "bic-cg")
  
  DAGscores<-unlist(example$DAGlogscore)
  maxDAG<-example$incidence[which(DAGscores==max(DAGscores))][[1]]
  colnames(maxDAG) <- rownames(maxDAG) <- colnames(data)
  
  return(list(DAGlogscores=DAGscores, DAG=maxDAG))
}

## Selecting Wave 2 data

B.W2 <- imputed.dat.k.all[[2]]

B.W2 <- B.W2 %>%
  select(SX,BMI,SE,AC,INC,FS,FH,ME1,FE1,BM1,BM2,RP1,DP1,FV,HF,HSD,SL,OD,GW,BWZ)

# Discretised variables
W2.d <- B.W2 %>% mutate(BMI=cut(BMI, breaks = quantile(BMI, c(0, 0.05, 0.85, 0.95, 1)),
                                labels=c("Underweight","HealthyWeight","Overweight","Obese"),
                                include.lowest = TRUE),
                        SE=cut(SE, breaks = quantile(SE, c(0,0.05, 0.25, 0.75, 0.95, 1)),
                               labels=c("Q1","Q2","Q3","Q4","Q5"),
                               include.lowest = TRUE),
                        INC=cut(INC, breaks = quantile(INC, c(0,0.05, 0.25, 0.75, 0.95, 1)),
                                labels=c("Q1","Q2","Q3","Q4","Q5"),
                                include.lowest = TRUE),
                        BM1=cut(BM1, breaks = c(0,18.5, 25.0, 30.0, Inf),
                                labels=c("Underweight","HealthyWeight","Overweight","Obese"),
                                include.lowest = TRUE),
                        BM2=cut(BM2, breaks = c(0,18.5, 25.0, 30.0, Inf),
                                labels=c("Underweight","HealthyWeight","Overweight","Obese"),
                                include.lowest = TRUE),
                        OD=cut(OD, breaks = quantile(OD, c(0,0.05, 0.25, 0.75, 0.95, 1)),
                               labels=c("Q1","Q2","Q3","Q4","Q5"),
                               include.lowest = TRUE),
                        BWZ=cut(BWZ, breaks = quantile(BWZ, c(0,0.05, 0.25, 0.75, 0.95, 1)),
                                labels=c("Q1","Q2","Q3","Q4","Q5"),
                                include.lowest = TRUE))

W2.f <- W2.d %>% mutate_if(sapply(W2.d, is.numeric), as.factor)    

W2.f <- as.data.frame(W2.f)

n <- ncol(B.W2)

# 1. Run pmcmc on the discretised data
dat.dag <- W2.f
source('LSACData_exp/makeblacklist.R')
partcat <- runPmcmc("bdecat", blacklist=blklist, W2.f)

# 2. Set up the Bayesian network structure from the learned DAG
bn_struct <- empty.graph(colnames(W2.f))
amat(bn_struct) <- as.matrix(partcat$DAG)
graphviz.plot(bn_struct)
options(digits=5)
bn_mod <- bn.fit(bn_struct, data = W2.f, method = "bayes")

# 4. Simulate random samples from the structure
set.seed(1124)

# Step 1: Generate a large dataset with 20,000 samples
B2.f_large <- lapply(1:nsim, function(i) rbn(bn_mod, 20000))

# Step 2: Filter down to a subset of 5,000 samples that includes all levels
B2.f_sim <- lapply(1:nsim, function(i) {
  sim_data <- B2.f_large[[i]]
  
  # Create a sample of 5,000 that includes all levels
  repeat {
    sampled_data <- sim_data[sample(1:nrow(sim_data), 5000), ]
    
    # Check if all levels are present in the sampled data
    all_levels_present <- all(sapply(1:ncol(sampled_data), function(j) {
      all(levels(sim_data[, j]) %in% sampled_data[, j])
    }))
    
    if (all_levels_present) break
  }
  
  return(sampled_data)
})
# Check if any levels do not present in data
# sapply(1:nsim, function(i) sum(sapply(1:ncol(B.W2), function(j) sum(table(B2.f_sim[[i]][,j])==0))))

# Run pMCMC on simulated datas using "bdecat" score
source('LSACData_exp/makeblacklist.R')
sim_partcat <- list()
numCores <- detectCores() - 1
cl <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl) 

sim_partcat <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(B2.f_sim[[i]])
  result$result <- runPmcmc("bdecat", blacklist=blklist, data)
  return(result)
}
stopCluster(cl)

save(sim_partcat, file = paste0("Cat-bdecat",".RData"))

# 4. Obtain DAGs and compare to true DAG
compDAGs_simcat <- lapply(1:nsim,
                          function(i) compareDAGs(sim_partcat[[i]]$DAG, partcat$DAG, cpdag = TRUE))

sapply(1:8, function (x) mean(as.numeric(lapply(1:nsim, function(i) compDAGs_simcat[[i]][x]))))

# 5. Run pMCMC on simulated datas using "bge" score

B2.f_gsim <- lapply(1:nsim, function(i) B2.f_sim[[i]] %>% mutate_if(sapply(B2.f_sim[[i]], is.factor), as.numeric))

simbge_part <- list()
numCores <- detectCores() - 1
cl <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl) 

simbge_part <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  dat.dag <- as.data.frame(B2.f_gsim[[i]])
  result$result <- runPmcmc("bge", blacklist=blklist, dat.dag)
  return(result)
}
stopCluster(cl)

save(simbge_part, file = paste0("Cat-bge",".RData"))

# Comparing with true DAGs
compDAGs_simbge <- lapply(1:nsim,
                          function(i) compareDAGs(simbge_part[[i]]$DAG, partcat$DAG, cpdag = TRUE))

sapply(1:8, function (x) mean(as.numeric(lapply(1:nsim, function(i) compDAGs_simbge[[i]][x]))))

# 6. Run pMCMC on simulated datas using CLG model 
B2.f_clgsim <- lapply(1:nsim, function(i) B2.f_sim[[i]] %>% 
                        mutate_if(names(B.W2) %in% c("BMI","INC","SE","BM1","BM2","BWZ", "OD"), as.numeric))

sim_clgcat <- list()
numCores <- detectCores() - 1
cl <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl) 

sim_clgcat <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  dat.dag <- as.data.frame(B2.f_clgsim[[i]])
  result$result <- runStrcmcmc(blacklist=blklist_clg, dat.dag)
  return(result)
}
stopCluster(cl)
save(sim_clgcat, file = paste0("Cat-CLG",".RData"))

# 7. Obtain DAGs and compare to true DAG
compDAGs_simcat <- lapply(1:nsim,
                          function(i) compareDAGs(sim_clgcat[[i]]$DAG, partcat$DAG, cpdag = TRUE))

sapply(1:8, function (x) mean(as.numeric(lapply(1:nsim, function(i) compDAGs_simcat[[i]][x]))))

####
# 1. Run pmcmc on the original data
B.W2 <- as.data.frame(B.W2)
dat.dag <- B.W2
source('LSACData_exp/makeblacklist.R') 
partbge <- runPmcmc("bge", blacklist=blklist, dat.dag)

# 2. Obtain DAG structure
bge.struct <- empty.graph(colnames(B.W2))
amat(bge.struct) <- as.matrix(partbge$DAG)
graphviz.plot(bge.struct)
options(digits=3)

# 3. Simulate from DAG structure
# Matching original categorical variables
bn.mod <- bn.fit(bge.struct, data = B.W2, method = "mle-g")
B2_sim.bge <- vector("list", nsim)  # Pre-allocate list

cat <- c("SX","AC","FS","FH","ME1","FE1","RP1","DP1","FV","HF","HSD","SL","GW")
cat_levels <- lapply(cat, function(j) levels(factor(B.W2[,j])))  # Precompute levels

for (i in 1:nsim) {
  temp <- rbn(bn.mod, 5000)
  
  # Apply cut only to categorical variables
  for (j in seq_along(cat)) {
    temp[, cat[j]] <- cut(temp[, cat[j]], 
                          breaks = quantile(temp[, cat[j]], probs = seq(0, 1, length.out = length(cat_levels[[j]]) + 1)),
                          labels = cat_levels[[j]],
                          include.lowest = TRUE)
  }
  
  # Convert factors to numeric for categorical variables
  B2_sim.bge[[i]] <- temp %>%
    mutate(across(all_of(cat), as.numeric))
}

# 4. Run pmcmc on simulated data - bge

sim_partbge <- list()
numCores <- detectCores() - 1
cl <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl) 

sim_partbge <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  dat.dag <- as.data.frame(B2_sim.bge[[i]])
  result$result <- runPmcmc("bge", blacklist=blklist, dat.dag)
  return(result)
}
stopCluster(cl)

save(sim_partbge, file = paste0("Cont-bge",".RData"))

gsim_bgedags <- lapply(1:nsim, function(i) sim_partbge[[i]]$DAG)

compDAGs_gsimbge <- lapply(1:nsim,
                           function(i) compareDAGs(gsim_bgedags[[i]], partbge$DAG, cpdag = TRUE))

sapply(1:8, function (x) mean(as.numeric(lapply(1:nsim, function(i) compDAGs_gsimbge[[i]][x]))))

plotdiffs(sim_partbge[[1]]$DAG,partbge$DAG)


# 6. bdecat version

catdiscrt <- function(data) {
  discrtdata <- data %>% mutate(BMI=cut(BMI, 
                                        breaks = quantile(BMI, c(0, 0.05, 0.85, 0.95, 1)),
                                        labels=c("Underweight","HealthyWeight","Overweight","Obese"),
                                        include.lowest = TRUE),
                                SE=cut(SE, breaks = quantile(SE, c(0,0.05, 0.25, 0.75, 0.95, 1)),
                                       labels=c("Q1","Q2","Q3","Q4","Q5"),
                                       include.lowest = TRUE),
                                INC=cut(INC, breaks = quantile(INC, c(0,0.05, 0.25, 0.75, 0.95, 1)),
                                        labels=c("Q1","Q2","Q3","Q4","Q5"),
                                        include.lowest = TRUE),
                                BM1=cut(BM1, breaks = c(0,18.5, 25.0, 30.0, Inf),
                                        labels=c("Underweight","HealthyWeight","Overweight","Obese"),
                                        include.lowest = TRUE),
                                BM2=cut(BM2, breaks = c(0,18.5, 25.0, 30.0, Inf),
                                        labels=c("Underweight","HealthyWeight","Overweight","Obese"),
                                        include.lowest = TRUE),
                                OD=cut(OD, breaks = quantile(OD, c(0,0.05, 0.25, 0.75, 0.95, 1)),
                                       labels=c("Q1","Q2","Q3","Q4","Q5"),
                                       include.lowest = TRUE),
                                BWZ=cut(BWZ, breaks = quantile(BWZ, c(0,0.05, 0.25, 0.75, 0.95, 1)),
                                        labels=c("Q1","Q2","Q3","Q4","Q5"),
                                        include.lowest = TRUE))
  discrtdata <- discrtdata %>% mutate_if(sapply(discrtdata, is.numeric), as.factor) 
  return(discrtdata)
}

B2_sim.bdecat <- lapply(1:nsim, function(i) catdiscrt(B2_sim.bge[[i]]))

source('LSACData_exp/makeblacklist.R')
sim_partbdecat <- list()
numCores <- detectCores() - 1
cl <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl) 

sim_partbdecat <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(B2_sim.bdecat[[i]])
  result$result <- runPmcmc("bdecat", blacklist=blklist, data)
  return(result)
}
stopCluster(cl)

save(sim_partbdecat, file = paste0("Cont-bdecat",".RData"))

# Compare to true dag
compDAGs_simbdecat <- lapply(1:nsim,
                             function(i) compareDAGs(sim_partbdecat[[i]]$DAG, partbge$DAG, cpdag = TRUE))

sapply(1:8, function (x) mean(as.numeric(lapply(1:nsim, function(i) compDAGs_simbdecat[[i]][x]))))

# 7. Run structure learning on simulated data - CLG

B2_sim.clg <- lapply(1:nsim, function(i) B2_sim.bdecat[[i]] %>% 
                        mutate_if(names(B.W2) %in% c("BMI","INC","SE","BM1","BM2","BWZ", "OD"), as.numeric))

sim_partclg <- list()
numCores <- detectCores() - 1
cl <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl) 

sim_partclg <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data_clg <- as.data.frame(B2_sim.clg[[i]])
  result$result <- runStrcmcmc(blacklist=blklist_clg, data_clg)
  return(result)
}
stopCluster(cl)
save(sim_partclg, file = paste0("Cont-CLG",".RData"))


sim_clgdags <- lapply(1:nsim, function(i) sim_partclg[[i]]$DAG)

compDAGs_simclg <- lapply(1:nsim,
                          function(i) compareDAGs(sim_clgdags[[i]], partbge$DAG, cpdag = TRUE))

sapply(1:8, function (x) mean(as.numeric(lapply(1:nsim, function(i) compDAGs_simclg[[i]][x]))))

#####
# 1. Obtain DAG from CLG model on the original data
B.W2_clg<- B.W2 %>% mutate_if(!names(B.W2) %in% c("BMI","INC","SE","BM1","BM2","BWZ"), factor)
dat.dag <- as.data.frame(B.W2_clg)
source('LSACData_exp/makeblacklist.R')
clgmod <- runStrcmcmc(blacklist=blklist_clg, data=dat.dag)

# 2. Obtain DAG structure
clg.struct <- empty.graph(colnames(B.W2))
amat(clg.struct) <- as.matrix(clgmod$DAG)
graphviz.plot(clg.struct)
options(digits=3)

# 3. Simulate from DAG structure
bn.mod <- bn.fit(clg.struct, data = dat.dag, method = "mle-cg")

set.seed(1124)
# Step 1: Generate a large dataset with 20,000 samples
B2_clg_large <- lapply(1:nsim, function(i) rbn(bn.mod, 20000))

# Step 2: Filter down to a subset of 5,000 samples that includes all levels
B2_sim.clg <- lapply(1:nsim, function(i) {
  sim_data <- B2_clg_large[[i]]
  
  # Create a sample of 5,000 that includes all levels
  repeat {
    sampled_data <- sim_data[sample(1:nrow(sim_data), 5000), ]
    
    # Check if all levels are present in the sampled data
    all_levels_present <- all(sapply(1:ncol(sampled_data), function(j) {
      all(levels(sim_data[, j]) %in% sampled_data[, j])
    }))
    
    if (all_levels_present) break
  }
  
  return(sampled_data)
})

# 4. Run structure learning on simulated data - CLG

sim_clg <- list()
numCores <- detectCores() - 1
cl <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl) 

sim_clg <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  data <- as.data.frame(B2_sim.clg[[i]])
  result$result <- runStrcmcmc(blacklist=blklist_clg, data)
  return(result)
}
stopCluster(cl)
save(sim_clg, file = paste0("Clg-CLG",".RData"))


sim_clgdags <- lapply(1:nsim, function(i) sim_clg[[i]]$DAG)

compDAGs_simclg <- lapply(1:nsim,
                          function(i) compareDAGs(sim_clgdags[[i]], clgmod$DAG, cpdag = TRUE))

sapply(1:8, function (x) mean(as.numeric(lapply(1:nsim, function(i) compDAGs_simclg[[i]][x]))))


# 5. Run structure learning on simulated data - Bge
B2.bgesim <- lapply(1:nsim, function(i) B2_sim.clg[[i]] %>% 
                      mutate_if(!names(B.W2) %in% c("BMI","INC","SE","BM1","BM2","BWZ"), as.numeric))

sim_bge <- list()
numCores <- detectCores() - 1
cl <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl) 

sim_bge <- foreach(i = seq_len(nsim), .combine = 'c') %dopar% {
  result <- multiResultClass()
  dat.dag <- as.data.frame(B2.bgesim[[i]])
  result$result <- runPmcmc("bge", blacklist=blklist, dat.dag)
  return(result)
}
stopCluster(cl)
save(sim_bge, file = paste0("Clg-bge",".RData"))

sim_bgedags <- lapply(1:nsim, function(i) sim_bge[[i]]$DAG)

compDAGs_simbge <- lapply(1:nsim,
                          function(i) compareDAGs(sim_bgedags[[i]], clgmod$DAG, cpdag = TRUE))

sapply(1:8, function (x) mean(as.numeric(lapply(1:nsim, function(i) compDAGs_simbge[[i]][x]))))

# 6. Run structure learning on simulated data - bdecat
discrt <- function(data) {
  discrtdata <- data %>% mutate(BMI=cut(BMI, 
                                        breaks = quantile(BMI, c(0, 0.05, 0.85, 0.95, 1)),
                                        labels=c("Underweight","HealthyWeight","Overweight","Obese"),
                                        include.lowest = TRUE),
                                SE=cut(SE, breaks = quantile(SE, c(0,0.05, 0.25, 0.75, 0.95, 1)),
                                       labels=c("Q1","Q2","Q3","Q4","Q5"),
                                       include.lowest = TRUE),
                                INC=cut(INC, breaks = quantile(INC, c(0,0.05, 0.25, 0.75, 0.95, 1)),
                                        labels=c("Q1","Q2","Q3","Q4","Q5"),
                                        include.lowest = TRUE),
                                BM1=cut(BM1, breaks = c(0,18.5, 25.0, 30.0, Inf),
                                        labels=c("Underweight","HealthyWeight","Overweight","Obese"),
                                        include.lowest = TRUE),
                                BM2=cut(BM2, breaks = c(0,18.5, 25.0, 30.0, Inf),
                                        labels=c("Underweight","HealthyWeight","Overweight","Obese"),
                                        include.lowest = TRUE),
                                BWZ=cut(BWZ, breaks = quantile(BWZ, c(0,0.05, 0.25, 0.75, 0.95, 1)),
                                        labels=c("Q1","Q2","Q3","Q4","Q5"),
                                        include.lowest = TRUE))
  #discrtdata <- discrtdata %>% mutate_if(sapply(discrtdata, is.numeric), as.factor) 
  return(discrtdata)
}

B2.bdecatsim <- lapply(1:nsim, function(i) discrt(B2_sim.clg[[i]]))

sim_bdecat <- list()
numCores <- detectCores() - 1
cl <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl) 

sim_bdecat <- foreach(i = seq_len(100), .combine = 'c') %dopar% {
  result <- multiResultClass()
  dat.dag <- as.data.frame(B2.bdecatsim[[i]])
  result$result <- runPmcmc("bdecat", blacklist=blklist, dat.dag)
  return(result)
}
stopCluster(cl)
save(sim_bdecat, file = paste0("Clg-bdecat",".RData"))

sim_bdecatdags <- lapply(1:nsim, function(i) sim_bdecat[[i]]$DAG)

compDAGs_simbdecat <- lapply(1:nsim, function(i) compareDAGs(sim_bdecatdags[[i]], clgmod$DAG, cpdag = TRUE))

sapply(1:8, function (x) mean(as.numeric(lapply(1:nsim, function(i) compDAGs_simbdecat[[i]][x]))))

