# LSAC data
load("~/Documents/BNProject/Data/imputed.dat.all.RData")

library(BiDAG)
library(bnlearn)
library(tidyverse)
library(tibble)
library(ggdag)
library(dagitty)
library(bnlearn)
library(reshape2)
library(igraph)
library(gRain)
library(dplyr)
source('blackpartition.R')

nsim <- 100

# Function to run pmcmc on LSAC data
runPmcmc <- function(score, blacklist, data){
  data <- as.data.frame(data)
  
  myScore<-scoreparameters(score, data)
  
  partfit<-partitionMCMC(myScore, blacklist=blklist, stepsave=1, iterations=100000)
  
  return(partfit)
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

# 1. Run pmcmc on the discretised data
dat.dag <- W2.f
source('makeblacklist.R') 
partcat <- runPmcmc("bdecat", blacklist=blklist, W2.f)

# 2. Set structure
bn_struct <- empty.graph(colnames(W2.f))
amat(bn_struct) <- partcat$DAG
graphviz.plot(bn_struct)
options(digits=3)
bn_mod <- bn.fit(bn_struct, data = W2.f, method = "mle")

# 3. Run pMCMC on simulated datas using "bdecat" score
options(digits = 3)
# Simulate random samples from the structure 
set.seed(97)
B2.f_sim <- lapply(1:nsim, function(i) rbn(bn_mod, 15000))

# Check if any levels do not present in data
# sapply(1:nsim, function(i) sum(sapply(1:ncol(B.W2), function(j) sum(table(B2.f_sim[[i]][,j])==0))))

sim_partcat <- list()
for (i in 1:nsim) {
  dat.dag <- B2.f_sim[[i]]
  source('makeblacklist.R') 
  sim_partcat[[i]] <- runPmcmc("bdecat", blacklist=blklist, dat.dag)
}

#save(sim_partcat, file = paste0("Cat-bdecat",".RData"))

# 4. Obtain DAGs and compare to true DAG
compDAGs_simcat <- lapply(1:nsim,
                          function(i) compareDAGs(sim_partcat[[i]]$DAG, partcat$DAG))

sapply(1:8, function (x) mean(as.numeric(lapply(1:nsim, function(i) compDAGs_simcat[[i]][x]))))

# 5. Run pMCMC on simulated datas using "bge" score

B2.f_gsim <- lapply(1:nsim, function(i) B2.f_sim[[i]] %>% mutate_if(sapply(B2.f_sim[[i]], is.factor), as.numeric))

simbge_part <- list()
for (i in 1:nsim) {
  dat.dag <- B2.f_gsim[[i]]
  source('makeblacklist.R') 
  simbge_part[[i]] <- runPmcmc("bge", blacklist=blklist, dat.dag)
}

#save(simbge_part, file = paste0("Cat-bge",".RData"))

# Comparing with true DAGs
compDAGs_simbge <- lapply(1:nsim,
                          function(i) compareDAGs(simbge_part[[i]]$DAG, partcat$DAG))

sapply(1:8, function (x) mean(as.numeric(lapply(1:nsim, function(i) compDAGs_simbge[[i]][x]))))


####
# 1. Run pmcmc on the original data
B.W2 <- as.data.frame(B.W2)
dat.dag <- B.W2
source('Data/makeblacklist.R') 
partbge <- runPmcmc("bge", blacklist=blklist, dat.dag)

# 2. Obtain DAG structure
bge.struct <- empty.graph(colnames(B.W2))
amat(bge.struct) <- partbge$DAG
graphviz.plot(bge.struct)
options(digits=3)

# 3. Simulate from DAG structure
# Matching original categorical variables
bn.mod <- bn.fit(bge.struct, data = B.W2, method = "mle")
B2_sim.bge <- list()

cat <- c("SX","AC","FS","FH","ME1","FE1","RP1","DP1","FV","HF","HSD","SL","GW")

for (i in c(1:nsim)) {
  temp <- B2_sim.bge[[i]] <-rbn(bn.mod, 7000)
  
  for (j in cat) {
    levs <- levels(factor(B.W2[,j]))
    
    B2_sim.bge[[i]][,j] <- cut(temp[,j], 
                               breaks = quantile(temp[,j], 1/length(levs)*c(0:length(levs))),
                               labels=levs,
                               include.lowest = TRUE)
  }
  
  B2_sim.bge[[i]] <- B2_sim.bge[[i]] %>% mutate_if(sapply(B2_sim.bge[[i]], is.factor), as.numeric) 
}


# 4. Run pmcmc on simulated data - bge

sim_partbge <- list()

for (i in c(1:nsim)) {
  dat.dag <- B2_sim.bge[[i]]
  source('Data/makeblacklist.R') 
  sim_partbge[[i]] <- runPmcmc("bge", blacklist=blklist, dat.dag)
}

gsim_bgedags <- lapply(1:nsim, function(i) sim_partbge[[i]]$DAG)

compDAGs_gsimbge <- lapply(1:nsim,
                           function(i) compareDAGs(gsim_bgedags[[i]], partbge$DAG))

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

sim_partbdecat <- list()

for (i in c(1:nsim)) {
  dat.dag <- B2_sim.bdecat[[i]]
  source('makeblacklist.R') 
  sim_partbdecat[[i]] <- runPmcmc("bdecat", blacklist=blklist, dat.dag)
}

save(sim_partbdecat, file = paste0("Cont-bdecat",".RData"))

# Compare to true dag
compDAGs_simbdecat <- lapply(1:nsim,
                             function(i) compareDAGs(sim_partbdecat[[i]]$DAG, partbge$DAG))

sapply(1:8, function (x) mean(as.numeric(lapply(1:nsim, function(i) compDAGs_simbdecat[[i]][x]))))