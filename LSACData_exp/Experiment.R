# LSAC data
load("imputed.dat.all.RData")

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

## Selecting Wave 2 for analysis

B.W2 <- imputed.dat.k.all[[2]]

B.W2 <- B.W2 %>%
  select(SX,BMI,SE,AC,INC,FS,FH,ME1,FE1,BM1,BM2,RP1,DP1,FV,HF,HSD,SL,OD,GW,BWZ)

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

#W2.cont <- W2.f %>% mutate_if(sapply(W2.f, is.factor), as.numeric)    

#W2.cont <- as.data.frame(W2.cont)


## Discretised dataset - Category option

bdecatScore <- scoreparameters("bge", B.W2)
partcat <- partitionMCMC(bdecatScore, blacklist=blklist, iterations=1000000)
save(partcat, file = paste0("Data/","BWave2_bge",".RData"))

bdecatScore <- scoreparameters("bdecat", W2.f)
partcat <- partitionMCMC(bdecatScore, blacklist=blklist, iterations=1000000)
save(partcat, file = paste0("Data/","BWave2_bdecat",".RData"))

# Plot
m.DAG = empty.graph(colnames(B.W2))
amat(m.DAG) = as.matrix(partcat$DAG)
graphviz.plot(m.DAG, layout = "dot")

