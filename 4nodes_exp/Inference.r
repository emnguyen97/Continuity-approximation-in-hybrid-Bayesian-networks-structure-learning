# Function to run partition MCMC from BiDAG and return results as a summary table
# Summary table of unique dags over 100000 iterations, with their respective index, adjacency matrix in strings, Dag score, and number of times they appear.

Inference <- function(score, data){
  data <- as.data.frame(data)
  
  myScore<-scoreparameters(score, data)
  
  partfit<-partitionMCMC(myScore, stepsave=1, iterations=100000)
  
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
  return(list(bgefreq, partfit$DAG))
}
