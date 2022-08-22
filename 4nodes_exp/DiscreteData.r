# Function to generate Discrete data fir 4 nodes

rbern <- function(n, p){
  sims <- sample(1:2, size = n, replace = TRUE, prob = c(1-p, p))
  return(sims)
}

DiscreteData <- function(data, nrep){
  data$A <- rbern(nrep, 0.5)
  data$C <- sample(1:4, size = nrep, replace = TRUE) 
  for (i in 1:nrep){
    if(data$A[i] == 1 && data$C[i] == 1) {
      data$B[i] <- rbern(1, 0.05)
      data$D[i] <- rbern(1, 0.95)
    } else if(data$A[i] == 1 && data$C[i] == 2) {
      data$B[i] <- rbern(1, 0.1)
      data$D[i] <- rbern(1, 0.9)
    } else if(data$A[i] == 1 && data$C[i] == 3) {
      data$B[i] <- rbern(1, 0.3)
      data$D[i] <- rbern(1, 0.7)
    } else if(data$A[i] == 1 && data$C[i] == 4) {
      data$B[i] <- rbern(1, 0.7)
      data$D[i] <- rbern(1, 0.3)
    } else if(data$A[i] == 2 && data$C[i] == 1) {
      data$B[i] <- rbern(1, 0.1)
      data$D[i] <- rbern(1, 0.9)
    } else if(data$A[i] == 2 && data$C[i] == 2) {
      data$B[i] <- rbern(1, 0.3)
      data$D[i] <- rbern(1, 0.7)
    } else if(data$A[i] == 2 && data$C[i] == 3) {
      data$B[i] <- rbern(1, 0.7)
      data$D[i] <- rbern(1, 0.3)
    } else {
      data$B[i] <- rbern(1, 0.95)
      data$D[i] <- rbern(1, 0.05)
    }
  }
  data <- data %>% mutate_if(sapply(data, is.numeric), as.factor)
  return(data)
}
