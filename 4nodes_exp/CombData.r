# Function to generate Combination data for 4 nodes

# S_dc
CombData <- function(data, nrep, nvar){
  data$A <- sample(1:2, size = nrep, replace = TRUE) #  has no parents
  data$C <- sample(1:4, size = nrep, replace = TRUE) #  has no parents
  data$B <- data$A * 2 + data$C * 3 + rnorm(nrep) 
  data$D <- data$A + rnorm(nrep) + data$C * 4
  return(data)
}

# S_cd
CombData_cd <- function(data, nrep, nvar, muA, muC){
  data$A <- rnorm(nrep,mean=muA, sd=1) #  has no parents
  data$C <- rnorm(nrep,mean=muC, sd=1) #  has no parents
  # Bern probability
  p <- function(a, b, c, muA, muC) {
    exp(b*(a + c - muA - muC))/(1 + exp(b*(a + c - muA - muC)))}
  for (i in 1:nrep){
    data$B[i] <- rbern(1, p(data$A[i], 2, data$C[i], muA, muC))
    data$D[i] <- rbern(1, p(data$A[i], -1.5, data$C[i], muA, muC))
  }
  return(data)
}

