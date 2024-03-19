## Generate a Directed Acyclic Graph (DAG) randomly with n nodes
rDAG <- function(n, p){
  adjmat <- matrix(0, n, n)
  colnames(adjmat) <- rownames(adjmat) <- LETTERS[1:n]
  adjmat[lower.tri(adjmat)] <- stats::rbinom(n = n*(n-1)/2, size = 1, prob = p)
  return(adjmat)
}

## Plot adjacency matrix
plotadj <- function(adjmat){
  n <- dim(adjmat)[1]
  e <- empty.graph(colnames(adjmat))
  amat(e) <- adjmat
  return(graphviz.plot(e))
}

##### S_cc Data Generating Functions
Scc_Data <- function(datmat, W_mat, D, nrep, nvar){
  # Nodes with/without parents
  wp <- which(colSums(datmat) != 0)
  np <- which(colSums(datmat) == 0)
  
  # Creating empty matrix
  data <- matrix(NA,nrow=nrep,ncol = nvar)
  colnames(data) <- colnames(datmat)
  data <- as.data.frame(data)
  
  # Outpoints are simulated from multivariate normal distribution
  for (i in np){
    data[,i] <- rnorm(nrep, sd=sqrt(D[i,i]))
  }
  
  # Second layer nodes are generated using true coefficients with parent nodes
  layer1 <- which(colSums(datmat[-np,-np]) == 0)
  
  for (j in names(layer1)){
    pa <- which(W_mat[,j] != 0)
    if (length(pa) > 1) {
      data[,j] <- as.vector(as.matrix(data[,pa]) %*% matrix(W_mat[,j][W_mat[,j]!=0])) + rnorm(nrep)
    } else {
      data[,j] <- data[,pa] * W_mat[,j][W_mat[,j]!=0] + rnorm(nrep)
    }
  }
  
  # Final layer nodes 
  nincld <- c(names(np), names(layer1))
  layer2 <- setdiff(colnames(datmat),nincld)
  
  for (j in layer2){
    pa <- which(W_mat[,j] != 0)
    if (length(pa) > 1) {
      data[,j] <- as.vector(as.matrix(data[,pa]) %*% matrix(W_mat[,j][W_mat[,j]!=0])) + rnorm(nrep)
    } else {
      data[,j] <- data[,pa] * W_mat[,j][W_mat[,j]!=0] + rnorm(nrep)
    }
  }
  
  return(data)
}

##### S_cd Data Generating Functions
# Bern probability
prb <- function(data, coeff, mu_pa) {
  exp(coeff*(sum(data) - sum(mu_pa)))/(1 + exp(coeff*(sum(data) - sum(mu_pa))))}

rbern <- function(n, p){
  sims <- sample(c("1","2"), size = n, replace = TRUE, prob = c(1-p, p))
  return(sims)
}

Scd_Data <- function(datmat, W_mat, D, nrep, nvar){
  # Nodes with/without parents
  wp <- which(colSums(datmat) != 0)
  np <- which(colSums(datmat) == 0)
  
  # Creating empty matrix
  data <- matrix(NA,nrow=nrep,ncol = nvar)
  colnames(data) <- colnames(datmat)
  data <- as.data.frame(data)
  
  # Outpoints are simulated from multivariate normal distribution
  #data[,np] <- mvtnorm::rmvnorm(n=nrep, sigma=Sigma[np,np])
  for (i in np){
    data[,i] <- rnorm(nrep, sd=sqrt(D[i,i]))
  }
  
  # Second layer nodes are generated using true coefficients with parent nodes
  layer1 <- which(colSums(datmat[-np,-np]) == 0)
  
  for (j in names(layer1)){
    pa <- which(W_mat[,j] != 0)
    if (length(pa) > 1) {
      data[,j] <- as.vector(as.matrix(data[,pa]) %*% matrix(W_mat[,j][W_mat[,j]!=0])) + rnorm(nrep)
    } else {
      data[,j] <- data[,pa] * W_mat[,j][W_mat[,j]!=0] + rnorm(nrep)
    }
  }
  
  # Final layer nodes are discrete RVs
  nincld <- c(names(np), names(layer1))
  layer2 <- setdiff(colnames(datmat),nincld)
  
  for (j in layer2){
    pa <- names(which(W_mat[,j] != 0))
    if (length(pa) > 1) {
      data[,j] <- as.factor(sapply(1:nrep, function(i) rbern(1, prb(data[,pa][i,], coeff=2, mu_pa=rep(0,length(pa))))))
    } else {
      data[,j] <- as.factor(sapply(1:nrep, function(i) rbern(1, prb(data[,pa][i], coeff=-2, mu_pa=rep(0,length(pa))))))
    }
  }
  
  return(data)
}

##### S_dc Data Generating Functions
rbern_num <- function(n, p){
  sims <- sample(c(1:2), size = n, replace = TRUE, prob = c(1-p, p))
  return(sims)
}

Dcrt_1pa <- function(data, pa1, child, p, flag=0, nrep){
  
  if (flag == 1){
    for (i in 1:nrep){
      if(data[,pa1][i] == 1) {
        data[,child][i] <- rbern_num(1, p)
      } else if(data[,pa1][i] == 2) {
        data[,child][i] <- rbern_num(1, min(1, p + 0.2))
      } else if(data[,pa1][i] == 3) {
        data[,child][i] <- rbern_num(1, min(1, p + 0.3))
      } else {
        data[,child][i] <- rbern_num(1, min(1, p + 0.4))
      }
    }
  } else {
    for (i in 1:nrep){
      if(data[,pa1][i] == 1) {
        data[,child][i] <- rbern_num(1, p)
      } else if(data[,pa1][i] == 2) {
        data[,child][i] <- rbern_num(1, max(0, p - 0.1))
      } else if(data[,pa1][i] == 3) {
        data[,child][i] <- rbern_num(1, max(0, p - 0.3))
      } else {
        data[,child][i] <- rbern_num(1, max(0, p - 0.5))
      }
    }
  }
  
  return(data[,child])
}

Dcrt_1pa_v2 <- function(data, pa1, child1, child2, nrep){
  
  for (i in 1:nrep){
    if(data[,pa1][i] == 1) {
      data[,child1][i] <- rbern_num(1, 0.3)
      data[,child2][i] <- rbern_num(1, 0.9)
      #data[,child3][i] <- rbern_num(1, 0.2)
    } else if(data[,pa1][i] == 2) {
      data[,child1][i] <- rbern_num(1, 0.5)
      data[,child2][i] <- rbern_num(1, 0.7)
      #data[,child3][i] <- rbern_num(1, 0.4)
    } else if(data[,pa1][i] == 3) {
      data[,child1][i] <- rbern_num(1, 0.7)
      data[,child2][i] <- rbern_num(1, 0.5)
      #data[,child3][i] <- rbern_num(1, 0.6)
    } else {
      data[,child1][i] <- rbern_num(1, 0.9)
      data[,child2][i] <- rbern_num(1, 0.3)
      #data[,child3][i] <- rbern_num(1, 0.8)
    }
  }
  
  return(data[,c(child1,child2)])
}

Dcrt_2pas <- function(data, pa1, pa2, flag=0, child, nrep){
  # Pa1 has 2 categories, Pa2 has 4 categories
  if (flag == 1){
    for (i in 1:nrep){
      if(data[,pa1][i] == 1 && data[,pa2][i] == 1) {
        data[,child][i] <- rbern_num(1, 0.2)
      } else if(data[,pa1][i] == 1 && data[,pa2][i] == 2) {
        data[,child][i] <- rbern_num(1, 0.4)
      } else if(data[,pa1][i] == 1 && data[,pa2][i] == 3) {
        data[,child][i] <- rbern_num(1, 0.6)
      } else if(data[,pa1][i] == 1 && data[,pa2][i] == 4) {
        data[,child][i] <- rbern_num(1, 0.8)
      } else if(data[,pa1][i] == 2 && data[,pa2][i] == 1) {
        data[,child][i] <- rbern_num(1, 0.9)
      } else if(data[,pa1][i] == 2 && data[,pa2][i] == 2) {
        data[,child][i] <- rbern_num(1, 0.7)
      } else if(data[,pa1][i] == 2 && data[,pa2][i] == 3) {
        data[,child][i] <- rbern_num(1, 0.5)
      } else {
        data[,child][i] <- rbern_num(1, 0.3)
      }
    }
  } else {
    for (i in 1:nrep){
      if(data[,pa1][i] == 1 && data[,pa2][i] == 1) {
        data[,child][i] <- rbern_num(1, 0.9)
      } else if(data[,pa1][i] == 1 && data[,pa2][i] == 2) {
        data[,child][i] <- rbern_num(1, 0.7)
      } else if(data[,pa1][i] == 1 && data[,pa2][i] == 3) {
        data[,child][i] <- rbern_num(1, 0.5)
      } else if(data[,pa1][i] == 1 && data[,pa2][i] == 4) {
        data[,child][i] <- rbern_num(1, 0.3)
      } else if(data[,pa1][i] == 2 && data[,pa2][i] == 1) {
        data[,child][i] <- rbern_num(1, 0.2)
      } else if(data[,pa1][i] == 2 && data[,pa2][i] == 2) {
        data[,child][i] <- rbern_num(1, 0.4)
      } else if(data[,pa1][i] == 2 && data[,pa2][i] == 3) {
        data[,child][i] <- rbern_num(1, 0.6)
      } else {
        data[,child][i] <- rbern_num(1, 0.8)
      }
    }
  }
  
  return(data[,child])
}

Sdc_Data <- function(datmat, W_mat, D, nrep, nvar){
  # Nodes with/without parents
  wp <- which(colSums(datmat) != 0)
  np <- which(colSums(datmat) == 0)
  # First layer with parents
  layer1 <- which(colSums(datmat[-np,-np]) == 0)
  nincld <- c(names(np), names(layer1))
  # Second layer with parents
  layer2 <- setdiff(colnames(datmat),nincld)
  
  # Creating empty matrix
  data <- matrix(0,nrow=nrep,ncol = nvar)
  colnames(data) <- colnames(datmat)
  data <- as.data.frame(data)
  
  # Simulate outpoints to be categorical variables
  data[,"J"] <- sample(c(1:2), size = nrep, replace = TRUE)
  data[,"E"] <- sample(c(1:4), size = nrep, replace = TRUE)
  
  # Simulate nodes from second layer other than D as categorical variables
  re <- "D"
  layer1 <- layer1[names(layer1) != re]
  for (j in names(layer1)){
    pa <- names(which(W_mat[,j] != 0))
    if (length(pa) == 1){
      prop <- sample(c(0.2, 0.9),1)
      if (prop == 0.2) {flag = 1} else {flag = 0}
      data[,j] <- Dcrt_1pa(data, pa1=pa, child = j, p = prop, flag=flag, nrep=nrep)
    }else {
      if (j == "F"){flag=0} else {flag=1}
      pa1 <- names(which(sapply(apply(data[,pa], 2, unique), function(x) length(x)) == 2))
      pa2 <- setdiff(pa, pa1)
      data[,j] <- Dcrt_2pas(data, pa1, pa2, flag=flag, child=j, nrep)
    }
  }
  #same_pa <- names(which(sapply(names(layer1), function(x) names(which(W_mat[,x] != 0)) == "J")))
  #data[,same_pa] <- Dcrt_1pa_v2(data, pa1="J", child1 = same_pa[1], child2 = same_pa[2], child3 = same_pa[3], nrep)
  #same_pa <- intersect(names(which(W_mat["J",] != 0)),names(layer1))
  #data[,same_pa] <- Dcrt_1pa_v2(data, pa1="J", child1 = same_pa[1], child2 = same_pa[2], nrep)
  
  #re <- setdiff(names(layer1), same_pa)
  #data[,re] <- Dcrt_1pa(data, pa1="E", child = re, p = runif(1,0.1,1), nrep)
  
  # Simulate nodes from third layer as continuous variables
  for (j in c(re, layer2)){
    pa <- which(W_mat[,j] != 0)
    if (length(pa) > 1) {
      data[,j] <- as.vector(as.matrix(data[,pa]) %*% matrix(W_mat[,j][W_mat[,j]!=0])) + rnorm(nrep)
    } else {
      data[,j] <- data[,pa] * W_mat[,j][W_mat[,j]!=0] + rnorm(nrep)
    }
  }
  
  return(data)
}

##### S_dd Data Generating Functions

Sdd_Data <- function(datmat, W_mat, nrep, nvar){
  # Nodes with/without parents
  wp <- which(colSums(datmat) != 0)
  np <- which(colSums(datmat) == 0)
  layer1 <- which(colSums(datmat[-np,-np]) == 0)
  nincld <- c(names(np), names(layer1))
  layer2 <- setdiff(colnames(datmat),nincld)
  
  # Creating empty matrix
  data <- matrix(0,nrow=nrep,ncol = nvar)
  colnames(data) <- colnames(datmat)
  data <- as.data.frame(data)
  
  # Simulate outpoints to be categorical variables
  data[,"J"] <- sample(c(1:4), size = nrep, replace = TRUE)
  data[,"E"] <- sample(c(0:1), size = nrep, replace = TRUE)
  
  # Simulate nodes from second layer as categorical variables
  #same_pa <- names(which(sapply(names(layer1), function(x) names(which(W_mat[,x] != 0)) == "J")))
  #data[,same_pa] <- Dcrt_1pa_v2(data, pa1="J", child1 = same_pa[1], child2 = same_pa[2], child3 = same_pa[3], nrep)
  #for (j in np) {
  #  same_pa <- intersect(names(which(W_mat[j,] != 0)),names(layer1))
  #  data[,same_pa] <- Dcrt_1pa_v2(data, pa1=j, child1 = same_pa[2], child2 = same_pa[1], nrep)
  #}
  for (j in names(layer1)){
    pa <- names(which(W_mat[,j] != 0))
    if (length(pa) == 1){
      prop <- sample(c(0.3, 0.9),1)
      if (prop == 0.3) {flag = 1} else {flag = 0}
      data[,j] <- Dcrt_1pa(data, pa1=pa, child = j, p = prop, flag=flag, nrep=nrep)
    }else {
      if (j == "F"){flag=0} else {flag=1}
      pa1 <- names(which(sapply(apply(data[,pa], 2, unique), function(x) length(x)) == 2))
      pa2 <- setdiff(pa, pa1)
      data[,j] <- Dcrt_2pas(data, pa1, pa2, flag=flag, child=j, nrep)
    }
  }
  
  #re <- setdiff(names(layer1), same_pa)
  #data[,re] <- Dcrt_1pa(data, pa1="E", child = re, p = 0.4, flag = 1, nrep)
  
  # for (j in names(layer1)){
  #   pa <- names(which(W_mat[,j] != 0))
  #   data[,j] <- Dcrt_1pa(data, pa1=pa, child = j, p = runif(1,0.1,1), nrep)
  # }
  
  # Simulate nodes from third layer to be categorical variables
  for (j in layer2){
    pa <- which(W_mat[,j] != 0)
    if (length(pa) > 1) {
      pa1 <- names(pa)[which(length(unique(data[,pa])) == 2)]
      pa2 <- setdiff(names(pa), pa1)
      if (j == "H"){flag=0} else {flag=1}
      data[,j] <- Dcrt_2pas(data, pa1, pa2, flag=flag, child=j, nrep)
    } else {
      prop <- sample(c(0.3, 0.9),1)
      if (prop == 0.3) {flag = 1} else {flag = 0}
      data[,j] <- Dcrt_1pa(data, pa1=pa, child = j, p = prop, flag = flag, nrep)
    }
  }
  
  return(data)
}


##### Discretization Policy
CombDiscretised_old <- function(dtype, data){
  if(dtype == '2l'){
    n <- as.numeric(gsub("l", "", dtype))
    d <- sapply(as.data.frame(data), function(x) unique(x%%1==0))
    data[!d] <- apply(as.data.frame(data)[!d], 2, 
                      function(x) cut(x, breaks = quantile(x, c(0,0.5, 1)),
                                      labels=paste0("Q",1:n),
                                      include.lowest = TRUE))
    data <- data %>% mutate_if(sapply(data, is.numeric), as.factor) %>% mutate_if(sapply(data, is.character), as.factor)
    return(data)
  } 
  if (dtype == '4l'){
    n <- as.numeric(gsub("l", "", dtype))
    d <- sapply(as.data.frame(data), function(x) unique(x%%1==0))
    data[!d] <- apply(as.data.frame(data)[!d], 2, 
                      function(x) cut(x, breaks = quantile(x, c(0,0.25,0.5,0.75,1)),
                                      labels=paste0("Q",1:n),
                                      include.lowest = TRUE))
    data <- data %>% mutate_if(sapply(data, is.numeric), as.factor) %>% mutate_if(sapply(data, is.character), as.factor)
    return(data)
  } 
}

Discretised <- function(data, nodes, method, ncategory){
  for (i in nodes){
    data[,i] <- arules::discretize(data[,i], method=method, breaks=ncategory)
  }
  return(data)
}
