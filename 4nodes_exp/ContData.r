# Function to generate Continuous data for 4 nodes

ContData <- function(data, nrep, nvar){
  data$A <- rnorm(nrep,mean=-3, sd=1) #  has no parents
  data$C <- rnorm(nrep,mean=6, sd=1) #  has no parents
  data$B <- data$A * 1.5 + data$C * 3 + rnorm(nrep) 
  data$D <- data$A * 2 + rnorm(nrep) + data$C * 1.5
  return(data)
}





## MOre general function to generate continuous data based on adj matrix

# Function to simulate nodes with no parents
NodesnoParents <- function(adjmat, dat, mu_vec, sd_vec, nrep, nvar){
  # Obtain nodes with no parents
  np <- which(colSums(adjmat) == 0)
  
  # Simulating nodes 
  for (i in 1:length(np)) {
    dat[,np[i]] <- rnorm(nrep, mean=mu_vec[np[i]], sd=sd_vec[np[i]])
  }
  return(dat)
}

# Function to simulate nodes with parents
NodeswParents <- function(n.adjmat, adjmat, W_mat, mu_vec, sd_vec, dat, nrep){
  # Finding nodes without parents
  #n.np <- colnames(n.adjmat[,which(colSums(n.adjmat) == 0)])
  n.np <- names(which(colSums(n.adjmat) == 0))
  # Get index from original adjmat
  n.indx <- sapply(1:length(n.np), function(j) grep(n.np[j], colnames(adjmat)))
  
  # Iterate through their original parent nodes and index
  for (i in 1:length(n.indx)) {
    p <- which(adjmat[,n.indx[i]]==1)
    
    if (length(p)>1){#having multiple parents
      # Getting parents names
      parents <- colnames(adjmat[,p])
      # Parents original index
      parents.indx <- sapply(1:length(parents), function(j) grep(parents[j], colnames(adjmat)))
    } else if (length(p)==1){#having only 1 parent
      parents <- colnames(adjmat)[p]
      parents.indx <- grep(parents, colnames(adjmat))
    }
    
  # Parents coefficients
    parents.coef <- W_mat[parents.indx, n.indx[i]]
    
  # Simulate coefficients
    lmat <- as.data.frame(sapply(1:length(p), function(x) parents.coef[x]*dat[,parents.indx[x]]))
    dat[,n.indx[i]] <- mu_vec[n.indx[i]] + apply(lmat, 1, sum) + rnorm(nrep, sd=sd_vec[n.indx[i]])
  }
  return(dat)
}


CData <- function(W_mat, mu_vec, sd_vec, nrep, nvar){
  ### Inputs ###
  ## W_mat: Weight matrix with nodes in alphabetical order
  ## mu_vec: mean vector for each nodes, in the order of adjmat
  ## nrep, nvar: number of observation, number of variables
  
  # Get Adjacency matrix 
  adjmat <- ifelse((as.matrix(W_mat)!= 0)=="TRUE",1,0)
  
  # Creating empty matrix
  dat <- matrix(NA,nrow=nrep,ncol = nvar)
  dat <- as_tibble(dat)
  colnames(dat) <- colnames(adjmat)
  
  # Simulating nodes without parents first
  dat <- NodesnoParents(adjmat, dat, mu_vec, sd_vec, nrep, nvar)
  
  # Nodes with no parents
  np <- which(colSums(adjmat) == 0)
  
  # Nodes with parents
  wp <- colnames(adjmat[,which(colSums(adjmat) != 0)])
  
  # Adj Matrix after removing nodes with parents
  mat <- n.adjmat <- adjmat[-np,-np]
  
  while (sum(mat) > 0) {
    
    # Find the next set of outpoints
    #n.np <- colnames(mat[,which(colSums(mat) == 0)])
    n.np <- names(which(colSums(mat) == 0))

    # Simulating those nodes with parents
    dat <- NodeswParents(mat, adjmat, W_mat, mu_vec, sd_vec, dat, nrep)
    
    # simulating remaining nodes
    re.nodes <- wp <- setdiff(wp,n.np)
    
    re.indx <- sapply(1:length(re.nodes), function(j) grep(re.nodes[j], colnames(adjmat)))
    
    # Adj Matrix after removing all simulated nodes
    r.indx <- sapply(1:length(re.nodes), function(j) grep(re.nodes[j], colnames(n.adjmat)))
    r.adjmat <- mat[-unique(col(mat))[-r.indx], -unique(col(mat))[-r.indx]]
    
    mat <- r.adjmat
  } 
  
  # Iterate through their original parent nodes and index
  for (i in 1:length(re.indx)) {
    p <- which(adjmat[,re.indx[i]]==1)
    
    if (length(p)>1){#having multiple parents
      # Getting parents names
      parents <- colnames(adjmat[,p])
      # Parents original index
      parents.indx <- sapply(1:length(parents), function(j) grep(parents[j], colnames(adjmat)))
    } else if (length(p)==1){#having only 1 parent
      parents <- colnames(adjmat)[p]
      parents.indx <- grep(parents, colnames(adjmat))
    }

    # Parents coefficients
    parents.coef <- W_mat[parents.indx, re.indx[i]]
    
    # Simulate coefficients
    lmat <- as.data.frame(sapply(1:length(p), function(x) parents.coef[x]*dat[,parents.indx[x]]))
    dat[,re.indx[i]] <- mu_vec[re.indx[i]] + apply(lmat, 1, sum) + rnorm(nrep, sd=sd_vec[re.indx[i]])
  }
  return(dat)
}
