adjMat2ggdag=function(adjMat,node.names){
  # this function convert adjacency matrix to ggdag
  # adjMat is the adjacency matrix. The ij element of adjacency matrix denotes the edge from the i^th row to J^th column
  # names is the names of variables
  library(stringr)
  library(lavaan)
  library(dagitty)
  s=c()
  n=nrow(adjMat)
  for( i in 1:n){
    for (j in 1:n){
      if (adjMat[i,j] !=0){
        tem=paste0(node.names[j],"~",node.names[i])
        s=c(s,tem)
      }
    }
  }
  lav.s=lavaanify(s)
  lav.s=subset(lav.s,op=="~")
  g=lavaanToGraph(lav.s)
  return(g)
}