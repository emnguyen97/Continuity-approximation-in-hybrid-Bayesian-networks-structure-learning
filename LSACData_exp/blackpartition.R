blackpartition=function(blacklist,dat){
  # this function convert the blacklist dataframe to a matrix
  b=matrix(0,nrow = ncol(dat), ncol = ncol(dat))
  for(i in 1: nrow(blacklist)){
    rowind=which(colnames(dat)==blacklist$from[i])
    colind=which(colnames(dat)==blacklist$to[i])
    b[rowind,colind]=1
  }
  rownames(b)=colnames(dat)
  colnames(b)=colnames(dat)
  return(b)
}