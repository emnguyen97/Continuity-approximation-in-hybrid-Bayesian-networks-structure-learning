library(dplyr)

extract.bn.fit =function(DAG,dat,bdecat = FALSE){
  # this function extract coefficients and corresponding std's of these coefficients.
  
  bn.temp=bnlearn::empty.graph(colnames(DAG))
  amat(bn.temp)=DAG
  
  if (bdecat == TRUE){
  fitted <- bn.fit(bn.temp,data = dat %>% mutate_if(sapply(dat, is.factor), as.numeric))
  } else {
    fitted = bn.fit(x=bn.temp,data = dat)
  }
  
  coef.mat=coef(fitted)
  coef.dat=data.frame(from=0,to = 0,coef = 0)
  cnt=0
  for(i in 1:length(coef.mat)){
    con.dist=coef.mat[[i]] # the conditional distribution
    if(length(names(con.dist))>1){
      for(j in 2:length(names(con.dist))){
        cnt=cnt+1
        coef.dat[cnt,"from"] = names(con.dist)[j]
        coef.dat[cnt,"to"] = names(coef.mat)[i]
        coef.dat[cnt,"coef"] = con.dist[j]
      }
    }
  }
  coef.dat=coef.dat[with(coef.dat, order(from,to)),] # rearrange the dataframe according to the column 'from' and 'to'
  
  # extract the std
  std.dat=data.frame(from=0,to = 0,std = 0)
  col.name=colnames(DAG)
  
  for(i in 1:dim(DAG)[2]){
    if(sum(DAG[,i])>0){
      ind = as.vector(which(DAG[,i]>0))
      pred=col.name[ind][1]
      if(length(ind)>1){
        for(j in 1:(length(ind)-1)){
          pred= paste(pred,col.name[ind][j+1],sep="+")
        }
      }
      
      model = lm( as.formula(paste0(col.name[i],"~",pred)), data=dat )
      o.std=summary(model)$coefficients[,2]
      
      for(j in 1:length(ind)){
        std.dat = rbind(std.dat,c(col.name[ind][j],col.name[i],o.std[j+1]))
      }
    }
  }
  std.dat=std.dat[-1,]
  
  
  para = dplyr::left_join(coef.dat,std.dat)
  return(para)
}