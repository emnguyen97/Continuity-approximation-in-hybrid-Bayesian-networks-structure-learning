# Run Partition MCMC from BiDAG package and return results
# type: 'Continuous' for continuous data; 'Combination', 'Combination_cd' for 2 combination data (discrete parents & continuous children + continuous parents & discrete children); 'Discrete' for discrete data.
# discretise: 'None' for no discretisation, '2l' and '4l' for discretising based on quantiles (2 and 4 levels), '2eq' and '4eq' for discretising based on equal intervals (2 and 4 levels).
# model: 'bge' score for continuous or 'bdecat' score for discrete.


partMCMC <- function(type, discretise, model, nvar, nsim){
  
  if (type == 'Continuous' & discretise == 'None' & model == 'bge'){
    # Results (in the form of freq table) for Cont data, no discretisation
    Cont_bge <- lapply(1:nsim, function(i) Inference("bge", datasets[[i]]))
    
    # Extracting frequency table for unique structures
    Freq_contbge <- lapply(1:nsim, function(i) Cont_bge[[i]][[1]])
    
    compDAGs_cont.bge <- lapply(1:nsim, 
                                function(i) compareDAGs(datmat, matrix(Freq_contbge[[i]]$StructureStr[[which.max(Freq_contbge[[i]]$DAGscores)]],nvar,nvar)))
    
    return(c(Cont_bge, compDAGs_cont.bge))
  }
  
  if (type == 'Continuous' & discretise == '2l' & model == 'bdecat'){
    # Results for Continuous data, Multinomial, discretisation 2l
    Cont.2l_bdecat <- lapply(1:nsim, function(i) Inference("bdecat", ContDiscretised('2l', datasets[[i]])))
    
    # Extracting frequency table for unique structures
    Freq_cont2lbdecat <- lapply(1:nsim, function(i) Cont.2l_bdecat[[i]][[1]])
    
    compDAGs_cont2l.bdecat <- lapply(1:nsim, 
                                     function(i) compareDAGs(datmat, matrix(Freq_cont2lbdecat[[i]]$StructureStr[[which.max(Freq_cont2lbdecat[[i]]$DAGscores)]],nvar,nvar)))
    
    return(c(Cont.2l_bdecat, compDAGs_cont2l.bdecat))
  }
  
  if (type == 'Continuous' & discretise == '4l' & model == 'bdecat'){
    # Results for Continuous data, Multinomial, discretisation 4l
    Cont.4l_bdecat <- lapply(1:nsim, function(i) Inference("bdecat", ContDiscretised('4l', datasets[[i]])))
    
    # Extracting frequency table for unique structures
    Freq_cont4lbdecat <- lapply(1:nsim, function(i) Cont.4l_bdecat[[i]][[1]])
    
    compDAGs_cont4l.bdecat <- lapply(1:nsim, 
                                     function(i) compareDAGs(datmat, matrix(Freq_cont4lbdecat[[i]]$StructureStr[[which.max(Freq_cont4lbdecat[[i]]$DAGscores)]],nvar,nvar)))
    
    return(c(Cont.4l_bdecat, compDAGs_cont4l.bdecat))
  }
  
  if (type == 'Continuous' & discretise == '2eq' & model == 'bdecat'){
    # Results for Continuous data, Multinomial, discretisation 2l equal space
    Cont.2eq_bdecat <- lapply(1:nsim, function(i) Inference("bdecat", ContDiscretised('2eq', datasets[[i]])))
    
    # Extracting frequency table for unique structures
    Freq_cont2eqbdecat <- lapply(1:nsim, function(i) Cont.2eq_bdecat[[i]][[1]])
    
    compDAGs_cont2l.bdecat <- lapply(1:nsim, 
                                     function(i) compareDAGs(datmat, matrix(Freq_cont2eqbdecat[[i]]$StructureStr[[which.max(Freq_cont2eqbdecat[[i]]$DAGscores)]],nvar,nvar)))
    
    return(c(Cont.2eq_bdecat, compDAGs_cont2l.bdecat))
  }
  
  if (type == 'Continuous' & discretise == '4eq' & model == 'bdecat'){
    # Results for Continuous data, Multinomial, discretisation 4l equal space
    Cont.4eq_bdecat <- lapply(1:nsim, function(i) Inference("bdecat", ContDiscretised('4eq', datasets[[i]])))
    
    # Extracting frequency table for unique structures
    Freq_cont4eqbdecat <- lapply(1:nsim, function(i) Cont.4eq_bdecat[[i]][[1]])
    
    compDAGs_cont4eq.bdecat <- lapply(1:nsim, 
                                      function(i) compareDAGs(datmat, matrix(Freq_cont4eqbdecat[[i]]$StructureStr[[which.max(Freq_cont4eqbdecat[[i]]$DAGscores)]],nvar,nvar)))
    
    return(c(Cont.4eq_bdecat, compDAGs_cont4eq.bdecat))
  }
  
  if (type == 'Continuous' & discretise == '2l' & model == 'bge'){
    # Transform data for inference
    data2l <- lapply(1:nsim, 
                     function(i) ContDiscretised('2l', datasets[[i]]) %>% mutate_if(sapply(ContDiscretised('2l', datasets[[i]]), is.factor), as.numeric))
    
    # Results for Continuous data, Gaussian, discretisation 2l 
    Cont.2l_bge <- lapply(1:nsim, function(i) Inference("bge", data2l[[i]]))
    
    # Extracting frequency table for unique structures
    Freq_cont2lbge <- lapply(1:nsim, function(i) Cont.2l_bge[[i]][[1]])
    
    compDAGs_cont2l.bge <- lapply(1:nsim, 
                                  function(i) compareDAGs(datmat, matrix(Freq_cont2lbge[[i]]$StructureStr[[which.max(Freq_cont2lbge[[i]]$DAGscores)]],nvar,nvar)))
    
    return(c(Cont.2l_bge, compDAGs_cont2l.bge))
  }
  
  if (type == 'Continuous' & discretise == '4l' & model == 'bge'){
    data4l <- lapply(1:nsim, 
                     function(i) ContDiscretised('4l', datasets[[i]]) %>% mutate_if(sapply(ContDiscretised('4l', datasets[[i]]), is.factor), as.numeric))
    
    # Results for Continuous data, Gaussian, discretisation 4l 
    Cont.4l_bge <- lapply(1:nsim, function(i) Inference("bge", data4l[[i]]))
    
    # Extracting frequency table for unique structures
    Freq_cont4lbge <- lapply(1:nsim, function(i) Cont.4l_bge[[i]][[1]])
    
    compDAGs_cont4l.bge <- lapply(1:nsim, 
                                  function(i) compareDAGs(datmat, matrix(Freq_cont4lbge[[i]]$StructureStr[[which.max(Freq_cont4lbge[[i]]$DAGscores)]],nvar,nvar)))
    
    return(c(Cont.4l_bge, compDAGs_cont4l.bge))
  }
  
  if (type == 'Continuous' & discretise == '2eq' & model == 'bge'){
    data2eq <- lapply(1:nsim, 
                      function(i) ContDiscretised('2eq', datasets[[i]]) %>% mutate_if(sapply(ContDiscretised('2eq', datasets[[i]]), is.factor), as.numeric))
    
    # Results for Continuous data, Gaussian, discretisation 2l equal space
    Cont.2eq_bge <- lapply(1:nsim, function(i) Inference("bge", data2eq[[i]]))
    
    # Extracting frequency table for unique structures
    Freq_cont2eqbge <- lapply(1:nsim, function(i) Cont.2eq_bge[[i]][[1]])
    
    compDAGs_cont2eq.bge <- lapply(1:nsim, 
                                   function(i) compareDAGs(datmat, matrix(Freq_cont2eqbge[[i]]$StructureStr[[which.max(Freq_cont2eqbge[[i]]$DAGscores)]],nvar,nvar)))
    
    return(c(Cont.2eq_bge, compDAGs_cont2eq.bge))
  }
  
  if (type == 'Continuous' & discretise == '4eq' & model == 'bge'){
    data4eq <- lapply(1:nsim, 
                      function(i) ContDiscretised('4eq', datasets[[i]]) %>% mutate_if(sapply(ContDiscretised('4eq', datasets[[i]]), is.factor), as.numeric))
    
    # Results for Continuous data, Gaussian, discretisation 4l equal space
    Cont.4eq_bge <- lapply(1:nsim, function(i) Inference("bge", data4eq[[i]]))
    
    # Extracting frequency table for unique structures
    Freq_cont4eqbge <- lapply(1:nsim, function(i) Cont.4eq_bge[[i]][[1]])
    
    compDAGs_cont4eq.bge <- lapply(1:nsim, 
                                   function(i) compareDAGs(datmat, matrix(Freq_cont4eqbge[[i]]$StructureStr[[which.max(Freq_cont4eqbge[[i]]$DAGscores)]],nvar,nvar)))
    
    return(c(Cont.4eq_bge, compDAGs_cont4eq.bge))
  }
  
  if (type %in% c('Combination', 'Combination_cd') & discretise == 'None' & model == 'bge'){
    # Results for Combination data, Gaussian, no discretisation
    if (type == 'Combination') {
      Comb_bge <- lapply(1:nsim, function(i) Inference("bge", Comb.dats[[i]]))
    } else if (type == 'Combination_cd'){
      Comb_bge <- lapply(1:nsim, function(i) Inference("bge", Comb_cd.dats[[i]]))
    }
    
    # Extracting frequency table for unique structures
    Freq_combbge <- lapply(1:nsim, function(i) Comb_bge[[i]][[1]])
    
    compDAGs_comb.bge <- lapply(1:nsim, 
                                function(i) compareDAGs(datmat, matrix(Freq_combbge[[i]]$StructureStr[[which.max(Freq_combbge[[i]]$DAGscores)]],nvar,nvar)))
    
    return(c(Comb_bge, compDAGs_comb.bge))
  }
  
  if (type %in% c('Combination', 'Combination_cd') & discretise == '2l' & model == 'bdecat'){
    # Results for Combination data, Multinomial, discretisation 2l
    if (type == 'Combination') {
      Comb.2l_bdecat <- lapply(1:nsim, function(i) Inference("bdecat", CombDiscretised('2l', Comb.dats[[i]])))
    } else if (type == 'Combination_cd'){
      Comb.2l_bdecat <- lapply(1:nsim, function(i) Inference("bdecat", CombDiscretised('2l', Comb_cd.dats[[i]])))
    }
    
    # Extracting frequency table for unique structures
    Freq_comb2lbdecat <- lapply(1:nsim, function(i) Comb.2l_bdecat[[i]][[1]])
    
    compDAGs_comb2l.bdecat <- lapply(1:nsim, 
                                     function(i) compareDAGs(datmat, matrix(Freq_comb2lbdecat[[i]]$StructureStr[[which.max(Freq_comb2lbdecat[[i]]$DAGscores)]],nvar,nvar)))
    
    return(c(Comb.2l_bdecat, compDAGs_comb2l.bdecat))
  }
  
  if (type %in% c('Combination', 'Combination_cd') & discretise == '4l' & model == 'bdecat'){
    # Results for Combination data, Multinomial, discretisation 4l
    if (type == 'Combination') {
      Comb.4l_bdecat <- lapply(1:nsim, function(i) Inference("bdecat", CombDiscretised('4l', Comb.dats[[i]])))
    } else if (type == 'Combination_cd'){
      Comb.4l_bdecat <- lapply(1:nsim, function(i) Inference("bdecat", CombDiscretised('4l', Comb_cd.dats[[i]])))
    }
    
    # Extracting frequency table for unique structures
    Freq_comb4lbdecat <- lapply(1:nsim, function(i) Comb.4l_bdecat[[i]][[1]])
    
    compDAGs_comb4l.bdecat <- lapply(1:nsim, 
                                     function(i) compareDAGs(datmat, matrix(Freq_comb4lbdecat[[i]]$StructureStr[[which.max(Freq_comb4lbdecat[[i]]$DAGscores)]],nvar,nvar)))
    
    return(c(Comb.4l_bdecat, compDAGs_comb4l.bdecat))
  }
  
  if (type %in% c('Combination', 'Combination_cd') & discretise == '2eq' & model == 'bdecat'){
    # Results for Combination data, Multinomial, discretisation 2l equal space
    if (type == 'Combination') {
      Comb.2eq_bdecat <- lapply(1:nsim, function(i) Inference("bdecat", CombDiscretised('2eq', Comb.dats[[i]])))
    } else if (type == 'Combination_cd'){
      Comb.2eq_bdecat <- lapply(1:nsim, function(i) Inference("bdecat", CombDiscretised('2eq', Comb_cd.dats[[i]])))
    }
    
    # Extracting frequency table for unique structures
    Freq_comb2eqbdecat <- lapply(1:nsim, function(i) Comb.2eq_bdecat[[i]][[1]])
    
    compDAGs_comb2eq.bdecat <- lapply(1:nsim, 
                                      function(i) compareDAGs(datmat, matrix(Freq_comb2eqbdecat[[i]]$StructureStr[[which.max(Freq_comb2eqbdecat[[i]]$DAGscores)]],nvar,nvar)))
    
    return(c(Comb.2eq_bdecat, compDAGs_comb2eq.bdecat))
  }
  
  if (type %in% c('Combination', 'Combination_cd') & discretise == '4eq' & model == 'bdecat'){
    # Results for Combination data, Multinomial, discretisation 4l equal space
    if (type == 'Combination') {
      Comb.4eq_bdecat <- lapply(1:nsim, function(i) Inference("bdecat", CombDiscretised('4eq', Comb.dats[[i]])))
    } else if (type == 'Combination_cd'){
      Comb.4eq_bdecat <- lapply(1:nsim, function(i) Inference("bdecat", CombDiscretised('4eq', Comb_cd.dats[[i]])))
    }
    
    # Extracting frequency table for unique structures
    Freq_comb4eqbdecat <- lapply(1:nsim, function(i) Comb.4eq_bdecat[[i]][[1]])
    
    compDAGs_comb4eq.bdecat <- lapply(1:nsim, 
                                      function(i) compareDAGs(datmat, matrix(Freq_comb4eqbdecat[[i]]$StructureStr[[which.max(Freq_comb4eqbdecat[[i]]$DAGscores)]],nvar,nvar)))
    
    return(c(Comb.4eq_bdecat, compDAGs_comb4eq.bdecat))
  }
  
  if (type %in% c('Combination', 'Combination_cd') & discretise == '2l' & model == 'bge'){
    if (type == 'Combination') {
      Comb2l <- lapply(1:nsim, function(i) CombDiscretised('2l', Comb.dats[[i]]) %>% mutate_if(sapply(CombDiscretised('2l', Comb.dats[[i]]), is.factor), as.numeric))
    } else if (type == 'Combination_cd'){
      Comb2l <- lapply(1:nsim, function(i) CombDiscretised('2l', Comb_cd.dats[[i]]) %>% mutate_if(sapply(CombDiscretised('2l', Comb_cd.dats[[i]]), is.factor), as.numeric))
    }
    
    # Results for Combination data, Gaussian, discretisation 2l
    Comb.2l_bge <- lapply(1:nsim, function(i) Inference("bge", Comb2l[[i]]))
    
    # Extracting frequency table for unique structures
    Freq_comb2lbge <- lapply(1:nsim, function(i) Comb.2l_bge[[i]][[1]])
    
    compDAGs_comb2l.bge <- lapply(1:nsim, 
                                  function(i) compareDAGs(datmat, matrix(Freq_comb2lbge[[i]]$StructureStr[[which.max(Freq_comb2lbge[[i]]$DAGscores)]],nvar,nvar)))
    
    return(c(Comb.2l_bge, compDAGs_comb2l.bge))
  }
  
  if (type %in% c('Combination', 'Combination_cd') & discretise == '4l' & model == 'bge'){
    if (type == 'Combination') {
      Comb4l <- lapply(1:nsim, function(i) CombDiscretised('4l', Comb.dats[[i]]) %>% mutate_if(sapply(CombDiscretised('4l', Comb.dats[[i]]), is.factor), as.numeric))
    } else if (type == 'Combination_cd'){
      Comb4l <- lapply(1:nsim, function(i) CombDiscretised('4l', Comb_cd.dats[[i]]) %>% mutate_if(sapply(CombDiscretised('4l', Comb_cd.dats[[i]]), is.factor), as.numeric))
    }
    
    # Results for Combination data, Gaussian, discretisation 4l
    Comb.4l_bge <- lapply(1:nsim, function(i) Inference("bge", Comb4l[[i]]))
    
    # Extracting frequency table for unique structures
    Freq_comb4lbge <- lapply(1:nsim, function(i) Comb.4l_bge[[i]][[1]])
    
    compDAGs_comb4l.bge <- lapply(1:nsim, 
                                  function(i) compareDAGs(datmat, matrix(Freq_comb4lbge[[i]]$StructureStr[[which.max(Freq_comb4lbge[[i]]$DAGscores)]],nvar,nvar)))
    
    return(c(Comb.4l_bge, compDAGs_comb4l.bge))
  }
  
  if (type %in% c('Combination', 'Combination_cd') & discretise == '2eq' & model == 'bge'){
    if (type == 'Combination') {
      Comb2eq <- lapply(1:nsim, function(i) CombDiscretised('2eq', Comb.dats[[i]]) %>% mutate_if(sapply(CombDiscretised('2eq', Comb.dats[[i]]), is.factor), as.numeric))
    } else if (type == 'Combination_cd'){
      Comb2eq <- lapply(1:nsim, function(i) CombDiscretised('2eq', Comb_cd.dats[[i]]) %>% mutate_if(sapply(CombDiscretised('2eq', Comb_cd.dats[[i]]), is.factor), as.numeric))
    }
    
    # Results for Combination data, Gaussian, discretisation 2l equal space
    Comb.2eq_bge <- lapply(1:nsim, function(i) Inference("bge", Comb2eq[[i]]))
    
    # Extracting frequency table for unique structures
    Freq_comb2eqbge <- lapply(1:nsim, function(i) Comb.2eq_bge[[i]][[1]])
    
    compDAGs_comb2eq.bge <- lapply(1:nsim, 
                                   function(i) compareDAGs(datmat, matrix(Freq_comb2eqbge[[i]]$StructureStr[[which.max(Freq_comb2eqbge[[i]]$DAGscores)]],nvar,nvar)))
    
    return(c(Comb.2eq_bge, compDAGs_comb2eq.bge))
  }
  
  if (type %in% c('Combination', 'Combination_cd') & discretise == '4eq' & model == 'bge'){
    if (type == 'Combination') {
      Comb4eq <- lapply(1:nsim, function(i) CombDiscretised('4eq', Comb.dats[[i]]) %>% mutate_if(sapply(CombDiscretised('4eq', Comb.dats[[i]]), is.factor), as.numeric))
    } else if (type == 'Combination_cd'){
      Comb4eq <- lapply(1:nsim, function(i) CombDiscretised('4eq', Comb_cd.dats[[i]]) %>% mutate_if(sapply(CombDiscretised('4eq', Comb_cd.dats[[i]]), is.factor), as.numeric))
    }
    
    # Results for Combination data, Gaussian, discretisation 4l equal space
    Comb.4eq_bge <- lapply(1:nsim, function(i) Inference("bge", Comb4eq[[i]]))
    
    # Extracting frequency table for unique structures
    Freq_comb4eqbge <- lapply(1:nsim, function(i) Comb.4eq_bge[[i]][[1]])
    
    compDAGs_comb4eq.bge <- lapply(1:nsim, 
                                   function(i) compareDAGs(datmat, matrix(Freq_comb4eqbge[[i]]$StructureStr[[which.max(Freq_comb4eqbge[[i]]$DAGscores)]],nvar,nvar)))
    
    return(c(Comb.4eq_bge, compDAGs_comb4eq.bge))
  }
  
  if (type == 'Categorical' & discretise == 'None' & model == 'bdecat'){
    # Results for Discrete data, no discretisation, Multinomial
    Dcrt_bdecat <- lapply(1:nsim, function(i) Inference("bdecat", dcrt.dats[[i]]))
    
    # Extracting frequency table for unique structures
    Freq_dcrtbdecat <- lapply(1:nsim, function(i) Dcrt_bdecat[[i]][[1]])
    
    compDAGs_Dcrt.bdecat <- lapply(1:nsim, 
                                   function(i) compareDAGs(datmat, matrix(Freq_dcrtbdecat[[i]]$StructureStr[[which.max(Freq_dcrtbdecat[[i]]$DAGscores)]],nvar,nvar)))
    
    return(c(Dcrt_bdecat, compDAGs_Dcrt.bdecat))
  }
  
  if (type == 'Categorical' & discretise == 'None' & model == 'bge'){
    dcrt.bge.dats <- lapply(1:nsim, 
                            function(i) dcrt.dats[[i]] %>% mutate_if(sapply(dcrt.dats[[i]], is.factor), as.numeric))
    # Results for Discrete data, no discretisation, Gaussian
    Dcrt_bge <- lapply(1:nsim, function(i) Inference("bge", dcrt.bge.dats[[i]]))
    
    # Extracting frequency table for unique structures
    Freq_dcrtbge <- lapply(1:nsim, function(i) Dcrt_bge[[i]][[1]])
    
    compDAGs_Dcrt.bge <- lapply(1:nsim, 
                                function(i) compareDAGs(datmat, matrix(Freq_dcrtbge[[i]]$StructureStr[[which.max(Freq_dcrtbge[[i]]$DAGscores)]],nvar,nvar)))
    
    return(c(Dcrt_bge, compDAGs_Dcrt.bge))
  }
  
}
