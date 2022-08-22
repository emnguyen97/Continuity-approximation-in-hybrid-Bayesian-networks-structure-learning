# Functions to discretised  continuous data
# dtype '2l' and '4l' mean discretise based on quantiles, into 2 and 4 levels respectively.
# dtype '2eq' and '4eq' mean discretise based on equal intervals, into 2 and 4 levels respectively.

ContDiscretised <- function(dtype, data){
  
  if(dtype == '2l'){
    n <- as.numeric(gsub("l", "", dtype))
    data <- apply(data, 2, function(x) cut(x, breaks = quantile(x, c(0,0.5, 1)),
                                              labels = 1:n,
                                              include.lowest = TRUE))
    data <- as.data.frame(data)
    data <- data %>% mutate_if(sapply(data, is.character), as.numeric) 
    data <- data %>% mutate_if(sapply(data, is.numeric), as.factor)
    return(data)
  } 
  
  
  if (dtype == '4l'){
    n <- as.numeric(gsub("l", "", dtype))
    data <- apply(data, 2, function(x) cut(x, breaks = quantile(x, c(0,0.25,0.5,0.75,1)),
                                           labels = 1:n,
                                           include.lowest = TRUE))
    data <- as.data.frame(data)
    data <- data %>% mutate_if(sapply(data, is.character), as.numeric) 
    data <- data %>% mutate_if(sapply(data, is.numeric), as.factor)
    return(data)
  } 
  if (dtype == '2eq'){
    n <- as.numeric(gsub("eq", "", dtype))
    data <- apply(data, 2, function(x) cut(x, breaks = c(min(x),(min(x)+max(x))/2,max(x)),
                                           labels = 1:n,
                                           include.lowest = TRUE))
    data <- as.data.frame(data)
    data <- data %>% mutate_if(sapply(data, is.character), as.numeric) 
    data <- data %>% mutate_if(sapply(data, is.numeric), as.factor)
    return(data)
  } 
  if (dtype == '4eq'){
    n <- as.numeric(gsub("eq", "", dtype))
    data <- apply(data, 2, function(x) cut(x, breaks = c(min(x), min(x)+(max(x)-min(x))/4, min(x)+(max(x)-min(x))/2,
                                                         min(x)+3/4*(max(x)-min(x)), max(x)),
                                           labels = 1:n,
                                           include.lowest = TRUE))
    data <- as.data.frame(data)
    data <- data %>% mutate_if(sapply(data, is.character), as.numeric) 
    data <- data %>% mutate_if(sapply(data, is.numeric), as.factor)
    return(data)
  }
}
