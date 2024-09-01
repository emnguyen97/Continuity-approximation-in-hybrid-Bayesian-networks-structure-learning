library(tidyr)
library(tibble)
library(dplyr)
beta= seq(-3,3,length.out=100) #50
sb=1
sa=1
n=100000 # for approximation of p11 or p10

# for infinity sample size

# S.cc
postratio.bge.cc= function(beta,sa,sb){
  r= 1/2 * log((sb^2+ beta^2*sa^2)/(sb^2))
  return(r)
}

mc.p11=function(sb,sa,beta){ # use numerical integration to compute p11
  ga = sb/sa/beta
  m=100
  x=seq(0,m,length.out=n)
  s=numeric(n)
  for(i in 1:n){
    s[i]= 1/sqrt(2*pi) *exp(-x[i]^2/2)*pnorm(x[i]/ga)
  }
  f=sum(s)*m/n
  return(f)
}

postratio.bde.cc= function(beta,sa,sb){
  p= mc.p11(sb,sa,beta)
  r= log(4)+ 2*p*log(p)+ (1-2*p)*log(0.5-p) 
  return(r)
}

Scc.inf.table = tibble() # store all the results for S_cc scenario

for(i in 1:length(beta)){
  r= postratio.bge.cc(beta[i],sa,sb) 
  tem.table = tibble(scenario = "cc",beta = beta[i],type="RAG",r = r)
  Scc.inf.table= rbind(Scc.inf.table,tem.table) # BGe
  
  r=postratio.bde.cc(beta[i],sa,sb)
  tem.table = tibble(scenario = "cc",beta = beta[i],type="DISC",r = r)
  Scc.inf.table= rbind(Scc.inf.table,tem.table) # BDe
  cat(postratio.bde.cc(beta[i],sa,sb),'\n')
}

library(ggplot2)
p = ggplot(Scc.inf.table,aes(x = beta,y = r, linetype = type, col=type))
p = p+ geom_line(lwd = 0.7) + scale_color_manual(values = c("#006600", "#000066")) 
p = p+labs(linetype='Strategy',color='Strategy') +xlab(expression(beta)) +ylab("Posterior ratio")
p = p + theme_light() 
p = p + theme(legend.position = c(0.5, 0.86))
p
# + geom_hline(yintercept=log(2) , color= 'blue',linetype = "dashed")
ggsave('2nodes_exp/Postratio.inf.scc.negbeta.pdf',plot = p,width = 6,height=4)


# Scd
mc.Ex2 = function(beta,sa){
  m=200
  x=seq(-m/2,m/2,length.out=n)
  s=numeric(n)
  for(i in 1:n){
    s[i]= 1/sqrt(2*pi) * exp(-x[i]^2/2)*exp(x[i]*sa*beta)/(1+ exp(x[i]*sa*beta))
  }
  f=sum(s)*m/n
  return(f)
}
mc.S12 = function(beta,sa){
  m=200
  x=seq(-m/2,m/2,length.out=n)
  s=numeric(n)
  for(i in 1:n){
    s[i]= x[i]*sa/sqrt(2*pi) * exp(-x[i]^2/2)*exp(x[i]*sa*beta)/(1+ exp(x[i]*sa*beta))
  }
  f=sum(s)*m/n
  return(f)
}
postratio.bge.cd= function(beta,sa){
  Sigma=matrix(0,2,2)
  Sigma[1,1]=sa^2
  Sigma[1,2]=mc.S12(beta,sa)
  Sigma[2,1]=Sigma[1,2]
  # E(X_2)
  ex2= mc.Ex2(beta,sa)
  Sigma[2,2] = ex2-ex2^2
  r= 1/2 * log((Sigma[1,1]*Sigma[2,2])/(Sigma[1,1]*Sigma[2,2] - Sigma[1,2]^2))
  return(r)
}

mc.p11.cd=function(sa,beta){ # use numerical integration to compute p11
  m=100
  x=seq(0,m,length.out=n)
  s=numeric(n)
  for(i in 1:n){
    s[i]= 1/sqrt(2*pi) *exp(-x[i]^2/2)*exp(x[i]*beta*sa)/(1+exp(x[i]*beta*sa))
  }
  f=sum(s)*m/n
  return(f)
}

postratio.bde.cd= function(beta,sa){
  p= mc.p11.cd(sa,beta)
  r= log(4)+ 2*p*log(p)+ (1-2*p)*log(0.5-p) 
  return(r)
}

Scd.inf.table = tibble() # store all the results for S_cd scenario

for(i in 1:length(beta)){
  r= postratio.bge.cd(beta[i],sa) 
  tem.table = tibble(scenario = "cd",beta = beta[i],type="RAG",r = r)
  Scd.inf.table= rbind(Scd.inf.table,tem.table) # BGe
  
  r=postratio.bde.cd(beta[i],sa)
  tem.table = tibble(scenario = "cd",beta = beta[i],type="DISC",r = r)
  Scd.inf.table= rbind(Scd.inf.table,tem.table) # BDe
  cat(postratio.bde.cc(beta[i],sa,sb),'\n')
}

library(ggplot2)
p = ggplot(Scd.inf.table,aes(x = beta,y = r, linetype = type, col=type))
p = p+ geom_line(lwd = 0.7) + scale_color_manual(values = c("#006600", "#000066")) 
p = p+labs(linetype='Strategy',color='Strategy') +xlab(expression(beta)) +ylab("Posterior ratio")
p = p + theme_light() 
p = p + theme(legend.position = c(0.5, 0.86))
p
# + geom_hline(yintercept=log(2) , color= 'blue',linetype = "dashed")
ggsave('2nodes_exp/Postratio.inf.scd.negbeta.pdf',plot = p,width = 6,height=4)

# Sdc
postratio.bge.dc= function(beta,sb){
  r= 1/2 * log(1+ (beta^2)/(4*sb^2))
  return(r)
}
mc.p11.dc=function(sb,beta){ # use numerical integration to compute p11
  m=100
  x=seq(-beta/2/sb,m,length.out=n)
  s=numeric(n)
  for(i in 1:n){
    s[i]= 1/sqrt(2*pi) *exp(-x[i]^2/2)
  }
  f=1/2*sum(s)*(m+beta/2/sb)/n
  return(f)
}

postratio.bde.dc= function(beta,sb){
  p= mc.p11.dc(sb,beta)
  r= log(4)+ 2*p*log(p)+ (1-2*p)*log(0.5-p) 
  return(r)
}

Sdc.inf.table = tibble() # store all the results for S_cd scenario

for(i in 1:length(beta)){
  r= postratio.bge.dc(beta[i],sb) 
  tem.table = tibble(scenario = "dc",beta = beta[i],type="RAG",r = r)
  Sdc.inf.table= rbind(Sdc.inf.table,tem.table) # BGe
  
  r=postratio.bde.dc(beta[i],sb)
  tem.table = tibble(scenario = "dc",beta = beta[i],type="DISC",r = r)
  Sdc.inf.table= rbind(Sdc.inf.table,tem.table) # BDe
  cat(r,'\n')
}


p = ggplot(Sdc.inf.table,aes(x = beta,y = r, linetype = type, col=type))
p = p+ geom_line(lwd = 0.7) + scale_color_manual(values = c("#006600", "#000066"))
p = p+labs(linetype='Strategy',color='Strategy') +xlab(expression(beta)) +ylab("Posterior ratio")
p = p + theme_light() 
p = p + theme(legend.position = c(0.5, 0.86))
p
# + geom_hline(yintercept=log(2) , color= 'blue',linetype = "dashed")
ggsave('2nodes_exp/Postratio.inf.sdc.negbeta.pdf',plot = p,width = 6,height=4)



# Sdd
p=1/2
beta= seq(-0.999,0.999,length.out=100) #50
postratio.bge.dd= function(beta,p){
  r= 1/2 * log( (1- (2*p-1)^2* beta^2)/(1-beta^2))
  return(r)
}

postratio.bde.dd= function(beta,p){
  p11= p*(0.5+0.5*beta)
  p10 = p*(0.5-0.5*beta)
  p01 = (1-p)*(0.5-0.5*beta)
  p00= (1-p)*(0.5+ 0.5*beta)
  p1. = p11+p10
  p.1 = p11+p01
  p0. = p01+p00
  p.0 = p10+p00
  r=p11* log(p11/(p1.* p.1)) + p10*log(p10/(p1.*p.0)) + p01* log(p01/(p0.*p.1)) + p00* log(p00/(p0.*p.0))
  return(r)
}

Sdd.inf.table = tibble() # store all the results for S_dd scenario

for(i in 1:length(beta)){
  r= postratio.bge.dd(beta[i],p) 
  tem.table = tibble(scenario = "dd",beta = beta[i],type="RAG",r = r)
  Sdd.inf.table= rbind(Sdd.inf.table,tem.table) # BGe
  
  r=postratio.bde.dd(beta[i],p)
  tem.table = tibble(scenario = "dd",beta = beta[i],type="DISC",r = r)
  Sdd.inf.table= rbind(Sdd.inf.table,tem.table) # BDe
  cat(r,'\n')
}
p = ggplot(Sdd.inf.table,aes(x = beta,y = r, linetype = type, col=type))
p = p+ geom_line(lwd = 0.7) + scale_color_manual(values = c("#006600", "#000066")) 
p = p+labs(linetype='Strategy',color='Strategy') +xlab(expression(beta)) +ylab("Posterior ratio")
p = p + theme_light() 
p = p + theme(legend.position = c(0.5, 0.86))
p
# + geom_hline(yintercept=log(2) , color= 'blue',linetype = "dashed")
ggsave('2nodes_exp/Postratio.inf.sdd.negbeta.pdf',plot = p,width =6,height=4)



