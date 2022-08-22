# Plotting pMCMC results

library(bnlearn)
library(stringr)
library(lavaan)
library(dagitty)
library(ggdag)
library(BiDAG)
library(ggraph)
library(tidygraph)
library(igraph)

## Implement partition MCMC for LSAC data

vartype = "BGe" # or "BGe"
iterations<-1000000
options(scipen=999)



source('adjMat2ggdag.R') # convert adjacency matrix to ggdag
source('extract.bn.fit.R') # extract coefficients and corresponding std's of these coefficients
source('blackpartition.R') # convert blacklist dataframe to a matrix
load('nodes.RData') # node labels


# "bge" results
# Need to load partition mcmc results
partcat <- bge_result$partbge

# dat.dag is the data used to run P-MCMCM.
dat.dag = B.W2

# sorting the dags output of partition MCMC.
top_down <- order((unlist(partcat$trace)),decreasing = TRUE) 
dags <- partcat$traceadd$incidence[top_down]
adjMat <- unique(dags)[[2]]

coef.dat <- extract.bn.fit(adjMat, dat.dag , bdecat=TRUE)

gr <- dplyr::mutate(coef.dat, col=as.numeric(coef>0)) %>%  tidygraph::as_tbl_graph()

new_gr <- tbl_graph(nodes = nodes, edges = igraph::as_data_frame(gr, what="edges"))

plots<- ggraph(new_gr,layout = "star") +   
  geom_edge_diagonal(aes(edge_width = 0.7),end_cap = circle(0.5,'cm'),arrow = arrow(angle = 30, length = unit(2, 'mm'),ends = "last", type = "closed")) +
  geom_node_point(aes(size=20,fill=type),shape=21) +
  geom_node_text(aes(label=name),colour="black") +
  scale_size(range = c(1,20),guide = FALSE)+
  scale_edge_width(range = c(0,1)) + 
  scale_fill_manual(values=c("darkseagreen2", "antiquewhite1","red"))+
  theme_gray() + theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(), axis.title.x=element_blank(),axis.title.y=element_blank()) + 
  theme_dag(legend.position = "none")

ggsave("BW2Plot.pdf", plots, device="pdf",width = 7.5, height = 7.5)


# "bdecat" results
# Need to load partition mcmc results

# dat.dag is the data used to run P-MCMCM.
dat.dag = W2.f

# sorting the dags output of partition MCMC.
top_down <- order((unlist(partcat$trace)),decreasing = TRUE) 
dags <- partcat$traceadd$incidence[top_down]
adjMat <- unique(dags)[[2]]

coef.dat <- extract.bn.fit(adjMat, dat.dag , bdecat=TRUE)

gr <- dplyr::mutate(coef.dat, col=as.numeric(coef>0)) %>%  tidygraph::as_tbl_graph()

new_gr <- tbl_graph(nodes = nodes, edges = igraph::as_data_frame(gr, what="edges"))

plots<- ggraph(new_gr,layout = "star") +   
  geom_edge_diagonal(aes(edge_width = 0.7),end_cap = circle(0.5,'cm'),arrow = arrow(angle = 30, length = unit(2, 'mm'),ends = "last", type = "closed")) +
  geom_node_point(aes(size=20,fill=type),shape=21) +
  geom_node_text(aes(label=name),colour="black") +
  scale_size(range = c(1,20),guide = FALSE)+
  scale_edge_width(range = c(0,1)) + 
  scale_fill_manual(values=c("darkseagreen2", "antiquewhite1","red"))+
  theme_gray() + theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(), axis.title.x=element_blank(),axis.title.y=element_blank()) + 
  theme_dag(legend.position = "none")

ggsave("W2fPlot.pdf", plots, device="pdf",width = 7.5, height = 7.5)

