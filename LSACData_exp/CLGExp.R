# LSAC data
load("RAG/LSACData_exp/imputed.dat.all.RData")
# Load the necessary functions
source('RAG/structureMCMC/combinations.R')
source('RAG/structureMCMC/scoretables.R')
source('RAG/structureMCMC/structurefns.R')
source('RAG/structureMCMC/samplefns.R')
source('RAG/structureMCMC/structureMCMC.R')
source('RAG/LSACData_exp/blackpartition.R')

library(BiDAG)
library(bnlearn)
library(tibble)
library(ggdag)
library(ggraph)
library(dagitty)
library(bnlearn)
library(reshape2)
library(igraph)
library(gRain)
library(dplyr)
library(tidygraph)

## Selecting Wave 2 for analysis

B.W2 <- imputed.dat.k.all[[2]]

B.W2 <- B.W2 %>%
  select(SX,BMI,SE,AC,INC,FS,FH,ME1,FE1,BM1,BM2,RP1,DP1,FV,HF,HSD,SL,OD,GW,BWZ)

data <- as.data.frame(B.W2)

# to do it for some names in a vector named 'col_names'
col_names <- names(data)[!names(B.W2) %in% c("BMI","INC","SE","BM1","BM2","BWZ")]
data[col_names] <- lapply(data[col_names] , factor)

# Structure learning

source('RAG/LSACData_exp/makeblacklist.R')

# Choose maximum number of parents
n <- ncol(data)
maxparents<-n-1 # Maximum number of parents to allow

# Fill up a matrix with possible parents

parenttable<-listpossibleparents(maxparents,c(1:n))
tablelength<-nrow(parenttable[[1]]) # size of the table

# Now need to score them!

#scoretable<-scorepossibleparents(parenttable,n) 

iterations<-100000 #number of iterations in the chain
moveprobs<-1 # having length 1 disallows the new edge reversal move
#if(!(length(moveprobs)==1)){print('Vector of move probabilities has the wrong length!')}

stepsave<-1 #stepsave<-iterations/1000 #how often to save the result

set.seed(97)

startDAG<-matrix(numeric(n*n),nrow=n) # starting DAG is empty say

revallowed<-1 # allow standard edge reversals

example<-structureMCMC(n,data,startDAG,iterations,stepsave,maxparents,parenttable,scoretable=NULL,revallowed,moveprobs,blacklist = blklist_clg,scoretype = "bic-cg")


# Plot DAG

source('RAG/LSACData_exp/adjMat2ggdag.R') # convert adjacency matrix to ggdag
source('RAG/LSACData_exp/extract.bn.fit.R') # extract coefficients and corresponding std's of these coefficients
source('RAG/LSACData_exp/blackpartition.R') # convert blacklist dataframe to a matrix
load('RAG/LSACData_exp/nodes.RData') # node labels


# sorting the dags output of partition MCMC.
top_down <- order((unlist(example$DAGlogscore)),decreasing = TRUE) 
dags <- example$incidence[top_down]
adjMat <- unique(dags)[[1]]


coef.dat <- extract.bn.fit(as.matrix(adjMat), dat.dag, discvar=TRUE)

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

ggsave("Clg-Plot.pdf", plots, device="pdf",width = 7.5, height = 7.5)

