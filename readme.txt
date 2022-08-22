This repository contains complimentary experimental implementations to the paper Structure Learning for Hybrid Bayesian Networks.

Three main experiments as described in the paper include:
(1) Networks with 2 nodes (2nodes_exp);
	- The naming for the files refers to the 4 scenarios described in the paper: Continuous, Discrete, Combination with continuous parent node and discrete child node, Combination with discrete parent node and continuous child node (referred to as CLG). Codes for data generation, implementing models, tabulating results are inclusive in the corresponding files for each scenarios. Fig2.R is the implementation for Figure 2 in the paper.

(2) Networks with 4 nodes (4nodes_exp);
	- Main implementation file is 4nodesExp.r; Other files are seperate code for data generation, implementing models, and tabulating results.

(3) Experiment on a real-world data (with data included) (LSACData_exp).
	- Main implementation files are Experiments.R to run partition MCMC on the data and plotDAG.R for Figure 4; Complementary_exp.R is for the complementary study (Table 6). Other files are utility  functions for plotting (converting adjacency matrix to ggdag, extracting coefficients, blacklists).


For all of the experiments, partition MCMC was used for the task of structure learning. The BGe score was used for experiments involving continuous data and the BDe score was used for experiments involving categorical data.
	


References:
Kuipers, J. and Moffa, G., 2017. Partition MCMC for inference on acyclic digraphs. Journal of the American Statistical Association, 112(517), pp.282-299.

