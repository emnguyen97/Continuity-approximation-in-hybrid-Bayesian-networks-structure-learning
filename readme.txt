This repository contains complementary experimental implementations to the "Continuity Approximation in Hybrid Bayesian Networks Structure Learning" paper.

Three main experiments, as described in the paper, include:
(1) Network with 2 nodes (2nodes_exp), including simulation for Lemma 1;
	- The naming for the files refers to the 4 scenarios described in the paper: Continuous, Discrete, Combination with continuous parent node and discrete child node, Combination with discrete parent node and continuous child node (referred to as CLG). Codes for data generation, implementing models, and tabulating results are included in the corresponding files for each scenario. Fig2.R is the implementation for Figure 2 in the paper.

(2) Network with 10 nodes (10nodes_exp);
	- The main implementation file is 10nodestest.R; Other files are separate codes for data generation, implementing models, and tabulating results.

(3) Experiment on real-world data (with data included) (LSACData_exp).
	The main implementation files are Experiments.R to run partition MCMC on the data and plotDAG.R for Figure 4; Complementary_exp.R is for the complementary study (Table 6). Other files are utility  functions for plotting (converting adjacency matrix to ggdag, extracting coefficients, blacklists).

(4) Additional experiment with a 4-node network (4nodes_exp);
	- The main implementation file is 4nodesExp.r; Other files are separate codes for data generation, implementing models, and tabulating results.


Partition and structure MCMC were used for the task of structure learning, including TABU and Hill-Climbing for comparison for a 10-node network experiment. The BGe score was used for experiments involving continuous data, and the BDe score was used for experiments involving categorical data.
	


References:
Kuipers, J. and Moffa, G., 2017. Partition MCMC for inference on acyclic digraphs. Journal of the American Statistical Association, 112(517), pp.282-299.
https://github.com/annlia/partitionMCMC/tree/master

