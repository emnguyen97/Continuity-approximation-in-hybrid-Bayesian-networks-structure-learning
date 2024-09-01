# LSAC data
load("LSACData_exp/imputed.dat.all.RData")

library(BiDAG)
library(bnlearn)
library(tibble)
library(ggdag)
library(dagitty)
library(reshape2)
library(igraph)
library(gRain)
library(dplyr)
library(caret)

## Selecting Wave 2 for analysis

B.W2 <- imputed.dat.k.all[[2]]

B.W2 <- B.W2 %>%
  select(SX,BMI,SE,AC,INC,FS,FH,ME1,FE1,BM1,BM2,RP1,DP1,FV,HF,HSD,SL,OD,GW,BWZ)

# Function to plot an adjacency matrix as a graph
plotadj <- function(adjmat) {
  n <- dim(adjmat)[1]
  e <- empty.graph(colnames(adjmat))
  amat(e) <- adjmat
  return(graphviz.plot(e))
}

# Assume 'data' is your input dataset
data <- dat.dag <- B.W2

# Source external R scripts for black partition and blacklist creation
source("LSACData_exp/blackpartition.R")
source("LSACData_exp/makeblacklist.R")

##### Data Processing ##########

# Add row indices to the data for reference
data <- data %>%
  mutate(row_index = row_number())

# Sample test data ensuring each category of 'FH' is represented proportionally
set.seed(2024)
test_indices <- data %>%
  group_by(FH) %>%
  sample_frac(size = 1/3) %>%
  pull(row_index) # Extract indices of sampled rows

# Create test and train datasets
test_set <- data %>%
  filter(row_index %in% test_indices) %>%
  select(-row_index)

train_set <- data %>%
  filter(!row_index %in% test_indices) %>%
  select(-row_index)

# Discretize the training and test datasets
# Function to discretize data using predefined breaks
discretizedata <- function(data, breaks_list) {
  # Discretize the specified variables using the provided breaks
  data_dscrt <- data %>%
    mutate(
      BMI = cut(BMI, 
                breaks = breaks_list$BMI,
                labels = c("Underweight", "HealthyWeight", "Overweight", "Obese"),
                include.lowest = TRUE),
      SE = cut(SE, 
               breaks = breaks_list$SE,
               labels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
               include.lowest = TRUE),
      INC = cut(INC, 
                breaks = breaks_list$INC,
                labels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
                include.lowest = TRUE),
      BM1 = cut(BM1, 
                breaks = c(0, 18.5, 25.0, 30.0, Inf),
                labels = c("Underweight", "HealthyWeight", "Overweight", "Obese"),
                include.lowest = TRUE),
      BM2 = cut(BM2, 
                breaks = c(0, 18.5, 25.0, 30.0, Inf),
                labels = c("Underweight", "HealthyWeight", "Overweight", "Obese"),
                include.lowest = TRUE),
      OD = cut(OD, 
               breaks = breaks_list$OD,
               labels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
               include.lowest = TRUE),
      BWZ = cut(BWZ, 
                breaks = breaks_list$BWZ,
                labels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
                include.lowest = TRUE)
    )
  
  # Convert all numeric columns to factors
  data.f <- data_dscrt %>%
    mutate_if(is.numeric, as.factor) %>%
    as.data.frame()
  
  return(data.f)
}

# Step 1: Calculate the quantile breaks on the training set
train_breaks <- list(
  BMI = quantile(train_set$BMI, c(0, 0.05, 0.85, 0.95, 1)),
  SE = quantile(train_set$SE, c(0, 0.05, 0.25, 0.75, 0.95, 1)),
  INC = quantile(train_set$INC, c(0, 0.05, 0.25, 0.75, 0.95, 1)),
  OD = quantile(train_set$OD, c(0, 0.05, 0.25, 0.75, 0.95, 1)),
  BWZ = quantile(train_set$BWZ, c(0, 0.05, 0.25, 0.75, 0.95, 1))
)

# Step 2: Discretize the training dataset using its own quantile breaks
train_dscrt <- discretizedata(train_set, train_breaks)

# Step 3: Discretize the test dataset using the training set's quantile breaks
test_dscrt <- discretizedata(test_set, train_breaks)


##### RAG ##########

# Obtain the estimated DAG using the BGe score and partition MCMC
# bgeScore <- scoreparameters(scoretype = "bge", data = train_set, bgepar = list(edgepf = floor(2*log(ncol(train_set)))))
bgeScore <- scoreparameters(scoretype = "bge", data = train_set)
partbge <- partitionMCMC(bgeScore, blacklist = blklist, iterations = 100000)

bge_MAP_estimate <- as.matrix(partbge$DAG)

# Plot the adjacency matrix of the estimated DAG
plotadj(bge_MAP_estimate)

# Get the parent nodes for BMI from the estimated DAG
BMI_pa <- names(bge_MAP_estimate[,'BMI'][bge_MAP_estimate[,'BMI'] != 0])

# Perform linear regression on training data using parents of BMI
BMI_dat_train <- train_set[, c("BMI", BMI_pa)]
reg <- lm(BMI ~ ., BMI_dat_train)

# Predict BMI for the test set
predictions <- predict(reg, newdata = test_set) 
predictions <- as.data.frame(predictions)
colnames(predictions) <- "BMI"

# Discretize the predicted BMI values using the training set quantiles
preds_dscrt <- cut(predictions$BMI,
                       breaks = train_breaks$BMI,
                       labels = c("Underweight", "HealthyWeight", "Overweight", "Obese"),
                       include.lowest = TRUE)

# Get the confusion matrix and classification metrics
RAG_confusion <- confusionMatrix(preds_dscrt, test_dscrt$BMI)

# Display the confusion matrix
print(RAG_confusion$table)

# Weighted F1 Score
weighted_f1 <- mean(RAG_confusion$byClass[,'F1'], na.rm = TRUE)
print(paste("RAG Weighted F1 Score:", weighted_f1))

# Overall metrics
print(RAG_confusion$overall)

calculate_weighted_accuracy_f1 <- function(actual, predicted) {
  # Ensure both are factors and have the same levels
  actual <- factor(actual)
  predicted <- factor(predicted, levels = levels(actual))
  
  # Generate the confusion matrix
  cm <- confusionMatrix(predicted, actual)
  
  # Initialize vectors to store accuracies and F1 scores for each class
  class_accuracies <- numeric(length = length(levels(actual)))
  class_f1_scores <- numeric(length = length(levels(actual)))
  
  # Calculate metrics for each class
  for (class in levels(actual)) {
    class_index <- which(levels(actual) == class)
    # True Positives
    TP = cm$table[class_index, class_index]
    # False Positives
    FP = sum(cm$table[class_index, ]) - TP
    # False Negatives
    FN = sum(cm$table[, class_index]) - TP
    # Calculate Accuracy for this class
    class_accuracies[class_index] = TP / (TP + FP + FN)
    
    # Calculate Precision and Recall for this class
    precision = ifelse(TP + FP > 0, TP / (TP + FP), 0)
    recall = ifelse(TP + FN > 0, TP / (TP + FN), 0)
    
    # Calculate F1 Score for this class
    class_f1_scores[class_index] = ifelse(precision + recall > 0, 
                                          2 * precision * recall / (precision + recall), 
                                          0)
  }
  
  # Calculate weights based on the actual distribution
  weights <- prop.table(table(actual))
  
  # Calculate Weighted Accuracy and F1
  weighted_accuracy <- sum(class_accuracies * weights, na.rm = TRUE)
  weighted_f1 <- sum(class_f1_scores * weights, na.rm = TRUE)
  
  # Return the weighted metrics
  return(list(
    Weighted_Accuracy = weighted_accuracy,
    Weighted_F1 = weighted_f1
  ))
}


calculate_weighted_accuracy_f1(test_dscrt$BMI,preds_dscrt)


##### DISC ##########

# Calculate BDe (Bayesian Dirichlet equivalent) score for categorical data in the training set
bdecatScore <- scoreparameters("bdecat", train_dscrt)

# Perform MCMC sampling to estimate the most probable DAG (Directed Acyclic Graph)
partcat <- partitionMCMC(bdecatScore, blacklist = blklist, iterations = 100000)

# Extract the DAG from the MCMC output
bde_MAP_estimate <- as.matrix(partcat$DAG)

# Plot the adjacency matrix of the estimated DAG
plotadj(bde_MAP_estimate)

# Identify the parent nodes of BMI in the estimated DAG
BMI_pa_dcrt <- names(bde_MAP_estimate[,'BMI'][bde_MAP_estimate[,'BMI'] != 0])

# Select the variables (BMI and its parents)
var <- c("BMI", BMI_pa_dcrt)

# Create an empty graph structure with only the relevant variables (BMI and its parents)
bdecatDAG <- empty.graph(var)

# Fill in the adjacency matrix with the structure from the estimated DAG
amat(bdecatDAG) <- bde_MAP_estimate[var, var]

# Fit the Bayesian network using the test dataset
bnfit <- bn.fit(bdecatDAG, test_dscrt[, var], method = "bayes")

# Initialize a vector to store the predictions
set.seed(2024)
preds <- vector("character", nrow(test_dscrt))

# Generate predictions for BMI based on the fitted Bayesian network
for (id in 1:nrow(test_dscrt)) {
  preds[id] <- sample(levels(test_dscrt$BMI), size = 1, prob = bnfit$BMI$prob[, test_dscrt$BM1[id]])
}

# Convert predictions to a factor with the same levels as the training set
preds_factor <- factor(preds, levels = levels(train_dscrt$BMI))

# Compute and display the confusion matrix to assess the classification performance
DISC_confusion_matrix <- caret::confusionMatrix(preds_factor, test_dscrt$BMI)

# Display the confusion matrix
print(DISC_confusion_matrix$table)

# Weighted F1 Score
print(paste("DISC Weighted F1 Score:", mean(DISC_confusion_matrix$byClass[,'F1'], na.rm = TRUE)))

# Overall metrics
print(DISC_confusion_matrix$overall)

calculate_weighted_accuracy_f1(test_dscrt$BMI,preds_factor)

#####CLG###################################

source('structureMCMC/combinations.R')
source('structureMCMC/scoretables.R')
source('structureMCMC/structurefns.R')
source('structureMCMC/samplefns.R')
source('structureMCMC/param_utils.R')
source('structureMCMC/structureMCMC.R')

runStrcmcmc <- function(blacklist, data){
  
  # Choose maximum number of parents
  n <- ncol(data)
  maxparents<-n-1 # Maximum number of parents to allow
  
  # Fill up a matrix with possible parents
  
  parenttable<-listpossibleparents(maxparents,c(1:n))
  
  iterations<-100000 #number of iterations in the chain
  moveprobs<-1 # having length 1 disallows the new edge reversal move
  #if(!(length(moveprobs)==1)){print('Vector of move probabilities has the wrong length!')}
  
  stepsave<-1 #stepsave<-iterations/1000 #how often to save the result
  
  set.seed(97)
  
  startDAG<-matrix(numeric(n*n),nrow=n) # starting DAG is empty say
  
  revallowed<-1 # allow standard edge reversals
  
  example<-structureMCMC(n,data,startDAG,iterations,stepsave,maxparents,parenttable,scoretable=NULL,revallowed,moveprobs,blacklist = blacklist,scoretype = "bic-cg")
  
  DAGscores<-unlist(example$DAGlogscore)
  maxDAG<-example$incidence[which(DAGscores==max(DAGscores))][[1]]
  colnames(maxDAG) <- rownames(maxDAG) <- colnames(data)
  
  return(list(DAGlogscores=DAGscores, DAG=maxDAG))
}

source('LSACData_exp/makeblacklist.R')

cont_nodes <- c("BMI","INC","SE","BM1","BM2","BWZ", "OD")
cat_nodes <- setdiff(colnames(train_set),cont_nodes)

train_clg <- train_set %>% mutate_if(names(train_set) %in% cat_nodes, as.factor)
test_clg <- test_set %>% mutate_if(names(test_set) %in% cat_nodes, as.factor)

result <- runStrcmcmc(blacklist=blklist_clg, train_clg)

# Extract the DAG from the MCMC output
clg_MAP_estimate <- result$DAG

# Get the parent nodes for BMI from the estimated DAG
BMI_pa_clg <- names(clg_MAP_estimate[,'BMI'][clg_MAP_estimate[,'BMI'] != 0])

# Perform linear regression on training data using parents of BMI
BMI_dat_train <- train_set[, c("BMI", BMI_pa_clg)]
reg <- lm(BMI ~ ., BMI_dat_train)

# Predict BMI for the test set
predictions_clg <- predict(reg, newdata = test_clg) 
predictions_clg <- as.data.frame(predictions_clg)
colnames(predictions_clg) <- "BMI"

# Discretize the predicted BMI values using the training set quantiles
preds_clg_dscrt <- cut(predictions_clg$BMI,
                       breaks = train_breaks$BMI,
                       labels = c("Underweight", "HealthyWeight", "Overweight", "Obese"),
                       include.lowest = TRUE)

# Get the confusion matrix and classification metrics
CLG_confusion <- confusionMatrix(preds_clg_dscrt, test_dscrt$BMI)

# Display the confusion matrix
print(CLG_confusion$table)

# Weighted F1 Score
print(paste("CLG Weighted F1 Score:", mean(CLG_confusion$byClass[,'F1'], na.rm = TRUE)))

# Overall metrics
print(CLG_confusion$overall)



