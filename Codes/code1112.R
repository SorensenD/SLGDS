# CODE1112
# READ SOUTH AFRICAN HEART DISEASE DATA
rm(list=ls()) # CLEAR WORKSPACE
library(sda)
library(tree)
library(glmnet)
library(randomForest)
# READ HEART DISEASE DATA FROM:
####################################################################
# 1. 
#  SAheart <- read.table(
#  "http://www-stat.stanford.edu/~tibs/ElemStatLearn/datasets/SAheart.data",  
#  sep=",",head=T,row.names=1)
##########  OR FROM  ###############################################
#  READ HEART DISEASE DATA FROM ESL
# 2. 
# SAheart <- 
#  read.table("https://hastie.su.domains/ElemStatLearn/datasets/SAheart.data",
#             sep=",",head=T)

#SAheart <- SAheart[,-1] # IF READ FROM HERE, DROP FIRST COLUMN
###  
###################################################################
 # IF DATA DOWNLOADED FROM 1 OR 2 
#  MUST CHANGE FAMILY HISTORY AND CHD TO FACTOR
# SAheart$famhist <- as.factor(SAheart$famhist)
 # CHANGE RESPONSE VARIABLE FROM INTEGER TO FACTOR
# SAheart$chd <- as.factor(ifelse(SAheart$chd == 1,"D","H"))
#######   OR FROM  ##################################
# HEART DATA CAN ALSO BE ACCESSED USING PACKAGE loon.data
# NO NEED TO CHANGE FAM HISTORY OR CHD TO FACTOR
# 3. 
library(loon.data) # MUST INSTALL PACKAGE loon.data
data("SAheart")
SAheart$chd <- factor(as.numeric(SAheart$chd),levels=c(2,1),labels=c("D","H"))
length(which(SAheart$chd=="H")) # NO HEART DISEASE
length(which(SAheart$chd=="D")) # HEART DISEASE

 ########################################################
# FUNCTION "ACCURACY":
accuracy = function(actual, predicted) {
  mean(actual == predicted)
}
# FIT TREE TO TRAINING AND VALIDATING HEART DATA
replicate <- 1
result <- matrix(data=NA, nrow=replicate,ncol=2)
# REPLICATION
for(i in 1:replicate){
  # SPLIT DATA INTO TRAINING / VALIDATING SET
  set.seed(31)
  train=sample(1:nrow(SAheart),nrow(SAheart)/2)
  treetrain <- SAheart[train,]
  treevalid <- SAheart[-train,]
  # FIT TREE TO THE TRAINING DATA
  treetr <- tree(SAheart$chd[train] ~. ,data=treetrain)
  trtrpred <- predict(treetr,treetrain,type="class")
  trvapred <- predict(treetr,treevalid,type="class")
  table(predicted = trvapred, actual = SAheart$chd[-train] )
  table(predicted = trtrpred, actual = SAheart$chd[train] )
  acval <- accuracy(actual=SAheart$chd[-train],predicted=trvapred)
  actst <- accuracy(actual=SAheart$chd[train],predicted=trtrpred)
  result[i,] <- c(1-acval,1-actst)
}
summary(treetr)
result
############################################################################
# CODE1112 (cont)
# PRUNE TREE FROM TRAINING RUN 
set.seed(33)
treetr_cv <- cv.tree(treetr, FUN = prune.misclass)
#treetr_cv
# index of tree with minimum error
min_idx = which.min(treetr_cv$dev)
#min_idx
# number of terminal nodes in that tree
treetr_cv$size[min_idx]
# misclassification rate of each tree
#treetr_cv$dev / length(train)
# IT APPEARS THAT A TREE WITH treetr_cv$size[min_idx] 
# TERMINAL NODES HAS THE SMALLER MISCLASSIFCATION
# EXECUTE PRUNE.MISCLASS: TREE WITH LOWEST C-V ERROR RATE
# SET best = TO NR TERMINAL NODES FROM THIS BEST TREE
tree_prune <- prune.misclass(treetr,best=treetr_cv$size[min_idx])
#summary(tree_prune)
# OBTAIN PREDICTIONS USING THIS PRUNED TREE
tree_prune_trn <- predict(tree_prune, SAheart[train,],type="class")
table(predict = tree_prune_trn, actual = SAheart$chd[train] )
1-accuracy (predict = tree_prune_trn, actual = SAheart$chd[train])
tree_prune_val <- predict(tree_prune, SAheart[-train,],type="class")
table(predict = tree_prune_val, actual = SAheart$chd[-train] )
1-accuracy (predict = tree_prune_val, actual = SAheart$chd[-train])
tree_prune
############################################################################
# CODE1112 (CONT)
# A LITTLE CALCULATION: PROPORTION OF H's IN LEFT BRANCH
# OF THE FIRST SPLIT OF THE PRUNED TREE
trtrlt <- treetrain[which(treetrain$age < 51.5),]
dim(trtrlt)
trtrltH <- trtrlt[which(trtrlt$chd == "H"),]
dim(trtrltH)[1]/dim(trtrlt)[1]
summary(tree_prune)

