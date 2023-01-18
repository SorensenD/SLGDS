# CODE1113
# RANDOM FOREST WITH HEART DATA
rm(list=ls()) # CLEAR WORKSPACE
library(sda)
library(tree)
library(glmnet)
library(randomForest)
# SAheart <- read.table(
#  "http://www-stat.stanford.edu/~tibs/ElemStatLearn/datasets/SAheart.data",
#   sep=",",head=T,row.names=1)
# CHANGE FAMILY HISTORY TO FACTOR
# SAheart$famhist <- as.factor(SAheart$famhist)
# CHANGE RESPONSE VARIABLE FROM INTEGER TO FACTOR
# SAheart$chd <- as.factor(ifelse(SAheart$chd == 1,"D","H"))

#########################################################
# HEART DATA CAN ALSO BE ACCESSED USING PACKAGE loon.data
library(loon.data) # MUST INSTALL PACKAGE loon.data
data("SAheart")
SAheart$chd <- factor(as.numeric(SAheart$chd),levels=c(2,1),labels=c("D","H"))

length(which(SAheart$chd=="H")) # NO HEART DISEASE
length(which(SAheart$chd=="D")) # HEART DISEASE
########################################################
# FUNCTION ACCURACY:
accuracy = function(actual, predicted) {
  mean(actual == predicted)
}
sumd <- data.frame()
nrep <- 1
reslt <- rep(NA,nrep)
mtry <- c(1,3,5,9)
mtry <- 3
set.seed(2)
for (m in mtry){
  cat("mtry ",m,"\n",sep="")
  for ( rep in 1:nrep ) {
    cat("Replicate ",rep,"\n",sep="")
    train=sample(1:nrow(SAheart),nrow(SAheart)/2)
    rf <- randomForest(SAheart$chd ~.,data=SAheart,
                       subset=train,mtry=m,importance=TRUE)
    predict <- predict(rf,SAheart[-train,])
    observed <- SAheart$chd[-train]
    t <- table(observed,predict)
    print(t)
    reslt[rep]<-1-accuracy(predicted=predict,actual=observed)
  }
  sumd <- rbind(
    sumd,c(m,min(reslt),mean(reslt),median(reslt),max(reslt),
           var(reslt)))
}
names(sumd) <- c("mtry","min","mean","median","max","var")
sumd[3]
#####################################################################
# CODE1113 (cont)
rf <- randomForest(SAheart$chd ~.,data=SAheart,mtry=3,importance=TRUE)
rf$mtry
rf$ntree
rf$confusion
####################################################################
# CODE1113 (cont)
importance(rf,type=1)
