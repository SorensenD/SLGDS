# CODE1113
# RANDOM FOREST WITH HEART DATA
rm(list=ls()) # CLEAR WORKSPACE
library(sda)
library(tree)
library(glmnet)
library(randomForest)
# sahd <- read.table(
#  "http://www-stat.stanford.edu/~tibs/ElemStatLearn/datasets/SAheart.data",
#   sep=",",head=T,row.names=1)
#########################################################
# HEART DATA CAN ALSO BE ACCESSED USING PACKAGE loon.data
library(loon.data) # MUST INSTALL PACKAGE loon.data
data("SAheart")
length(which(SAheart$chd=="No")) # NO HEART DISEASE
length(which(SAheart$chd=="Yes")) # HEART DISEASE
########################################################
# FUNCTION ACCURACY:
accuracy = function(actual, predicted) {
  mean(actual == predicted)
}
# READ THE HEART DATA FROM DISK
#sahd <- read.table(
#  "C:/Users/au223137/Dropbox/Rsessions/HeartDiseaseData/SAheartData",
#  header=TRUE)
#head(sahd)
#str(sahd)
#head(sahd)
#str(sahd)
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
