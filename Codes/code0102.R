# CODE0102
# READING SINGH ET EL 2002 DATA
rm(list=ls()) # CLEAR WORKSPACE
# Lasso solutions using package glmnet
#install.packages("glmnet", .libPaths()[1])
#install.packages("sda")
library("sda")
library(glmnet)
data(singh2002)
X<-singh2002$x
y<-ifelse(singh2002$y=="cancer",1,0)
n<-nrow(X)
Xlasso<-X
set.seed(3037)
train=sample(1:nrow(X),nrow(X)/2)
#train=sample(1:nrow(X),70)
test=(-train)
y.test=y[test]
y.train<-y[train]
#
# **********  FOR PREDICTION USING LASSO  *****************
repl <- 1 # NUMBER OF REPLICATES (RESAMPLES TRAINING / VALIDATING)
result <- matrix(data=NA, nrow=repl,ncol=4)
set.seed(3037)
for (i in 1:repl){
  if(i > 1){train <- sample(1:nrow(Xlasso),nrow(Xlasso)/2)}
  y.train <- y[train]
  y.test <- y[-train]
  # STEP 1: cross-validation; find best value of lambda
  # alpha=1: LASSO; alpha=0: RIDGE REGRESSION  
  cv.out=cv.glmnet(Xlasso[train,],y[train],alpha=1,
                   family="binomial",type = "class")
  #plot(cv.out)
  bestlam=cv.out$lambda.min
  #bestlam
  
  # Using best lambda, fit model on training data 
  #  to obtain final parameter estimates
  
  # STEP 2
  fm=glmnet(y=y[train],x=Xlasso[train,],alpha=1,lambda=bestlam,
            family="binomial",type.measure= "class")
  nzcf<-coef(fm)
  cf<-which(fm$beta[,1]!=0)
  if (length(cf) == 0){
    out <-c(i,length(cf))
    print(out)
    break
  }
  #length(cf) # NUMBER OF REGRESSION PARAMETERS IN THE FINAL MODEL
  # CONSTRUCT PREDICTIONS FROM OUTPUT OF fm
  #    1. VALIDATING DATA
  predglmnet<-fm$a0+Xlasso[-train,cf]%*%fm$beta[cf]
  probs <- exp(predglmnet)/(1+exp(predglmnet))
  predclass_test <- as.numeric(ifelse(probs > 0.5, "1", "0"))
  #    2. TRAINING DATA
  predglmnet<-fm$a0+Xlasso[train,cf]%*%fm$beta[cf]
  probs <- exp(predglmnet)/(1+exp(predglmnet))
  predclass_train <- as.numeric(ifelse(probs > 0.5, "1", "0"))
  result[i,] <- c(mean((predclass_train-y.train)^2),
                  mean((predclass_test-y.test)^2),bestlam,length(cf))
}
result
#NOTE: for prediction, GLMNET can be implemented more directly, 
# using in STEP2:
###############################################################
#  fm.predclass=predict(cv.out,s=bestlam,newx=Xlasso[test,],
#        family="binomial",type="class") 
#  mean((as.numeric(fm.predclass)-y.test)^2) # VALIDATION ERROR 
#       RATE (BASED ON CLASS LABELS)