# CODE0101SINGH
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
# RESAMPLES TRAIN/TEST DATA AND COMPUTES MSE
# FOR EACH RESAMPLE/REPLICATE
t1 <- seq(1,15,1) # CHOOSE THE FIRST 15 COLUMNS OF x
X1 <-X[,t1]
n <- length(t1) 
datarf <- data.frame(cbind(y,X))
nc <- 1 # EXAMPLE WITH 1 REPLICATE ONLY
res <- matrix(data=NA, nrow=nc,ncol=3)
for(i in 1:nc){
  if(i > 1){train <- sample(1:nrow(datarf),nrow(datarf)/2)}
  glm.fit <- glm(y[train] ~ X1[train,] , 
                 family=binomial(link="logit"))
  # CALCULATE PREDICTED LIABILITY FOR THE TRAINING (liabT) AND THE 
  # VALIDATING DATA (liabV)
  liabV <- X1[-train,1:n]%*%glm.fit$coefficients[2:(n+1)]+
    glm.fit$coefficients[1]
  liabT <- X1[train,1:n]%*%glm.fit$coefficients[2:(n+1)]+
    glm.fit$coefficients[1]
  # COMPUTE Pr(Y=1) BASED ON THESE LIABILITIES  
  probT <- exp(liabT)/(1+exp(liabT))
  probV <- exp(liabV)/(1+exp(liabV))
  # COMPUTE PREDICTED VALUES IN TRAINING AND VALIDATING DATA
  # ON THE 0/1 SCALE
  predT <- ifelse(probT > 0.5, "1", "0")
  predV <- ifelse(probV > 0.5, "1", "0")
  # COMPUTE MISCLASSIFICATION ERROR
  predclassT <- mean((as.numeric(predT) - y.train)^2)
  predclassV <- mean((as.numeric(predV) - y.test)^2)
  
  # IF CURIOUS COMPUTE LOG-LIKELIHOOD, DEVIANCE, 
  # AIC USING TRAINING DATA
  #  ll <- sum(y.train*liabT) - sum(log(1+exp(liabT)))
  #  dev <- -2*ll
  #  AIC <- dev + 2*(n+1)
  # ***********************************  
  
  res[i,] <- c(n,predclassT,predclassV)
}
res
