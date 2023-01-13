# CODE1002
#########  SIMULATE BINARY RECORDS FROM LOGISTIC MODEL 
### AND FIT LOGISTIC MODEL (2 COVARIATES) WITH METROPOLIS-HASTINGS 
#### CREATE TRAINING AND TESTING/VALIDATING DATA
rm(list=ls()) # CLEAR WORKSPACE
set.seed(77111)
require(graphics)
# THE CODE WILL USE THE PACKAGE MVTNORM; IT IS INSTALLED BELOW
#install.packages("mvtnorm", .libPaths()[1])
library(mvtnorm)
library(MASS)
#CHOOSE LENGTH OF CHAIN rep
rep<-6000
result<-matrix(data=NA,nrow=rep,ncol=6)
nindiv <- 2000
x1 <- rep(0,nindiv)
x2 <- rep(0,nindiv)
y <- rep(0,nindiv)
p<-0.2
x1 <- runif(nindiv,-1,1)
x2 <- runif(nindiv,-1,1)

b_0 <- log(p/(1-p))
b_1 <- 2
b_2 <- 2
z <- b_0 + b_1*x1 + b_2*x2
p1 <- exp(z)/(1+exp(z))
## GENERATE THE BINARY RECORDS
y <- rbinom(length(p1),1,p1)
# CONSTRUCT THE DATA SET 
dat1 <- matrix(c(y,x1,x2),nrow = nindiv, ncol = 3)
colnames(dat1) <- c("Y", "X1", "X2")
d<-data.frame(dat1)
#attach(d)
Y <- d$Y
X1 <- d$X1
X2 <- d$X2
set.seed(771)
train=sample(1:length(X1),length(X1)/2)
test=(-train)
y.test=Y[test]
y.train<-Y[train]
length(y.train)
###############################################################################
# CODE1002 (cont)
# CHOOSE TUNING PARAMETER LAMBDA AND COVARIANCE MATRIX C
# OF THE METROPOLIS-HASTINGS ALGORITHM
lambda<-0.015
c <- diag(c(1.5,0.6,0.6))
# INITIALISE THE MEAN OF THE TRIVARIATE DISTRIBUTION
theta<-c(-5,1.0,1.0)

## FUNCTION TO COMPUTE THE LOG-POSTERIOR
logpost <- function(data,theta)
{
  interm <- theta[1] + theta[2]*data$X1 + theta[3]*data$X2
  with(data=data,sum(Y*( interm )-log(1+exp(interm))))
}
#START M-H LOOP
ptm <- proc.time()
accept<-0
for (i in 1:rep)
{
  #  print(i)
  #SAMPLE PROPOSAL FOR THETA (Ytheta) FROM N(theta,lamdaC)
  Ytheta<- rmvnorm(1,mean=theta,sigma=lambda*c)
  logalfa<-logpost(d[train,],Ytheta) - logpost(d[train,],theta)
  unif<-runif(1)
  if (unif<exp(logalfa))
  {
    theta[1]<-Ytheta[1]
    theta[2]<-Ytheta[2]
    theta[3]<-Ytheta[3]
    interm <- theta[1] + theta[2]*X1[test] + theta[3]*X2[test]
    
    proby1 <- exp(interm)/(1+exp(interm))
    yhat <- rbinom(length(Y[test]),1,proby1)
    yhatBR <- ifelse(proby1 > 0.5, 1, 0)
    accept<-accept+1
  }
  else
  {
    interm <- theta[1] + theta[2]*X1[test] + theta[3]*X2[test]
    proby1 <- exp(interm)/(1+exp(interm))
    yhat <- rbinom(length(Y[test]),1,proby1)
    yhatBR <- ifelse(proby1 > 0.5, 1, 0)
  }
  misclas <- mean((yhat-Y[test])**2)
  misclasBR <- mean((yhatBR-Y[test])**2)
  brier <- mean((Y[test]-proby1)**2)
  vyhat <- var(yhat)
  logscYV <- sum(Y[test]*log(proby1)+(1-Y[test])*log(1-proby1))
  result[i, ]<-c(i,theta[1],theta[2],theta[3],misclas,misclasBR)
}
proc.time()-ptm
acceptratio <- accept/rep
# PRINT ACCEPTANCE RATIO OF THE JOINT UPDATING
acceptratio
# PRINT THE McMC ESTIMATES OF POSTERIOR MEANS
apply(result[501:rep,5:6],2,mean)
# PRINT 95% POSTERIOR INTERVALS FOR THE MISCLASSIFIATION RATES
misclasQ <- quantile(result[500:rep,5],c(0.025,0.975))
misclasQ
misclasBRQ <- quantile(result[500:rep,6],c(0.025,0.975))
misclasBRQ
# THE McMC ESTIMATES OF THE MEAN OF THE POSTERIOR DISTRIBUTION OF 
# REGRESSION PARAMETERS ARE
meanbetas <- apply(result[500:rep,2:4],2,mean)
meanbetas
# AND THE 95% POSTERIOR INTERVALS ARE
quantile(result[500:rep,2],c(0.025,0.975))
quantile(result[500:rep,3],c(0.025,0.975))
quantile(result[500:rep,4],c(0.025,0.975))
##############################################################################
# CODE1002 (cont)
# FIT MODEL BY ML
set.seed(771)
nrepl <- 1000
resLik <- matrix(data=NA,nrow=nrepl,ncol=3)
ptm <- proc.time()
for (i in 1:nrepl){
  train=sample(1:nrow(d),nrow(d)/2)
  f<-glm(Y ~.,data=d[train,],family="binomial")
  #  summary(f)
  #f<-glm(Y ~1,data=d[train,],family="binomial")
  #  summary(f)
  # PRED PROBABILITIES:
  predV <- predict(f,d[-train,],type="response") 
  BrierFreq <- mean((predV-d$Y[-train])^2) 
  # ASSIGN Y TO ITS CLASS ACCORDING TO BAYES RULE  
  yhatBR <- ifelse(predV > 0.5, 1, 0)
  mseBR <- mean((yhatBR-d$Y[-train])^2) # MSE_1
  ### OR SAMPLE Y FROM ITS PREDICTIVE DISTRIBUTION CONDITIONAL ON 
  ### ML ESTIMATES - THIS MAKES IT COMPARABLE TO THE McMC APPROACH
  yhatPpd <- rbinom(length(Y[-train]),1,predV)
  msePpd <- mean((yhatPpd-d$Y[-train])^2) # MSE_2
  resLik[i,] <- c(i,mseBR,msePpd)
}
proc.time()-ptm
apply(resLik[,2:3],2,mean)
quantile(resLik[,2],c(0.025,0.975))
quantile(resLik[,3],c(0.025,0.975))