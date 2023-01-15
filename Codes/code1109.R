# CODE1109
# SIMULATING BINARY DATA TO BE ANALYSED WITH KLR
#  USES THE X MATRIX FROM WHEAT IN BGLR
rm(list=ls()) # CLEAR WORKSPACE
library(BGLR)
data(wheat)
X <- wheat.X
set.seed(371)
####################################################
### USE BGLR MATRIX X
nindiv <-nrow(X)
nmark <- ncol(X)
###################################################
nloci<-20
p<-0.5
mu<-log(p/(1-p))

##### GENERATE LIABILITY ###################
va<-1.0 # additive variance of liability
Xc<-matrix(data=NA,nrow=nindiv,ncol= nmark)
be<-matrix(data=0.0,nrow=nmark,ncol=1) # parameter from true model
y<-rep(0,nindiv)
cm<-colMeans(X)
### CENTER AND SCALE X #################
for (i in 1:nmark)
{Xc[,i]<- (X[,i]-cm[i])/sd(X[,i])
}
QTLeff<-sqrt(va/nloci)# calculate the QTL effect so that the 
# total genetic variance is VA
IDq<-sample(1:nmark,nloci,replace=F) # from the nmark markers, 
# choose nloci as QTL
be[IDq]<-QTLeff # the only b's that are not zero are those 
# associated with QTL.
########### GENERATE PHENOTYPIC BINARY DATA ###################
xb<-Xc%*%be
pr <- exp(mu+xb)/(1+exp(mu+xb))
y <- rbinom(nindiv,1,pr)
#sum(y)/length(y) # OBSERVED PROPORTION OF 1'S IN SAMPLE
mean(y)
nitnr <- 10 # NUMBER OF N-R ITERATIONS
nrep <- 10 # NUMBER OF TRAINING / TESTING REPLICATES

#lambda <- 0.0 # ZERO PENALTY !!!!!!!!
lambda <- 0.4

newcostvnr <- rep(0,nrep)
res <- matrix(data=NA, nrow=nrep,ncol=9)
resulttvnr <- matrix(data=NA, nrow=nitnr,ncol=8)

msev <- rep(NA,nrep)
msevnr <- rep(NA,nrep)

#### A GAUSSIAN KERNEL ################
kgaus <- function(X,h){
  X <- scale(X,center=TRUE,scale=FALSE)
  S=sqrt(sum(apply(FUN=var, X=X,MARGIN=2)))
  X <- X/S
  D <- as.matrix(dist(X))^2
  K <- exp(-h*D)
}
############  CHOOSE KERNEL ##############
h <- 0.5
K <- kgaus(Xc,h)
#dim(K)
#qr(K)$rank
########################################
prob1 <- function(miu,alfa,K){
  pr <- exp(miu+K%*%alfa)/(1+exp(miu+K%*%alfa))
}
cost <- function(miu,alfa,K,y){-sum(y*(miu+K%*%alfa) - 
                                      log(1 + exp(miu+K%*%alfa))) + lambda*crossprod(alfa)}

######### ###########    NEWTON-RAPHSON   ####################
### FIT MODEL TO TRAINING DATA AND TEST IN VALIDATING DATA ###
set.seed(77131111)
ptm <- proc.time()
for (i in 1:nrep) {
  #    cat(i, "\n",sep="")
  train=sample(1:nrow(K),floor(0.5*nrow(K)))
  Ktrain <- K[train,train]
  Kval <- K[-train,train]
  ytrain <- y[train]
  yval <- y[-train]
  
  delta <- diag(c(0,rep(lambda,ncol(Ktrain))))
  M <- cbind(0,Ktrain)
  M <- rbind(0,M)
  ###### START VALUES FOR MIU AND ALFA ################
  miu <- 0.0
  alfa <- rep(0.0, ncol(Ktrain))
  W <- matrix(data = 0,nrow = ncol(Ktrain),ncol = ncol(Ktrain))
  for (j in 1:nitnr) {
    fdmiu <- -sum(ytrain - prob1(miu, alfa, Ktrain))
    fdalfa <-
      -Ktrain %*% (ytrain - prob1(miu, alfa, Ktrain)) + 
      lambda * Ktrain %*% alfa
    fd <- matrix(c(fdmiu, fdalfa), nrow = length(alfa) + 
                   1, ncol = 1)
    W <-
      diag(c(prob1(miu, alfa, Ktrain) * 
               (1 - prob1(miu, alfa, Ktrain))))
    Z <- cbind(1, Ktrain)
    zwz <- t(Z) %*% W %*% Z
    LHS <- zwz + lambda * M
    RHS <- -fd
    sol0 <- matrix(c(miu, alfa), nrow = 
                     length(alfa) + 1, ncol = 1)
    sol1 <- sol0 + solve(LHS,RHS)
    miu <- sol1[1, 1]
    alfa <- sol1[-1, 1]
    newcostvnr[j] <- cost(miu, alfa, Ktrain, ytrain)
    resulttvnr[j, ] <- c(j, newcostvnr[j], miu, alfa[1:5])
  }
  probval <- prob1(miu, alfa, Kval)
  y_predval <- as.numeric(ifelse(probval > 0.5, 1, 0))
  msev[i] <- mean((y_predval - yval) ^ 2)
  res[i,] <- c(i,j,newcostvnr[j],miu,alfa[1:5])
}
proc.time()-ptm
tail(resulttvnr[,1:6])
#tail(res)
summary(msev)
