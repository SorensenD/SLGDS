# CODE1315
#FULL-SIB CONTINUOUS DATA
rm(list=ls()) # Clear the workspace
set.seed(123771)
ptm<-proc.time()
require(graphics)
# INITIALISE PARAMETERS
mus<-10 # MEAN
vfs<-1 #VARIANCE BETWEEN FULL-SIBS
#vfs<-0.5 #VARIANCE BETWEEN FULL-SIBS
#vfs <- 0.1
# RESIDUAL VARIANCE
ves<-5
k <- ves/vfs
nf<-500 # NUMBER OF FULL-SIB FAMILIES
n<-2 # FULL-SIB FAMILY SIZE
nb <- 2 # NUMBER OF BREEDS
N<-nf*n # TOTAL NUMBER OF RECORDS
y<-matrix(data=0,nrow=nf*n,ncol=1)
z<-matrix(data=0,nrow=nf*n,ncol=1)
# GENERATE nf FULL-SIB EFFECTS fs
fs<-rnorm(nf,mean=0,sd=sqrt(vfs))
# BREED EFFECTS
br <- rep(0,nb)
br[1] <- 5
br[2] <- 8
# GENERATE nf*n RESIDUAL EFFECTS 
es<-rnorm(nf*n,mean=0,sd=sqrt(ves))
################################################
## GENERATING A FULL-SIB STRUCTURE
b <- rep(1:nb,each=N/2)
z <- rep(1:nf,each=n)
y <- br[b] + fs[z] + es
d <- data.frame(y,z)
################################################
d<-data.frame(y,z)
# GENERATE INCIDENCE MATRICES X & Z
family <- z
breed <- b
family <- as.factor(family)
breed <- as.factor(breed)
X<-model.matrix(~0+breed)
Z<-model.matrix(~0+family)
W <- cbind(X,Z)
LHS <- crossprod(W) # LHS OF MME
LHS[-(1:2),-(1:2)] <- LHS[-(1:2),-(1:2)]+diag(k,nrow=nrow(LHS)-2)
RHS <- crossprod(W,y) # RHS OF MME
SOL <- solve(LHS,RHS) # SOLUTION
HAT <- W%*%solve(LHS)%*%t(W)
V <- Z%*%t(Z)*vfs + diag(ves,nrow=length(y))
COVyyhat <- sum(diag(V%*%t(HAT)))
lambda <- 1/k
Vtilde <- (Z%*%t(Z)*lambda + diag(1,nrow=length(y)))
df <- sum(diag(Vtilde%*%HAT))
yhat <- HAT%*%y
MSEt <- mean((y-yhat)^2)
MSEt
optim1 <- 2*COVyyhat/length(y)
optim1
MSEv <- MSEt + optim1
MSEv
#####################################################
# MONTE CARLO ESTIMATE OF OPTIMISM
# CODE1315 (cont)
## SIMULATE DATA & STORE Y, YHAT
rep <- 1000
gemY <- matrix(data=NA,nrow=rep,ncol=length(y))
gemYHAT <- matrix(data=NA,nrow=rep,ncol=length(y))

br[1] <- SOL[1]
br[2] <- SOL[2]
for (i in 1:rep){
  fs<-rnorm(nf,mean=0,sd=sqrt(vfs))
  es<-rnorm(nf*n,mean=0,sd=sqrt(ves))
  y <- br[b] + fs[z] + es
  RHS <- crossprod(W,y)
  SOL <- solve(LHS,RHS)
  yhat <- W%*%SOL
  gemY[i,] <- y
  gemYHAT[i,] <- yhat
}
sumcov <- 0
## COMPUTE SUM(COV(Y,YHAT))
for (i in 1:length(y)){
  sumcov <- sumcov + cov(gemY[,i],gemYHAT[,i])
}
#####################################
##  A MORE EFFICIENT AND LESS TRANSPARENT CODE IS
# sumcov <- 
# sum(sapply(1:length(y),FUN=function(i)
# {cov(gemY[,i],gemYHAT[,i])}))
#####################################
sumcov
optim2 <- 2*sumcov/length(y)
optim2
