# CODE1108
# BAYESIAN KERNELISED REGRESSION
# WHEAT DATA FROM BGLR
# CODE ASSUMES K IS OF FULL RANK
# THEREFORE IT DOES NOT WORK WITH CENTERED G. 
# IT WORKS WITH GAUSSIAN KERNEL


rm(list=ls()) # CLEAR WORKSPACE
set.seed(37111)
library(BGLR)
data(wheat)
### USE BGLR MATRIX X
X <- wheat.X
y<- wheat.Y[,1]
nindiv<-length(y)
nmark<-ncol(X)

#### A GAUSSIAN KERNEL ################

kgaus <- function(X,h){
  X <- scale(X,center=TRUE,scale=FALSE)
  S=sqrt(sum(apply(FUN=var, X=X,MARGIN=2)))
  X <- X/S
  D <- as.matrix(dist(X))^2
  K <- exp(-h*D)
}
############  CHOOSE GAUSSIAN KERNEL ##############
#h <- 0.5
h <- 1
#h <- 3.0

K <- kgaus(X,h)
dim(K)
qr(K)$rank
G<-K
###################################################
# EIGEN DECOMPOSITION OF G
EVD <- eigen(G)
names(EVD)
head(EVD$values)
U <- EVD$vector
tU<-t(U)
val <- EVD$values
summary(val)
D <- diag(val,nrow=nindiv)
#Dp IS A VECTOR WITH NON-ZERO EIGENVALUES
Dp<-c(val[1:nindiv])
#INITIALISE Ve
Ve<-0.5
#INITIALISE Vg
Vg<-0.5
#INITIALISE k
k<-Ve/Vg
#INITIALISE VECTOR ALFA
alfa<-rep(0,nindiv)
# CHOOSE LENGTH OF GIBBS CHAIN
rep<-10000
#INITIALISE result
result<-matrix(data=NA,nrow=rep,ncol=5)
# START GIBBS CHAIN
ptm <- proc.time()
for (i in 1:rep)
{
  cat(i, "\n",sep="")
  # SAMPLE mu
  avmu<-sum(y-U%*%alfa)/nindiv
  varmu<-Ve/nindiv
  mu<-rnorm(1,mean=avmu,sd=sqrt(varmu))
  # SAMPLE alfa1 (VECTOR OF LENGTH nindiv)
  meanalfa1<-(Dp/(Dp+k))*tU%*%(y-mu)
  varalfa1<-((Dp)/(Dp+k))*Ve
  alfa1<-rnorm(nindiv,meanalfa1,sqrt(varalfa1))
  alfa<-alfa1
  # SAMPLE Vg
  # COMPUTE SCALE
  scVg<-sum(alfa1*alfa1*(1/Dp))
  Vg<-scVg/rchisq(1,nindiv-2)
  #Vg<-0.0001
  # SAMPLE Ve
  # COMPUTE SCALE
  ystar<-y-mu-U%*%alfa
  scVe<-sum(ystar*ystar)
  Ve<-scVe/rchisq(1,nindiv-2)
  k<-Ve/Vg
  ualfa <- U%*%alfa
  result[i,]<-c(i,mu,Vg,Ve,k)
  #  print(result[i,])
}
proc.time()-ptm
apply(result[1000:rep,2:5],2,mean)
# MEAN AND 95% POSTERIOR INTERVAL FOR LAMBDA
lambdamean <- mean(result[3000:rep,],5)
quantile(result[3000:rep,5],c(0.025,0.975))
###############################################
#CODE FOR THE MC VARIANCE BASED ON BATCHING
# result[,5] CORRESPONDS TO LAMBDA
y <- result[,5] # READS IN ALL DRAWS STORED IN RESULT
#choose number of batches b
b<-500
batch_size <- length(y)/b
batch_size
x<-matrix(y,ncol=b, byrow=FALSE)
avrb<-apply(x,2,mean)
mcvarb<-var(avrb)/length(avrb)
sqrt(mcvarb)
efchsizebatch<-var(y)/mcvarb
efchsizebatch
#############################################
### PLOT AUTOCORRELATION VERSUS LAG USING
### R-FUNCTION acf
require(graphics)
acf(y) ## AUTOCORRELATION OF McMC DRAWS
acf(avrb) ## AUTOCORRELATION OF BATCH MEANS
```