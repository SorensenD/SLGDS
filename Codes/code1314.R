# CODE1314
# PREDICTION EXERCISE 5 
rm(list=ls()) # CLEAR WORKSPACE
set.seed(123)
nindiv<-100
nmark <- 500
nt <- nindiv*nmark
# NUMBER QTL
nqtl <- 50
# GENERATE MARKER MATRIX FROM BINOMIAL DISTRIBUTION
X<-matrix(nrow=nindiv,ncol=nmark,rbinom(n=nt,size=2,p=.5))
#########################################################
# CHOOSE VALUE FOR MEAN mu AND GENOMIC VARIANCE vgs
mu <- 10
vgs<-10
# CHOOSE VALUE FOR ENVIRONMENTAL VARIANCE ves
ves<-20
her <- vgs/(vgs+ves)

btrue<-matrix(data=0.0,nrow=nmark,ncol=1) # parameter from 
#            true model
IDq<-sample(1:nmark,nqtl,replace=F) # from the nmark markers, 
#            choose nqtl as QTL 
QTLeff<-sqrt(vgs/nqtl)# calculate the QTL effect so that the 
#            total genetic variance is VA
btrue[IDq]<-QTLeff # the only b's that are not zero are those 
#           associated with QTL.
W <- matrix(data=NA,nrow= nindiv,ncol=nmark)
cm <- colMeans(X)
# CREATE MATRIX OF STANDARDISED MARKER GENOTYPE CODES
for (i in 1:nmark)
{
  W[,i] <-( X[,i]-cm[i]) / sd(X[,i])
}
# more efficiently, could use:
# W <- scale(X)
# GENERATE nindiv PHENOTYPES WITH MEAN 0, VAR=vgs+ves, 
#           HERITABILITY=vgs/(vgs+ves)
e<- rnorm(nindiv,mean=0,sd=sqrt(ves))
y <- mu + W%*%btrue+ e
k <- (ves/vgs)*nmark # ratio of residual to 
#           genomic variance Vb = vgs/nmark
train <- sample(1:nrow(W),floor(0.5*nrow(W)))
Xt <- W[train,]
yt <- y[train]
Xv <- W[-train,]
yv <- y[-train]
Zt <- cbind(1,Xt)
Zv <- cbind(1,Xv)
###############################################################
#####   ridge regression coefficient matrix, rhs & solution solt
RHSt <- crossprod(Zt,yt)
LHSt <- crossprod(Zt)
LHSt[-1,-1] <- LHSt[-1,-1]+diag(k,nrow=nrow(LHSt)-1)
solt <- solve(LHSt,RHSt)
# PREDICTION, CONDITIONAL ON ESTIMATED LOCATION PARAMETERS (solt)
predval <- Zv%*%solt # VALIDATING
predtrain <- Zt%*%solt # TRAINING
###############################################################
# CODE1314 (cont)
# COMPUTE SAMPLING DISTRIBUTION OF MSE, CONDITIONAL ON 
#            (mu_hat,b_hat) AND VARIANCES
rep <- 10000
res1 <- matrix(data=NA, nrow=rep,ncol=1)

meany <- predval
vary <- diag(ves,nrow=length(yv))
ptm <- proc.time()
for (i in 1:rep){
  yrep <- rnorm(length(yv),meany,sqrt(ves))
  mse1 <- mean((yrep-yv)^2)
  ztz <- (1/ves)*sum((yrep-yv)^2)
  res1[i,] <- mse1
}
proc.time()-ptm
meanmsev <- apply(res1,2,mean)
meanmsev
varmsev <- apply(res1,2,var)
varmsev
# hist(res1[,1])
ncp <- sum((yv-meany)^2)/ves
expmse <- (ves/length(yv))*(length(yv)+ncp)
expmse
varmse <- (2* length(yv) + 4*ncp)*(ves/length(yv))^2
varmse
expQF <- ves + (mean((yv-meany)^2))
expQF
###########################################################
# CODE1314 (cont)
# METHOD OF COMPOSITION: 
# (ACCOUNTING FOR UNKNOWN LOCATION PARAMETERS)
# 1. USING TRAINING DATA Yt, SAMPLE THETA* ~ THETA|Yt
# 2. SAMPLE VALIDATING DATA Yv* ~ Yv|THETA* 
# 3. COMPUTE VALIDATION MSEv = MEAN((Yv*-Yv)^2)
# 4. GOTO 1 UNTIL ENOUGH SAMPLES
rep <- 1000
res2 <- matrix(data=NA, nrow=rep,ncol=2)
theta <- solt
Cinv <- solve(LHSt)
ch <- chol(Cinv*ves)
varcov <- Cinv*ves
ptm <- proc.time()
for (i in 1:rep){
  #  print(i)
  theta <- solt + t(ch)%*%rnorm(length(theta),0,1)
  # DRAWS FROM THE VALIDATING DATA:  
  ystarval <- rnorm(length(yv),Zv%*%theta,sqrt(ves))
  # DRAWS FROM THE TRAINING DATA:  
  ystartrain <- rnorm(length(yt),Zt%*%theta,sqrt(ves))
  
  mse2val <- mean((ystarval-yv)^2) # VALIDATION MSE
  mse2train <- mean((ystartrain-yt)^2) # TRAINING MSE
  res2[i,] <- c(mse2val,mse2train)
}
proc.time()-ptm
# hist(res2[,1])
apply(res2,2,mean)
meanmse2val <- mean(res2[,1])
varmse2val <-  var(res2[,1])
interm <- Zv%*%Cinv%*%t(Zv)
expQF <- ves + (ves*sum(diag(interm)))/length(yv) + 
  mean((predval-yv)^2)
expQF
