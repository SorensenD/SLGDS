# CODE1001
rm(list=ls()) # CLEAR WORKSPACE
set.seed(123)

nindiv <- 2000
nmark <- 150
nsamples <- nindiv*nmark
# GENERATE COVARIATE MATRIX FROM BINOMIAL DISTRIBUTION
X<-matrix(nrow=nindiv,ncol=nmark,
          rbinom(n=nsamples,size=2,p=.5))
#########################################################
# CHOOSE VALUE FOR MEAN mu 
mu <- 10
# CHOOSE VALUE FOR ENVIRONMENTAL VARIANCE ves
ves<-20
b<-matrix(data=0.0,nrow=nmark,ncol=1) # b from operational model

et<- rnorm(nindiv,mean=0,sd=sqrt(ves))
b <- rnorm(nmark,mean=0,sd=2)
y <- mu + X %*%b + et
train <- sample(1:nrow(X),floor(0.5*nrow(X)))
Xt <- X[train,]
yt <- y[train]
Xv <- X[-train,]
yv <- y[-train]
Zt <- cbind(1,Xt)
Zv <- cbind(1,Xv)
##################################################
#####   coefficient matrix LHSt, rhs & solution solt
RHSt <- crossprod(Zt,yt)
LHSt <- crossprod(Zt)
solt <- solve(LHSt,RHSt)
e <- yt-Zt%*%solt
##################################################
rep <- 5000 # NUMBER OF DRAWS USING COMPOSITION
ystartrain <- matrix(data=NA,nrow=length(yt),ncol=1)
ystarval <- matrix(data=NA,nrow=length(yt),ncol=1)

resMSE <- matrix(data=NA,nrow=rep,ncol=3)
res <- matrix(data=NA,nrow=rep,ncol=(nmark+2))

scale <- sum(e^2)
Cinv <- solve(LHSt)
ch <- chol(Cinv)
ptm<-proc.time()

for (i in 1:rep){
  #  print(i)
  df <- length(yt)-nmark-2 
  # DRAW RESIDUAL VARIANCE
  varstar <- scale/rchisq(1,df)
  resid <- rnorm(length(solt),0,1)
  # DRAW LOCATION PARAMETERS
  bstar <- solt + t(ch)%*% resid*sqrt(varstar) 
  ystartrain <- Zt%*%bstar + 
    rnorm(length(ystartrain),0,sd=sqrt(varstar))
  # # DRAW VALIDATING DATA  
  ystarval <- Zv%*%bstar + 
    rnorm(length(ystarval),0,sd=sqrt(varstar)) 
  msevalystar <- mean((ystarval-yv)^2) # MSE_3
  msetrainystar <- mean((ystartrain-yt)^2)
  msevalbstar <- mean((Zv%*%bstar-yv)^2) # MSE_2
  resMSE[i,] <- c(msetrainystar,msevalystar,msevalbstar)
  res[i,] <- c(varstar,bstar)
}
proc.time()-ptm
av <- apply(res,2,mean)
bstarhat <- av[2:(nmark+2)]
ystarhatval <- Zv%*%bstarhat
mse_1 <- mean((ystarhatval-yv)^2) # POINT PREDICTOR
mserrors <- apply(resMSE,2,mean)
mse_2 <- mserrors[3]
mse_3 <- mserrors[2]
ci_2 <- quantile(resMSE[,3],c(0.025,0.975))
ci_3 <- quantile(resMSE[,2],c(0.025,0.975))
#ci_2
#setwd("C:/Users/au223137/Dropbox/Rsessions/SummerCourse/OutputR")
#pdf("C:/Users/au223137/Dropbox/Rsessions/SummerCourse16/
#    OutputR/mse3.pdf")
hist(resMSE[,2],main=NA,xlab=NA)
#dev.off()
#setwd("C:/Users/au223137/Dropbox/Rsessions/SummerCourse/OutputR")
#pdf("C:/Users/au223137/Dropbox/Rsessions/
#     SummerCourse16/OutputR/mse2.pdf")
hist(resMSE[,3],main=NA,xlab=NA)
#dev.off()
