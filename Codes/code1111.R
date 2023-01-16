# CODE1111
# TRAINING AND TESTING OF GRAIN YIELD -
# WHEAT INBRED LINES FROM BGLR PACKAGE
rm(list=ls()) # CLEAR WORKSPACE
set.seed(37111)
library(BGLR)
data(wheat)
####################################################
### USE BGLR MATRIX X
x <- wheat.X
y<- wheat.Y[,4]
#### A GAUSSIAN KERNEL ################

kgaus <- function(X,h){
  X <- scale(X,center=TRUE,scale=FALSE)
  S=sqrt(sum(apply(FUN=var, X=x,MARGIN=2)))
  X <- X/S
  D <- as.matrix(dist(X))^2
  K <- exp(-h*D)
}
############################################
##### A LINEAR KERNEL #################
klin1 <- function(X){
  X=scale(X,center=TRUE,scale=FALSE)
  S=sqrt(sum(apply(FUN=var, X=X,MARGIN=2)))
  X=X/S
  K=tcrossprod(X)
}
############  CHOOSE GAUSSIAN KERNEL ##############

#h <- 0.5
h <- 1
# <- 3

K <- kgaus(x,h)
#dim(K)
#qr(K)$rank

############## CHOOSE LINEAR KERNEL ###############
#K <- klin1(x)
#dim(K)
#qr(K)$rank
##################################################
# READ NUMBER OF TRAINING / VALIDATING SPLITS nitval
nitval <- 20
# READ REGULARISATION PARAMETER lambda
lambda <- 0.3

result <- matrix(data=NA, nrow=nitval,ncol=5)

ptm <- proc.time()
for (j in 1:nitval) {
  cat(j, "\n",sep="")
  train = sample(1:nrow(x), floor(0.5 * nrow(x)))
  xt <- x[train, ]
  yt <- y[train]
  xv <- x[-train, ]
  yv <- y[-train]
  YHATt <- matrix(nrow = length(yt), ncol = 1)
  YHATv <- matrix(nrow = length(yv), ncol = 1)
  Ktrain <- K[train, train]
  Kval <- K[-train, train]
  Xt <- cbind(1, Ktrain)
  RHSt <- crossprod(Xt, yt)
  LHSt <- crossprod(Xt)
  LHSt[-1, -1] <- LHSt[-1, -1] + Ktrain * lambda
  #  diag(LHSt) <- diag(LHSt) + c(0, rep(1e-8,ncol(Xt)-1))
  solt <- solve(LHSt, RHSt)
  YHATt <- solt[1] + Ktrain %*% solt[-1]
  YHATv <- solt[1] + Kval %*% solt[-1]
  mset <- mean((YHATt - yt) ^ 2)
  msev <- mean((YHATv - yv) ^ 2)
  cort <- cor(YHATt, yt)
  corv <- cor(YHATv, yv)
  result[j, ] <- c(j, mset, msev, cort, corv)
}
proc.time() - ptm
summary(result[,3])
summary(result[,5])
