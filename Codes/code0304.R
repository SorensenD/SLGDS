# CODE0304
# GRADIENT DESCENT WITH A LINEAR NODEL
rm(list=ls()) # Clear the workspace
set.seed(195021)
N<-100
x<-seq(from=0,to=5,length=N)
signal<-10 + 0.2*x
error<-rnorm(N)
y<-signal+error
one <- rep(1,N)
X <- cbind(one,x)
LHS <- crossprod(X) # LEFT HAND SIDE
RHS <- crossprod(X,y) # RIGHT HAND SIDE
bhat<- solve(LHS,RHS) # SOLUTION TO THE LEAST SQUARES EQUATIONS

################################
nit <- 200
alfa <- 0.002
miu <- matrix(data=NA, nrow=nit+1,ncol=1)
b <- matrix(data=NA, nrow=nit+1,ncol=1)
c <- matrix(data=NA, nrow=nit+1,ncol=1)
cost <- function(miu,b){sum(y-miu-b*x)^2}

miu[1] <- 5
b[1] <- 1
c[1] <- cost(miu[1],b[1])
for(i in 1:nit) {
  fdmiu <- -sum(y - miu[i] - b[i] * x)
  fdbeta <- -sum((y - miu[i] - b[i] * x) * x)
  fd <- matrix(c(fdmiu, fdbeta), nrow = 2, ncol = 1)
  sol0 <- matrix(c(miu[i], b[i]), nrow = 2, ncol = 1)
  alfa <- 0.002
  sol1 <- sol0 - alfa * fd
  miu[i + 1] <- sol1[1, 1]
  b[i + 1] <- sol1[2, 1]
}
## CHECK ######################
beta <- c(miu[i],b[i]) 
beta # GRADIENT DESCENT SOLUTION
bhat # LEAST SQUARES SOLUTION
cost(miu[i],b[i]) # COST FUNCTION EVLUATED AFTER nit ITERATIONS