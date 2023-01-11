# CODE0305
# STOCHASTIC GRADIENT DESCENT WITH A LINEAR MDEL
rm(list=ls()) # Clear the workspace
set.seed(195021)
N <- 1000
x<-seq(from=0,to=5,length=N)
signal<-10 + 0.2*x
error<-rnorm(N)
y<-signal+error
one <- rep(1,N)
X <- cbind(one,x)
LHS <- crossprod(X) # LEFT HAND SIDE
RHS <- crossprod(X,y) # RIGHT HAND SIDE
bh<- solve(LHS,RHS) # SOLUTION TO THE LEAST SQUARES EQUATIONS

################################
nit <- 10
alfa_0 <- 0.0038
c <- matrix(data=NA, nrow=nit+1,ncol=1)
cost <- function(miu,b){sum(y-miu-b*x)^2}

miu <- 5
b <- 1
c[1] <- cost(miu,b)
for(j in 1:nit){
  alfa <- (1-(j/nit))*alfa_0 + ((j/nit)*alfa_0*0.01)
  for(i in 1:length(y)) {
    #  cat(i, "\n",sep="")
    fdmiu <- -2*(y[i] - miu - b * x[i])
    fdbeta <- -2*((y[i] - miu - b * x[i]) * x[i])
    fd <- matrix(c(fdmiu, fdbeta), nrow = 2, ncol = 1)
    sol0 <- matrix(c(miu, b), nrow = 2, ncol = 1)
    sol1 <- sol0 - alfa * fd
    miu <- sol1[1, 1]
    b <- sol1[2, 1]
  }
}
## CHECK ######################
beta <- c(miu,b) 
beta # STOCHASTIC GRADIENT DESCENT SOLUTION
bh # LEAST SQUARES SOLUTION
