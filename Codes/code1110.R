# CODE1110
# FIT THE NEURAL NETWORK TO SIMULATED BINARY PHENOTYPES
# USING THE WHEAT INBRED LINES WITH THE 1,279 GENETIC MARKERS

rm(list=ls()) # CLEAR WORKSPACE
set.seed(37111)

library(BGLR)
data(wheat)
####################################################
### USE BGLR MATRIX X

X <- wheat.X
nindiv <-nrow(X)
nmark <- ncol(X)
###################################################
# NUMBER OF LOCI AFFECTING THE SIMULATED DATA
nloci<-20
p<-0.5
mu<-log(p/(1-p)) # BIAS TERM FOR SIMULATED DATA

##### INITIALISE PARAMETERS AND ALLOCATE MATRICES ##############
va<-1.0 # additive variance of liability
ve<-1.0 # environmental variance
Xc<-matrix(data=NA,nrow=nindiv,ncol= nmark)
be<-matrix(data=0.0,nrow=nmark,ncol=1) # parameter of true model
y<-rep(0,nindiv)
cm<-colMeans(X)
for (i in 1:nmark)
{Xc[,i]<- (X[,i]-cm[i])/sd(X[,i])
}
QTLeff<-sqrt(va/nloci)# QTL effect so that the total
# genetic variance is VA
IDq<-sample(1:nmark,nloci,replace=F) # from the nmark markers,
# choose nloci as QTL
be[IDq]<-QTLeff # the only b's that are not zero are those
# associated with QTL.
########### GENERATE BINARY DATA y ###################
xb<-Xc%*%be
pr <- exp(mu+xb)/(1+exp(mu+xb))
y <- rbinom(nindiv,1,pr)
######################################################
df = data.frame(cbind(X,y))
m <- length(y)
# BIG X!
X <- t(X)
p <- nrow(X)

Y <- matrix(y,nrow=1,ncol=length(y))

n_1 <- p # NUMBER OF FEATURES IN INPUT DATA
# READ NUMBER OF NEURONS IN LAYER 2
n_2 <- 5
# READ NUMBER OF NEURONS IN LAYER 3
n_3 <- 1
#########################################
### FUNCTIONS:
# SIGMOID FUNCTION
sigm <- function(par){
  1/(1+exp(-par))
}
# COST FUNCTION EXCLUDING REGULARISATION TERM
cost <- function(A,Y){-(tcrossprod(Y,log(A))+
                          tcrossprod((1-Y),(log(1-A))))/m}
#############################################

# READ REGULARISATION PARAMETER delta

delta <- 0.02

# READ GD LEARNING RATE (HERE LABELLED gamma)
gamma <- 0.08

# eps: range of initial values of elements of W_1: (-eps,eps)
eps <- 0.85
# READ NUMBER OF GRADIENT DESCENT ITERATIONS
nit <- 5000
# READ NUMBER OF TRAINING / VALIDATING REPS
nitval <- 20

resultval <- matrix(data=NA, nrow=nitval,ncol=3)
result <- matrix(data=NA, nrow=nit,ncol=8+n_2)

ptm<-proc.time()

for (j in 1:nitval) {
  # INITIALISE MATRIX OF WEIGHT W_1 (n_2 x n_1), n_1=p=rows of X
  # INITIALISE MATRIX OF WEIGHT W_2 (n_3 x n_2)
  W_1 <- matrix(nrow=n_2,ncol=n_1,runif(n_2*n_1,-eps,eps))
  W_2 <- matrix(nrow=n_3,ncol=n_2,runif(n_3*n_2,-eps,eps))
  train=sample(1:ncol(X),floor(0.5*ncol(X)))
  Yt <- matrix(Y[train],nrow=1,ncol=length(Y[train]))
  Yv <- matrix(Y[-train],nrow=1,ncol=length(Y[-train]))
  Xt <- X[ ,train]
  Xv <- X[,-train]
  b_1 <- matrix(0,nrow=n_2,ncol=ncol(Xt))
  b_2 <- matrix(0,nrow=n_3,ncol=ncol(Xt))
  for (i in 1:nit) {
    cat("j=",j, " ","i=",i, "\n", sep = "")
    # FORWARD PROPAGATION
    Z_2 <- W_1 %*% Xt + b_1
    A_2 <- sigm(Z_2) # SIGMOID FUNCTION
    #     A_2 <- pmax(Z_2,0.01*Z_2) # Leaky ReLU function
    Z_3 <- W_2 %*% A_2 + b_2
    A_3 <- sigm(Z_3)
    #   A_3 <- Z_3 # USE THE IDENTITY FUNCTION FOR CONTINUOUS DATA
    # BACK PROPAGATION
    DZ_3 <- A_3 - Yt
    DW_2 <- (DZ_3 %*% t(A_2) / m) + delta * W_2
    Db_2 <- mean(DZ_3)
    DZ_2 <- t(W_2) %*% DZ_3 * A_2 * (1 - A_2) # sigmoid function
    # Leaky ReLU function:
    #   DZ_2 <- t(W_2)%*%DZ_3 * ifelse(Z_2 > 0,1,0.01)  
    
    DW_1 <- (DZ_2 %*% t(Xt) / m) + delta * W_1
    Db_1 <- apply(DZ_2, 1, mean)
    # GRADIENT DESCENT ON TRAINING DATA Xt
    W_1 <- W_1 - gamma * DW_1
    W_2 <- W_2 - gamma * DW_2
    b_1 <- b_1 - gamma * Db_1
    b_2 <- b_2 - gamma * Db_2
    # BELOW: ADD PENALTY TERM TO THE LOSS FUNCTION
    newcost <- cost(A_3, Yt) + (delta/2)*(sum(W_2^2)+sum(W_1^2))
    result[i,] <- c(i,newcost,DW_2[1:5],W_1[1],W_2)
    
  }
  ytrain <- as.numeric(ifelse(A_3 > 0.5, 1, 0))
  msetrain <- mean((ytrain - Yt) ^ 2)
  #  print(table(yh, Y))
  # VALIDATION STAGE WITH DATA Xv
  # A LITTLE TWIST: IN LINES 400-404 b_1 & b_2 MUST BE ADJUSTED
  # BECAUSE/IF # TRAINING RECORDS < # VALIDATING RECORDS !!!!!!
  a <- floor(0.5*ncol(X))
  b <- 0.5*ncol(X)
  bind <- function(b){cbind(b,b[,1])}
  b_1 <- if(a<b) {bind(b_1)}
  b_2 <- if(a<b) {bind(b_2)}
  Z_2 <- W_1 %*% Xv + b_1
  A_2 <- sigm(Z_2)
  Z_3 <- W_2 %*% A_2 + b_2
  A_3 <- sigm(Z_3)
  yval <- as.numeric(ifelse(A_3 > 0.5, 1, 0))
  mseval <- mean((yval - Yv) ^ 2)
  resultval[j, ] <- c(j, msetrain,mseval)
  
}
proc.time() - ptm
print(table(yval,Yv))
plot(resultval[,2],type="l",ylim=c(min(resultval[,2]),
                                   max(resultval[,3])))
lines(resultval[,3],col="red")
summary(resultval[,3])

# FIT neuralnet ON TRAINING DATA AND 
# TEST ON VALIDATING DATA
# REPEAT nrepnn TIMES
set.seed(371111)
library(neuralnet)

nrepnn <- 20

resmsenn <- matrix(data=NA,nrow=nrepnn,ncol=3)

for (i in 1:nrepnn){
  cat(i, "\n",sep="")
  train=sample(1:ncol(X),floor(0.5*ncol(X))) # ASSUMES X
  #   HAS BEEN TRANSPOSED!
  # SETTING BELOW hidden = c(5,6) FITS TWO LAYERS OF 5 NEURONS
  # IN LAYER 1 AND 6 NEURONS IN LAYER 2.  hidden =c(5) FITS
  # A SINGLE LAYER WITH 5 NEURONS
  nn = neuralnet(y~.,data=df[train,],hidden = c(5),linear.output =
                   FALSE,act.fct="logistic")
  pnv <- as.numeric(predict(nn,df[-train,]))
  predv <- ifelse(pnv > 0.5, 1, 0)
  msev <- mean((pnv-Y[-train])^2)
  pnt <- as.numeric(predict(nn,df[train,]))
  predt <- ifelse(pnt > 0.5, 1, 0)
  mset <- mean((pnt-Y[train])^2)
  resmsenn[i,] <- c(i,msev,mset)
}
plot(resmsenn[,2],type="l",ylim=c(min(resmsenn[,3]),
                                  max(resmsenn[,2])))
lines(resmsenn[,3],col="red")
summary(resmsenn[,2])
