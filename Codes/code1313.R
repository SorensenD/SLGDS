# CODE1313
# BAYES PROBLEMS II Exercise 5. GENOMIC MODEL
# TWO OPTIONS FOR PRIOR ASSUMPTIONS OF VARIANCE COMPONENTS (SEE LINE 34)
# DATA BASED ON GENOMIC MODEL; OBTAIN SVD OF WW'(1/m)
########################################################################
#set.seed(12345)
#nindiv<-500
#nmark<-1000
#vgs<-10
#ves<-25
#rep<-4000
# RESULTS IN AN AVERAGE MC POSTERIOR MEAN OF Vg EQUAL TO APPROX 10
######################################################################
rm(list=ls()) # CLEAR WORKSPACE
set.seed(12345)
nindiv<-500
nmark<-1000
#nindiv <- 1000
#nmark <- 10000
nt <- nindiv*nmark
# GENERATE MARKER MATRIX FROM BINOMIAL DISTRIBUTION
X<-matrix(nrow= nindiv,ncol= nmark,rbinom(n=nt,size=2,p=.5))
stdev <- matrix(data=NA,nrow= nmark,ncol=1)
W <- matrix(data=NA,nrow= nindiv,ncol=nmark)
U <- matrix(data=NA,nrow= nindiv,ncol= nindiv)
G<-matrix(data=NA,nrow= nindiv,ncol= nindiv)
cm <- colMeans(X)
#CHOOSE VALUE FOR GENOMIC VARIANCE vgs
vgs<-10
#CHOOSE VALUE FOR ENVIRONMENTAL VARIANCE ves
ves<-25

########################################################################
### CODE r=1 IF VARIANCE COMPONENTS ASSIGNED INPROPER UNIFORM PRIORS
### ELSE: VARIANCE COMPONENTS ASSIGNED SCALED INVERTED CHI SQUARE DISTRIBUTIONS
#########################################################################
r <- 2
if(r!=1){
  ## HYPERPARAMETERS OF THE SCALED INVERTED PRIORS FOR va AND ve
  nua <- 4.1
  nue <- 4.1
  Sa <- 10
  Se <- 30
  ## THIS GENERATES PRIOR MODES FOR va AND ve equal to 10*(4.1/6.1) = 6.7
  ## and 10*(4.1/6.1) = 20.2
}
########################################################################


# CREATE MATRIX OF STANDARDISED MARKER GENOTYPE CODES
for (i in 1:nmark)
{
  W[,i] <-( X[,i]-cm[i]) / sd(X[,i])
}
# GENERATE nindiv GENOMIC VAL FROM N(0,(1/nmark)WW'*vgs) 
#NOTE: MARKER VALUES ARE DRAWS FROM N(0,I sqrt(vgs/nmark))
g <- (1/sqrt(nmark))*W%*%rnorm(nmark,mean=0,sd=sqrt(vgs))
# GENERATE nindiv PHENOTYPES WITH MEAN 0, 
#  VAR=vgs+ves, HERITABILITY=vgs/(vgs+ves)
e<- rnorm(nindiv,mean=0,sd=sqrt(ves))
y <- g+ e
# GENOMIC RELATIONSHIP MATRIX G
G <- (1/nmark)*W%*%t(W)
# SVD OF G
EVD <- eigen(G)
names(EVD)
#head(EVD$values)
U <- EVD$vector
tU<-t(U)
val <- EVD$values
val[length(y)] <-0
D <- diag(val,nrow=nindiv)
#Dp IS A VECTOR WITH NON-ZERO EIGENVALUES
Dp<-c(val[1:nindiv-1])
#INITIALISE Ve
Ve<-5
#INITIALISE Vg
Vg<-5
#INITIALISE k
k<-Ve/Vg
#INITIALISE VECTOR ALFA
alfa<-rep(0,nindiv)
# CHOOSE LENGTH OF GIBBS CHAIN
rep<-4000
#INITIALISE result
result<-matrix(data=NA,nrow=rep,ncol=7)
# START GIBBS CHAIN
ptm<-proc.time()

for (i in 1:rep)
{
   print (i)
  # SAMPLE mu
  avmu<-sum(y-U%*%alfa)/nindiv
  varmu<-Ve/nindiv
  mu<-rnorm(1,mean=avmu,sd=sqrt(varmu))
  #mu<-0
  # SAMPLE alfa1 (VECTOR OF LENGTH nindiv-1)
  meanalfa1<-(Dp/(Dp+k))*tU[1:nindiv-1,]%*%(y-mu)
  varalfa1<-((Dp)/(Dp+k))*Ve
  alfa1<-rnorm(nindiv-1,meanalfa1,sqrt(varalfa1))
  alfa<-c(alfa1,0)
  ystar<-y-mu-U%*%alfa
  
  # SAMPLE Vg and Ve
  if(r==1) #  IMPROPER UNIFOR PRIORS FOR Vg, Ve
  { scVg<-sum(alfa1*alfa1*(1/Dp))
  scVe<-sum(ystar*ystar)
  Vg<-scVg/rchisq(1,nindiv-3)
  Ve<-scVe/rchisq(1,nindiv-2)
  }
  if (r!=1) # SCALED INVERTED CHI SQUARE PRIORS FOR Vg,Ve
  { scVg<-sum(alfa1*alfa1*(1/Dp)) + nua*Sa
  scVe<-sum(ystar*ystar) + nue*Se
  Vg<-scVg/rchisq(1,nindiv-1+nua)
  Ve<-scVe/rchisq(1,nindiv + nue)
  }
  
  # SAMPLE Vg and Ve

  k<-Ve/Vg
  result[i,]<-c(i,mu,Vg,Ve,Vg/(Vg+Ve),1/k,mean(alfa*alfa))
  #  print(result[i,])
}
proc.time()-ptm
apply(result,2,mean)
# FUNCTION LOGLIK TO CONSTRUCT THE LOGLIKELIHOOD
# TO COMPARE THE BAYESIAN RESULTS WITH OPTIM
#NOTE: k IN THE LOGLIKELIHOOD IS Vg/Ve
loglik<-function(data,par)
{
  mu<-par[1]
  ve<-par[2]
  k<-par[3]
  dkiinv<-diag(1/((val*k)+1),nrow=nindiv,ncol=nindiv)
  ll<- -0.5*(length(y)*log(ve)+sum(log((val*k)+1))+
               (1/ve)*t(y-mu)%*%U%*%dkiinv%*%tU%*%(y-mu))
  return(-ll)
}
result1<-optim(par=c(0.5,12,2),loglik,data=y,hessian=TRUE)
result1$par
vgml <- result1$par[2]*result1$par[3]
solve(result1$hessian)
meanmu <- mean(result[,2])
meanmu
#apply(result[,2:6],2,mean)
cimu <- quantile(result[,2],c(0.025,0.975))
cimu
meanvg <- mean(result[,3])
meanvg
civg <- quantile(result[,3],c(0.025,0.975))
civg
meanve <- mean(result[,4])
meanve
cive <- quantile(result[,4],c(0.025,0.975))
cive
meanher <- mean(result[,5])
meanher
ciher <- quantile(result[,5],c(0.025,0.975))
ciher
meank <- mean(result[,6])
meank
cik <- quantile(result[,6],c(0.025,0.975))
cik
meanalfasq <- mean(result[,7])
meanalfasq
cialfasq <- quantile(result[,7],c(0.025,0.975))
cialfasq

