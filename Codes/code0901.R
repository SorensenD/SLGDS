# CODE0901

# THIS CODE COLLECTS VARIOUS PIECES OF CODE FROM THE BOOK (CHAPTER 9)


# 1. GENERATES DATA y (AND MATRIX OF GENOTYPIC MARKERS X)
# 2. PERFORMS A GWAS ANALYSIS USING BONFERRONI CORRECTION
# 3. PERFORMS AN McMC ANLAYSIS AND COMPUTES BAYESIAN FDR
# 4. COMPUTES CLASSICAL FDR-BH
# 5. OBTAINS MSE OF PREDICTION USING LASSO (WITH PACKAGE GLMNET)
# 6. PERFORMS A PENALISED LOGISTIC REGRESSIOn WITH GRADIENT DESCENT
# 7. PERFORMS A PENALISED LOGISTIC REGRESSIOn WITH NEWTON-RAPHSON

rm(list=ls()) # CLEAR WORKSPACE
#install.packages("modeest")
library(modeest)

set.seed(30337)
##############################
# 1. INITIALISE AND CREATE DATA

# FREQUENCY OF 1'S IN DATA 
p0 <- 0.25

mu <- qnorm(p0)

va<-1.0 # additive variance of liability
ve<-1.0 # environmental variance
herprobit <- va/(va+ve)

nqtl<- 50
nindiv <- 2000

nmark<-1500

# Choose length of McMC chain
rep<-1000
#rep<-4000
#rep <- 10000


# Set initial proportion of QTL
piqtl<-nqtl/nmark
# parameters a and b of the beta prior
sc1=1.5
sc2=10
IDq<-sample(1:nmark,nqtl,replace=F) # from the nmark markers, choose nqtl as QTL
X<-matrix(nrow= nindiv,ncol= nmark,rbinom(nindiv*nmark,size=2,p=.5))

be<-matrix(data=0.0,nrow=nmark,ncol=1) # parameter from true model
b<-matrix(data=0.0,nrow=nmark,ncol=1) # parameter from instrumental model
theta1<-matrix(data=0.0,nrow=nmark,ncol=1) # parameter from instrumental model
xmg<-matrix(nrow=nindiv,ncol=1,rbinom(nindiv,size=2,p=0.01))
GamM<-matrix(data=NA,nrow=nmark,ncol=2)
GamP<-matrix(data=NA,nrow=nmark,ncol=1)
Xc<-matrix(data=NA,nrow=nindiv,ncol=nmark)
Xlasso<-matrix(data=NA,nrow=nindiv,ncol=nmark)
Xt<-matrix (data=NA,nrow=nindiv,ncol=nqtl) # MATRIX OF QTL GENOTYPES 
xmgc<-matrix(data=NA,nrow=nindiv,ncol=1)
xmgc<- xmg-mean(xmg)
xx<-rep(0,nmark)
zz<-rep(0,nmark)
avr<-rep(0,nmark)
truefd<-rep(0,nmark)

fd<-rep(0,nmark)
fdloc<-rep(0,nmark)
delta<-matrix(data=0,nrow=nmark,ncol=1)
deltaff<-matrix(data=0,nrow=nmark,ncol=1)

res0<-matrix(data=0.0,nrow=nmark,ncol=1)
res1<- matrix(data=0.0,nrow=nmark,ncol=1)
#delta<-matrix(numeric(0),nrow=nmark,ncol=1)
result<-matrix(data=NA,nrow=rep,ncol=11)
resultBFDR<-matrix(data=NA,nrow=rep,ncol=nmark)
resultprobtheta<-matrix(data=NA,nrow=rep,ncol=nmark)

alfa<-matrix(data=0,nrow=nmark,ncol=1)
postprob<-rep(0,nmark)
postprobff<-rep(0,nmark)

msetest<-rep(0,rep)
msevnr <- rep(0,rep)
msevnrtheta <- rep(0,rep)
genomicvalue<- rep(0,nindiv)
ypred<-matrix(data=NA,nrow=rep,ncol=nindiv/2)
cmX<-colMeans(X)
#Center matrix X
for (i in 1:nmark)
{Xc[,i]<- X[,i]-cmX[i]
}
#qr(Xc)$rank
meanXc <- apply(Xc,2,mean)
varXc<-apply(Xc,2,var)
sumvar<-sum(varXc[IDq]) # Compute 2*Sum(p(1-p)), where the Sum is over the nqtl QTL
QTLeff<-sqrt(va/sumvar) # calculate the QTL effect so that the total genetic variance is VA
be[IDq] <- QTLeff
Xt<-Xc[,IDq]
Sb<- (va/sumvar)*(6.1/4.1)
# Initialise scale parameter of scaled inverse chi-square distributions for conditional variance V
Sv<- 1*(6.1/4.1)

liab<-Xc%*%be+rnorm(nindiv,mean=0,sd=1) # liability PROBIT (excluding mu)

liab<-liab+mu # include mu
xb<-Xc%*%be
p1 <-pnorm(mu+xb) # PROBIT MODEL
y <- rbinom(nindiv,1,p1) # BINARY RECORDS
mean(y) # PROPORTION OF 1'S
sum(varXc[IDq]*QTLeff^2)
var(Xc%*%be)
############# END OF GENERATION OF DATA ####################################
############################################################################
# 2. GWAS USING A LOGISTIC REGRESSION WITH BONFERRONI CORRECTION
GWAS=matrix(nrow=ncol(Xc),ncol=4)
colnames(GWAS)=c('estimate','SE','t-value','p-value')
for(i in 1:ncol(Xc)){
  fm=glm(y~Xc[,i],family = binomial(link = "probit"))
  
  #fm=glm(y~Xc[,i])
  GWAS[i,]= summary(fm)$coef[2,]
}
#setwd("C:/Users/au223137/Dropbox/Rsessions/SummerCourse/OutputR")
#pdf("C:/Users/au223137/Dropbox/Rsessions/SummerCourse/OutputR/GWAS2000-1500Bonf.pdf")
plot(-log10(GWAS[,4]),type='o',ylab='-log10-pValue',xlab='Marker ID',cex=0.5,col=4,cex.lab=1.3)
abline(h=-log10(0.05/nmark),lty=2,col=2,lwd=2)
#dev.off()
cat('Bonferoni',-log10(0.05/nmark),'\n')
GWASdetct<-which(-log10(GWAS[,4]) > -log10(0.05/nmark))
GWASdetctUnc<-which(-log10(GWAS[,4]) > -log10(0.05))
discovBonf <- length(GWASdetct) # DISCOVERY SET BONFERRONI
discovUnc <- length(GWASdetctUnc) # DISCOVERY SET UNCORRECTED
################################################################################
# 3. CODE FOR THE GIBBS SAMPLER
# USES DATA y AND MARKER MATRIX Xc AS INPUT

V <- 1
# Initialise Vg, genomic variance, which initially is set equal to the additive genetic variance
Vg<-va
sumvarqtl<-sum(apply(as.matrix(Xc[,IDq]),2,var))
# start Gibbs chain
#marker effect variance is initially set initially equal to genomic variance divided by 
# the sum of the variance of the columns of Xc corresponding to QTLs
Vb<-Vg/sumvarqtl
maxlogods<- 0
# **************************************
#train=sample(1:nrow(Xc),nrow(Xc)/2) # FOR PREDICTION
train=sample(1:nrow(Xc),nrow(Xc)) # FOR DETECTION

test=(-train)
y.test=y[test]
y.train<-y[train]
Xtest <- Xc[test,]
Xtrain <- Xc[train,]
nindivt<-nrow(Xc[train,])
# Compute col vector xx with sum of squared terms involved in t(Xc)%*%Xc for the training data
for (i in 1:nmark){xx[i]<-crossprod(Xc[train,i])}
ptm<-proc.time()
# **************************************
u<-mean(liab)
# *****************************
res1<-liab-u
# ****************************

for (i in 1:rep)
{
  print(i)
  # sample liabilities liab from truncated normals
  meanliab <- u + Xc[train,]%*%b
  sd <- 1
  interm <- y.train*pnorm(0,mean=meanliab,sd=sd)+runif(length(y.train))*
    (pnorm(0,mean=meanliab,sd=sd)*(1-y.train) + (1-pnorm(0,mean=meanliab,sd=sd))*y.train)
  liab <- qnorm(interm,mean=meanliab,sd=sd)
  res1 <- liab - meanliab
  # sample u - the intercept 
  res1<-res1+u
  avu<-sum(res1)/nindivt
  varu<-V/nindivt
  u<-rnorm(1,mean=avu,sd=sqrt(varu))
  k<-V/Vb
  res1<-res1-u
  for (j in 1:nmark){
    res0<-res1+Xc[train,j]*b[j,1]
    rss0<-sum(res0*res0)
    if(delta[j,1] == 1) 
    {
      alfahat<- (t(Xc[train,j])%*%(res0))/(xx[j]+k)
      varalfa<-(1/(xx[j]+k))*V
      alfa[j,1]<- rnorm(1,mean=alfahat,sd=sqrt(varalfa))
    } else {
      alfa[j,1]<-rnorm(1,mean=0,sd=sqrt(Vb)) 
    }
    # sample deltas
    res1temp<-res0-Xc[train,j]*alfa[j,1]
    rss1temp<-sum(res1temp*res1temp)
    logodds<- (1/(2*V))*(rss0-rss1temp)-(log(1-piqtl)-log(piqtl))
    theta1[j]<-exp(logodds)/(1+exp(logodds))
    deltaff[j,1] <- rbinom(1,1,theta1[j])
    
    uniform<-runif(1)
    logunif<-log(uniform/(1-uniform))
    if(logunif<=logodds){delta[j,1]<-1}else{delta[j,1]<-0}
    # update the regression    
    b[j,1]<-delta[j,1]*alfa[j,1]
    postprob[j]<-postprob[j]+delta[j,1]
    res1<-res0-as.matrix(Xc[train,j])*b[j,1]
    postprobff[j] <- postprobff[j] + deltaff[j,1]
  }
  
  # Update nqtl = sum(delta[,1])
  rss1<-sum(res1*res1)
  nqtl<-sum(delta[,1])
  # sample piqtl
  piqtl<-rbeta(1,nqtl+sc1,nmark-nqtl+sc2)
  
  # update variance of a marker effect (Vb)
  df<- 4.1+nmark
  Vb<- (sum(alfa*alfa) + Sb*4.1)/rchisq(1,df)
  # update conditional variance V (that is, Var[y|u,X,b] of the instrumental model)
  #  df<- 4.1+nindivt
  #  V<- (rss1 + Sv*4.1)/rchisq(1,df)
  genomicvalue<- Xc[test,]%*%b # Pred genomic values for the current iteration for testing data
  probval <- pnorm(u+genomicvalue)
  msetest[i]<-mean((probval-y.test)^2) # MSE using testing data based on probabilities
  y_predval <- as.numeric(ifelse(probval > 0.5, 1, 0))
  #  y_predvaltheta <- as.numeric(ifelse(theta1 > 0.5,1,0))
  msevnr[i] <- mean((y_predval - y.test) ^ 2) # MSE using testing data based on Ypredval
  #  msevnrtheta[i] <- mean((y_predvaltheta - y.test) ^ 2) # MSE using testing data based on Ypredvaltheta
  resultprobtheta[i,] <- theta1
  result[i,]<- c(i,piqtl,nqtl/nmark,Vb,V,delta[IDq[1:3]],u,msetest[i],msevnr[i])
}
proc.time()-ptm

postprob <- postprob/rep
postprobff <- postprobff/rep

#postprobability<-postprob/rep
#setwd("C:/Users/au223137/Dropbox/Rsessions/SummerCourse/OutputR")
#pdf("C:/Users/au223137/Dropbox/Rsessions/SummerCourse/OutputR/PostProb2000-1500.pdf")
plot(postprob,type='o',ylab='PostProb',xlab='Marker ID',cex=.5,col=4,cex.lab=1.3)
#axis(1,seq(0,nmark,200))

abline(h=0.8,lty=2,col=2,lwd=2)
####################################################################################
# FUNCTION TO GENERATE HISTOGRAM OF BAYESIAN FDR USING GIBBS OUTPUT
# INPUT:
# 1. resultprobtheta: A FILE WITH DRAWS FROM (7.57)
# 2. fd: THE TRUE PROPORTION OF FALSE DISCOVERIES OBSERVED IN THE SAMPLE (see line 402)
# 3. prob: THE USER-USED THRESHOLD THAT DEFINES THE DISCOVERY SET
# 4. postprob: VECTOR OF POSTERIOR PROBABILITIES 
#   (OF BELONGING TO THE SLAB COMPONENT OF THE MIXURE) FOR EACH GENETIC MARKER
fdrhist <- function(resultprobtheta,truefd,prob){
  discset <- which(postprob > prob)
  fdisc <- apply(1-(resultprobtheta[,discset]),1,mean)
  # MEAN AND POSTERIOR INTERVAL FOR BAYES FDR:  
  avfdis <- mean(fdisc)
  quantilefdis <- quantile(fdisc,c(0.025,0.975))
  #setwd("C:/Users/au223137/Dropbox/Rsessions/SummerCourse/OutputR")
  #pdf("C:/Users/au223137/Dropbox/Rsessions/SummerCourse/OutputR/histfdrdisc08.pdf")
  #pdf("C:/Users/au223137/Dropbox/Rsessions/SummerCourse/OutputR/histfdrdisc056.pdf")
 # HISTOGRAM OF McMC ESTIMATE OF POSTERIOR DISTRIBUTION OF BAYES FDR:  
  hist(fdisc,breaks=40,xlab='McMC-Bayes FDR',main=NULL, freq=FALSE,cex.lab=1.3)
  ### TRUE FDR: 
  abline(v=truefd[length(discset)],col="red",lwd=2) # TRUE FDR REALISED IN SAMPLE
 # dev.off()
  return <- c(avfdis,quantilefdis,length(discset),truefd[length(discset)]*length(discset))
}
out <- fdrhist(resultprobtheta,truefd,0.8)
out <- fdrhist(resultprobtheta,truefd,0.56)

out
#################################################################################
#################################################################################
#   4. FDR-BH (BenjaminiHochberg)
# USES OUTPUT FROM THE GWAS ANALYSIS AS INPUT (P-VALUES)
bID<-seq(1:nmark)
trueH1<-IDq
trueH0<-bID[-IDq]
# ALFA HERE IS THE q IN BOOK 
# (chosen value of expected proportion of false discoveries)
#alfa<-0.05
alfa<-0.15 

adjpv<-vector()
qvalue<-vector()
qvbh<-vector()
pval<-GWAS[,4]
bhat<-GWAS[,1] # NOT NEEDED
compldata<-data.frame(bID,bhat,pval)
sortcd<-compldata[order(compldata$pval),]
for(i in nmark:1){
  adjpv[i]<-(i/nmark)*alfa
}
compldata<-data.frame(sortcd,adjpv)
for (i in nmark:1){
  if(compldata$pval[i] <= compldata$adjpv[i]){
    print(i)
    break
  }
}
sizediscovset<-i # SIZE OF DISCOVERY SET
signif<-compldata[1:i,]
nslb<-i+1
notsignif<-compldata[nslb:nmark,]
# NO. TRUE FALE DISCOVERIES
sizediscovset
fdiscov<-length(setdiff(signif$bID,trueH1))
fdiscov # NO. TRUE FALSE DISCOVERIES
fdiscov/sizediscovset # TRUE OBSERVED PROPORTION OF FD IN SAMPLE 
# EXPECTED PROPORTION ALFA (CONCEPTUALLY AVERAGED OVER SAMPLES)
alfa
#################################################################
#################################################################
# 5. Lasso solutions using package glmnet
#install.packages("glmnet", .libPaths()[1])
#ptm<-proc.time()
# USES DATA y AND MARKER MATRIX Xc AS INPUT
set.seed(3337)
library(glmnet)
n<-nrow(Xc)
Xlasso<-Xc
mu_y <- mu
train=sample(1:nrow(Xc),nrow(Xc)/2)
test=(-train)
y.test=y[test]
y.train<-y[train]
Xtest <- Xc[test,]
Xtrain <- Xc[train,]
# STEP 1
cv.out=cv.glmnet(Xlasso[train,],y[train],alpha=1,intercept=TRUE,family="binomial",type="class")
plot(cv.out)
bestlam <- cv.out$lambda.min
bestlam
# NUMBER OF NON-ZERO COVARIATES
length(which(as.vector(coef(cv.out,s=bestlam))!=0))
# STEP 2
# OBTAIN PREDICTIONS BASED ON CLASS LABELS
fm.predclass <- predict(cv.out,s=bestlam,newx=Xtest,family="binomial",type="class")
# OBTAIN PREDICTIONS BASED ON CLASS PROBABILITIES (Brier score)
fm.predresp <- predict(cv.out,s=bestlam,newx=Xtest,family="binomial",type="response")
# ERROR RATE (CLASS LABELS)
mean((as.numeric(fm.predclass)-y.test)^2)
# ERROR RATE (RESPONSE LABELS)
mean((as.numeric(fm.predresp)-y.test)^2)
# ERROR RATE NULL MODEL (liability = mu + e)
mean(y.test)
################################################################################
################################################################################
# 6. PENALISED LOGISTIC REGRESSION - GRADIENT DESCENT

# USES DATA y AND MARKER MATRIX Xc AS INPUT

nrepgd <- 2000 # NUMBER OF GD ITERATIONS (USE AT LEAST 5000 OR BETTER, 10000)

replic <- 1 # TO ESTIMATE VARIATION OF MSE OVER REPLICATES
X <- Xc
beta <- rep(0.0,ncol(Xc))
miu <- 0.0

#lambda <- 0.0 # ZERO PENALTY !!!!!!!!
lambda <- 0.4

gama <- 0.00085# Learning rate
c <- rep(0,nrepgd)

newcostv <- rep(0,nrepgd)

resulttv <- matrix(data=NA, nrow=nrepgd,ncol=8)

res <- matrix(data=NA, nrow=replic,ncol=9)

msev <- rep(NA,replic)

#################################################
# FUNCTION TO COMPUTE Pr[Y = 1]
prob1 <- function(miu,beta,X){
  pr <- exp(miu+X%*%beta)/(1+exp(miu+X%*%beta))
}
# FUNCTION TO COMPUTE THE LOSS FUNCTION
cost <- function(miu,beta,X,y)
{-sum(y*(miu+X%*%beta) - log(1 + exp(miu+X%*%beta)))
  + crossprod(beta)*(lambda/2)}
c[1]<- cost(miu,beta,X,y)
##################################################################
#########               GRADIENT DESCENT        ##################
### FIT MODEL TO TRAINING DATA AND TEST IN VALIDATING DATA ###

set.seed(771311)
numone <- sum(y)
numzero <- length(y) - numone
nindiv <- length(y)
ptm <- proc.time()
for (i in 1:replic){
  cat(i, "\n",sep="")
  train=sample(1:nrow(X),floor(0.5*nrow(X))) # FOR PREDICTION
  #  train=sample(1:nrow(Xc),nrow(Xc)) # FOR ESTIMATION
  
  Xtrain <- X[train,]
  Xval <- X[-train,]
  ytrain <- y[train]
  yval <- y[-train]
  miu <- 0.0
  beta <-  rep(0.0,ncol(Xtrain))
  for(j in 1:nrepgd){
    cat(j, "\n",sep="")
    fdmiu <- -sum(ytrain - prob1(miu, beta,Xtrain))
    fdbeta <- -t(Xtrain) %*%
      (ytrain - prob1(miu, beta,Xtrain)) + lambda * beta
    fd <- matrix(c(fdmiu, fdbeta), nrow = length(beta)+1, ncol = 1)
    sol0 <- matrix(c(miu, beta), nrow = length(beta)+1, ncol = 1)
    sol1 <- sol0 - gama * fd
    miu <- sol1[1,1]
    beta <- sol1[-1,1]
    newcostv[j] <- cost(miu, beta,Xtrain,ytrain)
    resulttv[j,] <- c(j,newcostv[j],miu,beta[1:5])
  }
  probval <- prob1(miu,beta,Xval)
  y_predval <- as.numeric(ifelse(probval > 0.5,1,0))
  msev[i] <- mean((y_predval-yval)^2)
  res[i,] <- c(i,j,newcostv[j],miu,beta[1:5])
}
proc.time()-ptm
tail(resulttv)

tail(res)
summary(msev) # SUMMARY OF MISCLASSIFICATION ACROSS REPLICATES
#############################################################################
#############################################################################
# 7. PENALISED LOGISTIC REGRESSION - NEWTON-RAPHSON
# USES DATA y AND MARKER MATRIX Xc AS INPUT

# INITILISATION
nitnr <- 15 # N-R ITERATIONS PER REPLICATE
nrep <- 2 # REPLICATES (DISPLAYS VARIANCE OF MSE) 
lambda <- 0.1
miu <- 0
b<-  rep(0.1,ncol(X))
newcostvnr <- rep(0,nrep)
res <- matrix(data=NA, nrow=nrep,ncol=11)
resulttvnr <- matrix(data=NA, nrow=nitnr,ncol=8)
msevnr <- rep(NA,nrep)
msebrier <- rep(NA,nrep)
############################################################
prob1 <- function(miu,b,Xc){
  pr <- exp(miu+Xc%*%b)/(1+exp(miu+Xc%*%b))
}
cost <- function(miu,b,Xc,y){-sum(y*(miu+Xc%*%b)-log(1+exp(miu+Xc%*%b))) 
  + lambda*crossprod(b)}

#   FIT MODEL TO TRAINING DATA AND TEST IN VALIDATING DATA ###
set.seed(7713)
ptm <- proc.time()

for (i in 1:nrep) {
  cat(i, "\n",sep="")
   train=sample(1:nrow(Xc),floor(0.5*nrow(Xc)))
#  train=sample(1:nrow(Xc),nrow(Xc))
  
  Xtrain <- Xc[train,]  
  Xval <- Xc[-train,]
  ytrain <- y[train]
  yval <- y[-train]
  delta <- diag(c(0,rep(lambda,ncol(Xtrain))))
  M <- cbind(0,diag(1,nrow=ncol(Xtrain)))
  M <- rbind(0,M)
  ###### START VALUES FOR MIU AND ALFA ################
  miu <- 0.0
  b <- rep(0.0, ncol(Xtrain))
  ###############################################
  W <- matrix(data = 0,nrow = ncol(Xtrain),ncol = ncol(Xtrain))
  for (j in 1:nitnr) {
    cat(j, "\n",sep="")
    fdmiu <- -sum(ytrain - prob1(miu, b, Xtrain))
    fb <-
      -t(Xtrain) %*% (ytrain - prob1(miu, b, Xtrain)) + lambda *  b
    fd <- matrix(c(fdmiu, fb), nrow = length(b) + 1, ncol = 1)
    W <-
      diag(c(prob1(miu, b, Xtrain) * (1 - prob1(miu, b, Xtrain))))
    Z <- cbind(1, Xtrain)
    zwz <- t(Z) %*% W %*% Z
    #    zwz <- crossprod(Z,W)%*%Z
    #    zwz <- t(Z)%*%crossprod(W,Z)
    LHS <- zwz + lambda * M
    RHS <- -fd
    sol0 <- matrix(c(miu, b), nrow = length(b) + 1, ncol = 1)
    sol1 <- sol0 + solve(LHS,RHS)
    
    miu <- sol1[1, 1]
    b <- sol1[-1, 1]
    newcostvnr[j] <- cost(miu, b, Xtrain, ytrain)
    resulttvnr[j, ] <- c(j, newcostvnr[j], miu, b[1:5])
  }
  probval <- prob1(miu, b, Xval)
  probtrain <- prob1(miu, b, Xtrain)
  msebrier[i]<-mean((probval-yval)^2) # MSE using validating data based on probabilities
  
  y_predval <- as.numeric(ifelse(probval > 0.5, 1, 0))
  y_predtrain <- as.numeric(ifelse(probtrain > 0.5, 1, 0))
  
  msevnr[i] <- mean((y_predval - yval) ^ 2) # MSE using validating data BASED ON Y(0,1) 
  res[i,] <- c(i,j,newcostvnr[j],miu,b[1:5],msevnr[i],msebrier[i])  
}
proc.time()-ptm
tail(resulttvnr)
tail(res)
summary(msevnr) # SUMMARY OF MISCLASSIFICATION ACROSS REPLICATES


