# CODE1003
# FITS A SPIKE AND SLAB MODEL TO GAUSSIAN DATA
# CALCULATES BAYESIAN LEAVE-ONE-OUT CROSS-VALIDATIONS

rm(list=ls()) # CLEAR WORKSPACE

set.seed(3033711)
va <- 10 # additive genetic variance

ve<-30 # environmental variance

her<-va/(va+ve)

nindiv <- 2500 #(nindiv is the # individuals in training and in validating data)
nmark<-5000
Nqtl <- 25
mu_y<-0

# Set initial proportion of QTL
pi<-Nqtl/nmark
# Choose length of McMC chain
rep <- 2000

# parameters a and b of the beta prior
sc1=1.5
sc2=10
 X<-matrix(nrow= nindiv,ncol= nmark,rbinom(nindiv*nmark,size=2,p=.5))

be<-matrix(data=0.0,nrow=nmark,ncol=1) # parameter from true model
b<-matrix(data=0.0,nrow=nmark,ncol=1) # parameter from instrumental model
theta1<-matrix(data=0.0,nrow=nmark,ncol=1) # parameter from instrumental model
Xc<-matrix(data=NA,nrow=nindiv,ncol=nmark)
Xt<-matrix (data=NA,nrow=nindiv,ncol=Nqtl) # MATRIX OF QTL GENOTYPES 
xx<-rep(0,nmark)
avr<-rep(0,nmark)
fd<-rep(0,nmark)
fdloc<-rep(0,nmark)
truefd<-rep(0,nmark)

delta<-matrix(data=0,nrow=nmark,ncol=1)
res0<-matrix(data=0.0,nrow=nmark,ncol=1)
res1<- matrix(data=0.0,nrow=nmark,ncol=1)
result<-matrix(data=NA,nrow=rep,ncol=(2*nmark)+10)
resultprobtheta<-matrix(data=NA,nrow=rep,ncol=nmark)
resultMSE<-matrix(data=NA,nrow=rep,ncol=4)
keepZvb <- matrix(data=NA,nrow=rep, ncol=nindiv)
keepZtb <- matrix(data=NA,nrow=rep, ncol=nindiv)


alfa<-matrix(data=0,nrow=nmark,ncol=1)
postprob<-rep(0,nmark)
mseval<-rep(0,rep)
msetrain <- rep(0,rep)


genomicvalval<- rep(0,nindiv)
genomicvaltrain<- rep(0,nindiv)

ypred<-matrix(data=NA,nrow=rep,ncol=nindiv)

cmX<-colMeans(X)
#Center matrix X
for (i in 1:nmark)
{Xc[,i]<- X[,i]-cmX[i]
}
#qr(Xc)$rank
meanXc <- apply(Xc,2,mean)
#meanXc <- colMeans(Xc)
varXc<-apply(Xc,2,var)
IDq<-sample(1:nmark,Nqtl,replace=F) # from the nmark markers, choose Nqtl as QTL 
sumvar<-sum(varXc[IDq]) # Compute 2*Sum(p(1-p)), where the Sum is over the Nqtl
# QTL
QTLeff <- sqrt(va/sumvar) # calculate the QTL effect so that the total 
# genetic variance is VA
be[IDq]<-QTLeff # the only b's that are not zero are those associated with QTL.
Xt<-Xc[,IDq]

# Initialise scale parameters of scaled inverse chi-square distributions for 
# marker effect variance, Vg/Nqtl
# Using the mode, choose scale Sg = (va/Nqtl)*(vg+2/vg), and set degrees of 
# freedom vg = 4.1
Sb<- (va/sumvar)*(6.1/4.1)
# Initialise scale parameter of scaled inverse chi-square distributions for 
# conditional variance V
Sv<- (va+ve)*(1-her)*(6.1/4.1)

et <- rnorm(nindiv,0,sqrt(ve))
ev <- rnorm(nindiv,0,sqrt(ve))
var(et)
#######################################################
######### RESIDUAL VARIANCE LOG-NORMALLY DISTRIBUTED mu=1.2, sd=0.2 
####### these parameters generate a mean of approx 4.95 and var of approx 30
#et <- rlnorm(nindiv,mean=1.2,sd=sqrt(0.8))
#ev <- rlnorm(nindiv,mean=1.2,sd=sqrt(0.8))
#########################################################

xbe <- Xc%*%be
var(xbe)
yt <- mu_y + xbe + et
yv <- mu_y + xbe + ev
#y<-y+mu_y
sum(varXc[IDq]*QTLeff^2)
var(Xc%*%be)

################################################################################
##########################################################
train <- seq(1:nrow(Xc))
val <- train
y.train <- yt
y.val <- yv
#########################################################
# PERFORM A GWAS ON THE TRAINING DATA
GWAS=matrix(nrow=ncol(Xc),ncol=4)
colnames(GWAS)=c('estimate','SE','t-value','p-value')
for(i in 1:ncol(Xc)){
  fm=lm(y.train~Xc[train,i])
  GWAS[i,]= summary(fm)$coef[2,]
}
plot(-log10(GWAS[,4]),type='o',ylab='-log10-pValue',xlab='Marker Label',cex=.5,
     col=4,cex.lab=1.3)

abline(h=-log10(0.05/nmark),lty=2,col=2,lwd=2)
points(IDq,-log10(GWAS[IDq,4]),pch="*",cex=2,col="orange")
cat('Bonferoni',-log10(0.05/nmark),'\n')
GWASdetctBon<-which(-log10(GWAS[,4]) > -log10(0.05/nmark))
GWASdetctUnc<-which(-log10(GWAS[,4]) > -log10(0.05))

GWASdetctBon<-which(-log10(GWAS[,4]) > -log10(0.05/nmark))

true<-IDq
idmark<-seq(1:nmark)
fals<-idmark[-IDq]
discov<-GWASdetctBon
length(discov)
accepth0<-which(-log10(GWAS[,4]) < -log10(0.05/nmark))

fdisc<-length(setdiff(discov,true))

fdisc # number of false discoveries


######################################
# GIBBS SAMPLING ALGORITHM
# Fit the instrumental model y = 1u + Xb + r
# INITIALISATION
# Initialise conditional variance V = Var(y)[1 - Her], V = Var(r)
V<-(va+ve)*(1-her)
# Initialise Vg = Var(y)*Her, genomic variance, which initially is set equal to 
# the additive genetic variance
Vg<-(va+ve)*her
sumvarqtl<-sum(apply(Xc[,IDq],2,var))
# start Gibbs chain
#marker effect variance is initially set initially equal to genomic variance 
# divided by the sum of the variance of the columns of Xc corresponding to QTLs
Vb<-Vg/sumvarqtl
maxlogods<- 0
#######################################################

#########################################################
vary <- var(y.train)
nindivt<-nrow(Xc[train,])
# Compute col vector xx with sum of squared terms involved in t(Xc)%*%Xc for 
# the training data
for (i in 1:nmark){xx[i]<-crossprod(Xc[train,i])}
ptm<-proc.time()
# **************************************
u<-mean(y.train)
# *****************************
res1<-y.train-u
# ****************************

for (i in 1:rep)
{
 print(i)
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
# sample alfa_s from the FCD
      
      alfa[j,1]<- rnorm(1,mean=alfahat,sd=sqrt(varalfa))
    }
    else
    {
      alfa[j,1]<-rnorm(1,mean=0,sd=sqrt(Vb)) 
    }
# sample deltas
    res1temp<-res0-Xc[train,j]*alfa[j,1]
    rss1temp<-sum(res1temp*res1temp)
    logodds<- (1/(2*V))*(rss0-rss1temp)-(log(1-pi)-log(pi))
    theta1[j]<-exp(logodds)/(1+exp(logodds))
    uniform<-runif(1)
    logunif<-log(uniform/(1-uniform))
    if(logunif<=logodds){delta[j,1]<-1}else{delta[j,1]<-0}
# update the regression    
    b[j,1]<-delta[j,1]*alfa[j,1]
    res1<-res0-Xc[train,j]*b[j,1]
    postprob[j]<-postprob[j]+delta[j,1]
  }
  Vgenomic <- var(Xc[train,]%*%b)
# Update nqtl = sum(delta[,1])
  rss1<-sum(res1*res1)
  nqtl<-sum(delta[,1])
# sample pi
  pi<-rbeta(1,nqtl+sc1,nmark-nqtl+sc2)

# update variance of a marker effect (Vb)
  df<- 4.1+nmark
  Vb<- (sum(alfa*alfa) + Sb*4.1)/rchisq(1,df)
# update conditional variance V (that is, Var[y|u,X,b] of the instrumental model)
  df<- 4.1+nindivt
  V<- (rss1 + Sv*4.1)/rchisq(1,df)
  resultprobtheta[i,] <- theta1
  ind <- which(theta1 > 0.5) # FIND POST PROB > 0.5
  nzero <- which(b>0)
  genomicvalval<- Xc[val,]%*%b # Pred genomic values for the current iteration 
  # for validating data
  genomicvaltrain<- Xc[train,]%*%b # Pred genomic values for the current 
  # iteration for training data
  ystarval <- u+ genomicvalval + rnorm(length(y.val),0,sqrt(V))
  msetrain2<-mean((u+genomicvaltrain-y.train)^2) # MSEt2 using training data
  mseval2<-mean((u+genomicvalval-y.val)^2) # MSEv2 using validating data
  
  mseval3 <- mean((ystarval-y.val)^2) # MSEv3 using validating data
  msetrain3 <- mean((ystarval-yt)^2) # MSEv3 using training data
  
  ##############################################################################
  
  
  result[i,]<- c(i,pi,Vgenomic,nqtl,Vb,V,delta[IDq[1:3]],u,b[,1],theta1)
  resultMSE[i,] <- c(msetrain2,mseval2,mseval3,msetrain3)
  ##############################################
  # STORE u+Xb TO COMPUTE McMC ESTIMATE OF Avr(Var(Xi b_star))
  keepZvb[i,] <- u+genomicvalval
  keepZtb[i,] <- u+genomicvaltrain
}
proc.time()-ptm
postprob<-postprob/rep
plot(postprob,type='o',ylab='PostProb',xlab='Marker Label',cex=.5,col=4,
     cex.lab=1.3)

abline(h=0.5,lty=2,col=2,lwd=2)
points(IDq,postprob[IDq],pch="*",cex=2,col="orange")
################################################################################
### Genomic variance 
plot(result[,3])
mean(result[,3])
################################################################################
##################    COMPUTE LOOCV BAYES ######################################
sumloo <- 0
sumdg <- 0
sumavrb <- 0
sumXb <- 0
for(i in 1:nindiv){
 w <- 1/dnorm(y.train[i],keepZvb[,i],sqrt(result[,6])) 
 w <- w/sum(w)
 index <-sample(1:rep,rep,replace=T,prob=w)
 yhatminusi <- rnorm(rep,keepZvb[index,i],sqrt(result[index,6])) 
 sumloo <- sumloo + (y.train[i]-yhatminusi)^2
 expvalueXbminusi <- w%*%keepZvb[,i]
 sumdg <- sumdg + (y.train[i]-expvalueXbminusi)^2
 sumavrb <- sumavrb + (y.train[i]-mean(keepZvb[index,i]))^2 # (y - Xb^hat)^2; 
 # Xb^hat: post mean
 sumXb <- sumXb + (y.train[i]-keepZvb[index,i])^2 # (y - Xb)^2; b~[b|y_i]
 
}
mselooBayes <- sumloo/length(y.train) # estimate based on y^hat = y^*~[y|y_t]

hist(mselooBayes,breaks=20,xlab=NULL,main=NULL)
summary(mselooBayes)
quantile(mselooBayes,c(0.025,0.975))
mselooXb <- sumXb/length(y.train) # estimate based on y^hat = xb, b~[b|y_t]
summary(mselooXb)
quantile(mselooXb,c(0.025,0.975))
hist(mselooXb,breaks=20,xlab=NULL,main=NULL)


mseloodg <- sumdg/(length(y.train)) # estimate of mseloodg (Gianola-Schön)
mseloodg
mselooavrb <- sumavrb/length(y.train) # McMC estimate of mseloodg (Gianola-Schön)
mselooavrb
hist(mselooBayes,breaks=20,xlab=NULL,main=NULL)
abline(v=mean(resultMSE[10:rep,3]),col='red',lwd=3)
################################################################################


