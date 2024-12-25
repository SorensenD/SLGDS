# CODE0908

# FITS A SPIKE AND SLAB PROBIT MODEL AND BELOW 
# FITS A NORMAL APPROXIMATION TO THE BINARY DATA


rm(list=ls()) # CLEAR WORKSPACE
#install.packages("modeest")
library(modeest)

set.seed(30337)
##############################
# 1. INITIALISE AND CREATE DATA


## MARKER ALLELE FREQUENCY:
pm <- 0.5
## QTL ALLELE FREQUENCY:
#pQTL <- 0.05
pQTL <- 0.2
#pQTL <- 0.5
#mu <- qnorm(pm)
## FREQUENCY OF 1's IN POPULATION
pY <- 0.2
mu <- qnorm(pY)


va<-0.5 # additive variance of liability
ve<-1.0 # environmental variance
herprobit <- va/(va+ve)

nqtl<- 25
#nindiv <- 2500 #2.500 in training and 2500 in validating
#nindiv <- 5000 #2.500 in training and 2500 in validating
nindiv <- 1000
########################################################
# SPECIFY THE PROPORTION OF va EXPLAINED BY THE nqtl LOCI
#rq <- 0.4
rq <- 1
########################################################


nmark <- 5000
# Choose length of McMC chain
rep <- 400
# Set initial proportion of QTL
piqtl<-nqtl/nmark
# parameters a and b of the beta prior
sc1=1.5
sc2=10
IDq<-sample(1:nmark,nqtl,replace=F) # from the nmark markers, choose nqtl as QTL
X<-matrix(data=NA,nrow= nindiv,ncol= nmark)

X[,-IDq] <-matrix(nrow= nindiv,ncol= (nmark-nqtl),rbinom(nindiv*(nmark-nqtl),size=2,p=pm))
X[,IDq] <-matrix(nrow= nindiv,ncol= nqtl,rbinom(nindiv*nqtl,size=2,p=pQTL))
be<-matrix(data=0.0,nrow=nmark,ncol=1) # parameter from true model
b<-matrix(data=0.0,nrow=nmark,ncol=1) # parameter from fitted model
sumb<-matrix(data=0.0,nrow=nmark,ncol=1) # auxiliary vector used to calculate McMC E(b|y_t)
avrb <- matrix(data=0.0,nrow=nmark,ncol=1) # MCMC estimates E(b|y)

theta1<-matrix(data=0.0,nrow=nmark,ncol=1) # parameter from fitted model
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
result<-matrix(data=NA,nrow=rep,ncol=10)
resultBFDR<-matrix(data=NA,nrow=rep,ncol=nmark)
resultprobtheta<-matrix(data=NA,nrow=rep,ncol=nmark)
mse<-matrix(data=NA,nrow=rep,ncol=7)

keepypredval <- matrix(data=NA,nrow=rep, ncol=nindiv)
keepypredtrain <- matrix(data=NA,nrow=rep, ncol=nindiv)
keepystarval <- matrix(data=NA,nrow=rep, ncol=nindiv)
keepystartrain <- matrix(data=NA,nrow=rep, ncol=nindiv)
keepystarvalbhat <- matrix(data=NA,nrow=rep, ncol=nindiv)

logscoreall <- matrix(data=NA, nrow=rep,ncol=1)
logscorenull <- matrix(data=NA, nrow=rep,ncol=1)

storeprobval <- matrix(data=NA, nrow=rep,ncol=nindiv)
probvalmean <- rep(0,nindiv)
storeprobtrain <- matrix(data=NA, nrow=rep,nco=nindiv)
drawfrom24 <- matrix(data=NA,nrow=rep, ncol=1)
drawfrom27 <- matrix(data=NA,nrow=rep, ncol=1)
drawfrom28 <- matrix(data=NA,nrow=rep, ncol=1)

alfa<-matrix(data=0,nrow=nmark,ncol=1)
postprob<-rep(0,nmark)
postprobff<-rep(0,nmark)

mse2briertest<-rep(0,rep)
mse2briertrain <- rep(0,rep)
mse2test <- rep(0,rep)
mse22test <- rep(0,rep)

mse2train <- rep(0,rep)
mse3test <- rep(0,rep)
mse3train <- rep(0,rep)
msetestbhat <- rep(0,rep)

msevnrtheta <- rep(0,rep)
ypred<-matrix(data=NA,nrow=rep,ncol=nindiv)
cmX<-colMeans(X)
###########################################
#Center and scale matrix X
for (i in 1:nmark)
{Xc[,i]<- (X[,i]-cmX[i])/sd(X[,i])
}
############################################
#qr(Xc)$rank
meanXc <- apply(Xc,2,mean)
varXc<-apply(Xc,2,var)
sumvar<-sum(varXc[IDq]) # Compute 2*Sum(p(1-p)), where the Sum is over the nqtl QTL
QTLeff<-sqrt(rq*va/sumvar) # calculate the QTL effect so that the total genetic variance is VA
be[IDq] <- QTLeff
Xt<-Xc[,IDq]
Sb<- (va/sumvar)*(6.1/4.1)
# Initialise scale parameter of scaled inverse chi-square distributions for conditional variance V
Sv<- 1*(6.1/4.1)
#############################################
## RESIDUAL VARIANCE OF THE LIABILITY:
resliab <- 1 + (1-rq)*va
##############################################
liab<-Xc%*%be+rnorm(nindiv,mean=0,sd=sqrt(resliab)) # liability PROBIT (excluding mu)

liab<-liab+mu # include mu
xb<-Xc%*%be
p1 <-pnorm(mu+xb,mean=0,sd=sqrt(resliab)) # PROBIT MODEL
#####################################
## Prob(p1 > 0.5) ####
length(which(p1>0.5))/length(p1)
######################################
y.train <- rbinom(nindiv,1,p1) # TRAINING BINARY RECORDS
y.val <- rbinom(nindiv,1,p1) # VALIDATING BINARY RECORDS
z <- c(y.train,y.val)
#######################################################
y.val2 <- as.numeric(ifelse(p1 > 0.5, 1, 0))
######################################################
mean(y.train) # PROPORTION OF 1'S TRAIN
mean(y.val) # PROPORTION OF 1'S VAL
sum(varXc[IDq]*QTLeff^2)
var(Xc%*%be)

############# END OF GENERATION OF DATA ####################################
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
nindivt <- nindiv

# Compute col vector xx with sum of squared terms involved in t(Xc)%*%Xc for the training data
for (i in 1:nmark){xx[i]<-crossprod(Xc[,i])}
ptm<-proc.time()
# **************************************
u<-mean(liab)
# *****************************
res1<-liab-u
# ****************************
sumb <- 0
for (i in 1:rep)
{
  print(i)
  # sample liabilities liab from truncated normals
  meanliab <- u + Xc%*%b
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
    res0<-res1+Xc[,j]*b[j,1]
    rss0<-sum(res0*res0)
    if(delta[j,1] == 1) 
    {
      alfahat<- (t(Xc[,j])%*%(res0))/(xx[j]+k)
      varalfa<-(1/(xx[j]+k))*V
      alfa[j,1]<- rnorm(1,mean=alfahat,sd=sqrt(varalfa))
    } else {
      alfa[j,1]<-rnorm(1,mean=0,sd=sqrt(Vb)) 
    }
    # sample deltas
    res1temp<-res0-Xc[,j]*alfa[j,1]
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
    res1<-res0-as.matrix(Xc[,j])*b[j,1]
    postprobff[j] <- postprobff[j] + deltaff[j,1]
  }
  
  # Update Nqtl = sum(delta[,1])
  rss1<-sum(res1*res1)
  Nqtl<-sum(delta[,1])
  # sample piqtl
  piqtl<-rbeta(1,Nqtl+sc1,nmark-Nqtl+sc2)
  
  # update variance of a marker effect (Vb)
  df<- 4.1+nmark
  Vb<- (sum(alfa*alfa) + Sb*4.1)/rchisq(1,df)
    genomicvaluetest<- Xc%*%b # Pred genomic values for the current iteration - testing data
  probval <- pnorm(u+genomicvaluetest) # vector of Prob(y=1|x,b) for validating/testing data
  probtrain <- probval
  avrerrorvar <- mean(probtrain*(1-probtrain)) # Avr Var(y|b,x) - training data
  mse2briertest[i]<-mean((probval-y.val)^2) # MSE_v2 using validating data based on probabilities
  mse2briertrain[i]<-mean((probval-y.train)^2) # MSE_v2 using validating data based on probabilities
  y_predval <- as.numeric(ifelse(probval > 0.5, 1, 0))
  y_predtrain <- y_predval
  y_starval <- rbinom(length(probval),size=1,p=probval)
  y_startrain <- rbinom(length(probtrain),size=1,p=probtrain)
  mse2test[i] <- mean((y_predval - y.val) ^ 2) # TEST MSE using testing data based on y_predval
  mse2train[i] <- mean((y_predtrain - y.train) ^ 2) # TTRAINING MSE using training data based on y_predtrain
  mse3test[i] <- mean((y_starval - y.val) ^ 2) # TEST MSE using testing data based on y_starval
  mse3train[i] <- mean((y_startrain - y.train) ^ 2) # TEST MSE using testing data based on y_startrain
  mse22test[i] <- mean((y_predval - y.val2) ^ 2)
  resultprobtheta[i,] <- theta1
  result[i,]<- c(i,piqtl,Nqtl/nmark,Vb,V,delta[IDq[1:3]],u,avrerrorvar)
  mse[i,] <- c(mse2train[i],mse2test[i],mse3train[i],mse3test[i],mse2briertest[i],mse2briertrain[i],
               mse22test[i])
###########################################################################################
  ##############################################
  # KEEP y_predval AND y_predtrain
  keepypredval[i,] <- y_predval
  keepypredtrain[i,] <- y_predtrain
  # ################
  ##############################################
  # KEEP y_starval AND y_startrain
  keepystarval[i,] <- y_starval
  keepystartrain[i,] <- y_startrain
    #####################################################
  sumb <- sumb + b
  ######################################################################################
  storeprobval[i,] <- probval # Store Prob(y=1|x,b) for validating/testing data
  storeprobtrain[i,] <- probtrain # Store Prob(y=1|x,b) for training data
  logscoreall[i,] <- sum(y.val*log(probval)+(1-y.val)*log((1-probval)))
  logscorenull[i,] <- sum(y.val*log(pnorm(u))+(1-y.val)*log((1-pnorm(u))))
 }
proc.time()-ptm
##############################################################################
colmeanprob <- apply(storeprobval[100:rep,],2,mean)
summary(colmeanprob)
##########################################################################
##############  PROPORTION OF 1S AMONG 50 HIGHEST AND 50 LOWEST PROBABILITIES
sortprob <- sort(colmeanprob)
pmax <- sortprob[951:1000]
idmax <- which(colmeanprob > 0.5983372)
length(idmax)
sum(y.val[idmax])
# 77
pmin <- sortprob[1:50]
idmin <- which(colmeanprob < 0.08433434)
length(idmin)
sum(y.val[idmin])
# 1
################  logscores  #####################################
probnullmod <- mean(y.val)
logliknull <- -sum(y.val*log(probnullmod)+(1-y.val)*log((1-probnullmod)))
loglikall <- -sum(y.val*log(colmeanprob)+(1-y.val)*log((1-colmeanprob)))
logliknull
loglikall
############################# BRIER SCORE ########################
mean(mse[100:rep,5])
quantile(mse[,5],c(0.025,0.975))
############################################################
###########################  PROP OF MISCLASSIFICATIONS (MSE2) #####
mean(mse[100:rep,2])
quantile(mse[,2],c(0.025,0.975))
#################################################################################
###############################################################################
###############   FIT THE NORMAL APPROXIMATION   ##############################
##################################################################################
### IN ORDER TO USE THE SAME DATA AS USED WITH THE PROBIT MODEL:
##########   EXECUTE THE NORMAL APPROXIMATION FROM HERE #####################
########  AFTER EXECUTING OptimismBinary.r (THE PROBIT MODEL) ################
keepmualpha <- matrix(data=NA,nrow=rep, ncol=nindiv)
keepyvpred <- matrix(data=NA,nrow=rep, ncol=nindiv)
resultMSEbrier <- matrix(data=NA,nrow=rep,ncol=2)
resultMSE3 <- matrix(data=NA,nrow=rep,ncol=2)
meff <- matrix(data=NA,nrow=nindiv,ncol=1)
resultridge<-matrix(data=NA,nrow=rep,ncol=7)

storeyvpred <- matrix(data=NA, nrow=rep,ncol=nindiv)
storeystarval <- matrix(data=NA, nrow=rep,ncol=nindiv)
storeyvprob <- matrix(data=NA, nrow=rep,ncol=nindiv)
storeyvprobnull <- matrix(data=NA, nrow=rep,ncol=1)
W <- Xc
#########################################################################
## HYPERPARAMETERS OF THE SCALED INVERTED PRIORS FOR va AND ve
nua <- 4.1
nue <- 4.1
Sa <- 10
Se <- 30
## THIS GENERATES PRIOR MODES FOR va AND ve equal to 10*(4.1/6.1) = 6.7
## and 10*(4.1/6.1) = 20.2
#######################################################################
y <- y.train
yv <- y.val
mean(y)
var(y)
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
msev3 <- rep(0,rep)
#rep <- 40000
#INITIALISE resultrudge
resultridge<-matrix(data=NA,nrow=rep,ncol=7)
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
  alfa1<-rnorm((nindiv-1),meanalfa1,sqrt(varalfa1))
  alfa<-c(alfa1,0)
  # SAMPLE Vg
  # COMPUTE SCALE
  #  scVg<-sum(alfa1*alfa1*(1/Dp))
  scVg<-sum(alfa1*alfa1*(1/Dp)) + nua*Sa
  #  Vg<-scVg/rchisq(1,nindiv-3)
  Vg<-scVg/rchisq(1,(nindiv-1+nua))
  
  #Vg<-0.0001
  # SAMPLE Ve
  # COMPUTE SCALE
  Ualfa <- U%*%alfa
  ystar<-y-mu-Ualfa
  mualfa <- mu + Ualfa
  #  scVe<-sum(ystar*ystar)
  
  scVe<-sum(ystar*ystar) + nue*Se
  #  Ve<-scVe/rchisq(1,nindiv-2)
  Ve<-scVe/rchisq(1,(nindiv + nue))
  
  #Ve<-25
  k<-Ve/Vg
  uund <- qnorm(mu)
  yvprobnull <- mu
  yvprob <- pnorm(uund + Ualfa/dnorm(uund))
  y_starval <- rbinom(length(yvprob),size=1,p=yvprob)
  
  ypredval <- as.numeric(ifelse(yvprob > 0.5, 1, 0))
  msevbrier <- mean((yv-yvprob)^2)
  msetbrier <- mean((y-yvprob)^2)
  msevyval <- mean((yv-ypredval)^2)
  msetyval <- mean((y-ypredval)^2)
  resultridge[i,]<-c(i,mu,Vg,Ve,Vg/(Vg+Ve),1/k,mean(alfa*alfa))
  resultMSEbrier[i,] <- c(msevbrier,msetbrier)
  resultMSE3[i,] <- c(msevyval,msetyval)

  storeyvpred[i,] <- ypredval
  storeyvprob[i,] <- yvprob
  storeyvprobnull[i,] <- yvprobnull
  storeystarval[i,] <- y_starval
}
proc.time()-ptm
################################
##############################################################################
colmeanprob_na <- apply(storeyvprob[1:rep,],2,mean)
summary(colmeanprob_na)
##########################################################################
##############  PROPORTION OF 1S AMONG 50 HIGHEST AND 50 LOWEST PROBABILITIES
sortprobna <- sort(colmeanprob_na)
pmaxna <- sortprobna[951:1000]
idmaxna <- which(colmeanprob_na > 0.5816794)
length(idmaxna)
sum(y.val[idmaxna])
# 77
pminna <- sortprobna[1:50]
idminna <- which(colmeanprob_na < 0.2079941)
length(idminna)
sum(y.val[idminna])
# 1
################  logscores  #####################################
probnullmod <- mean(yv)

logliknull <- -sum(y.val*log(probnullmod)+(1-y.val)*log((1-probnullmod)))
logliknull
loglikall <- -sum(yv*log(colmeanprob_na)+(1-yv)*log((1-colmeanprob_na)))
loglikall
#############################################################################
############################# BRIER SCORE ########################
mean(resultMSEbrier[,1])
quantile(resultMSEbrier[,1],c(0.025,0.975))
############################################################
###########################  PROP OF MISCLASSIFICATIONS (MSE2) #####
mean(resultMSE3[,1])
quantile(resultMSE3[,1],c(0.025,0.975))
###################################################################

