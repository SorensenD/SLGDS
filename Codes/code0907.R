# CODE0907

# 1. GENERATES BINARY DATA y (AND MATRIX OF GENOTYPIC MARKERS X)
# 2. PERFORMS A GWAS ANALYSIS USING BONFERRONI CORRECTION
# 3. PERFORMS AN McMC ANLAYSIS WITH A SPIKE AND SLAB MODEL 
#    AND COMPUTES BAYESIAN FDR
# 4. COMPUTES CLASSICAL BENJAMINI AND HOCHBERG FDR-BH

rm(list=ls()) # CLEAR WORKSPACE
#install.packages("modeest")
library(modeest)

set.seed(30337) # FOR rq=1
#set.seed(303371) # FOR rq=0.4

##############################
# 1. INITIALISE AND CREATE DATA


## MARKER FREQUENCY:
pm <- 0.5
## QTL FREQUENCY:
pQTL <- 0.05
#mu <- qnorm(pm)
## FREQUENCY OF 1's IN POPULATION
pY <- 0.05
mu <- qnorm(pY)


va<-0.5 # additive variance of liability
ve<-1.0 # environmental variance
herprobit <- va/(va+ve)

nqtl<- 10
nindiv <- 2500 #2.500 in training and 2500 in validating
#nindiv <- 5000 #2.500 in training and 2500 in validating

########################################################
# SPECIFY THE PROPORTION OF va EXPLAINED BY THE nqtl LOCI
#rq <- 0.4
#rq <- 0.5
rq <- 1
########################################################

# Choose number of markers
nmark<-15000
#nmark <- 1500
# Choose length of McMC chain
#rep <- 5000
rep<-100 # SET NUMBER OF REPLICATES TO A SMALL NUMNER TO TEST THE CODE
#rep <- 10000


# Set initial proportion of QTL
piqtl<-nqtl/nmark
# parameters a and b of the beta prior
sc1=1.5
sc2=10
IDq<-sample(1:nmark,nqtl,replace=F) # from the nmark markers, choose nqtl as QTL
X<-matrix(data=NA,nrow= nindiv,ncol= nmark) # Initialise matricx of genetic markers
#### CREATE MATRIX OF GENETIC MARKERS X[,-IDq] AND OF THE nqtl QTL LOCI X[,IDq]
X[,-IDq] <-matrix(nrow= nindiv,ncol= (nmark-nqtl),rbinom(nindiv*(nmark-nqtl),
                                                         size=2,p=pm))
X[,IDq] <-matrix(nrow= nindiv,ncol= nqtl,rbinom(nindiv*nqtl,size=2,p=pQTL))
### INITIALISE MATRICES AND VECTORS 
be<-matrix(data=0.0,nrow=nmark,ncol=1) # parameter from true model
b<-matrix(data=0.0,nrow=nmark,ncol=1) # parameter from instrumental model
sumb<-matrix(data=0.0,nrow=nmark,ncol=1) # auxiliary vector used to calculate 
#  McMC E(b|y_t)
avrb <- matrix(data=0.0,nrow=nmark,ncol=1) # MCMC E(b|y)

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
result<-matrix(data=NA,nrow=rep,ncol=10)
resultBFDR<-matrix(data=NA,nrow=rep,ncol=nmark)
resultprobtheta<-matrix(data=NA,nrow=rep,ncol=nmark)
mse<-matrix(data=NA,nrow=rep,ncol=6)

keepypredval <- matrix(data=NA,nrow=rep, ncol=nindiv)
keepypredtrain <- matrix(data=NA,nrow=rep, ncol=nindiv)
keepystarval <- matrix(data=NA,nrow=rep, ncol=nindiv)
keepystartrain <- matrix(data=NA,nrow=rep, ncol=nindiv)
keepystarvalbhat <- matrix(data=NA,nrow=rep, ncol=nindiv)
storeb <-matrix(data=NA, nrow=rep,ncol=nmark)
storevarxb <- matrix(data=NA, nrow=rep,ncol=1)
predsel <- matrix(data=NA,nrow=rep,ncol=nindiv)
varxbsel <- matrix(data=NA,nrow=rep,ncol=1)
storeprobval <- matrix(data=NA, nrow=rep,ncol=nindiv)
probvalmean <- rep(0,nindiv)

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
msebriersel <- rep(0,rep)
msebrierselhigh <- rep(0,rep)
msevnrtheta <- rep(0,rep)
ypred<-matrix(data=NA,nrow=rep,ncol=nindiv)
### END INITIALISATION
cmX<-colMeans(X)
## CENTER MATRIX X
for (i in 1:nmark)
{Xc[,i]<- X[,i]-cmX[i]
}
#qr(Xc)$rank
meanXc <- apply(Xc,2,mean)
varXc<-apply(Xc,2,var)
sumvar<-sum(varXc[IDq]) # Compute 2*Sum(p(1-p)), where the Sum is over the nqtl 
sumvart <- sum(varXc)
##   QTL
QTLeff<-sqrt(rq*va/sumvar) # calculate the QTL effect so that the total genetic 
##   variance is VA
be[IDq] <- QTLeff
Xt<-Xc[,IDq]
#Sb<- (va/sumvar)*(6.1/4.1)
Sb <- (va/(sumvart*(nqtl/nmark)))*(6.1/4.1)
# Initialise scale parameter of scaled inverse chi-square distributions for 
##   conditional variance V
#Sv<- 1*(6.1/4.1)
#############################################
## RESIDUAL VARIANCE OF THE LIABILITY:
resliab <- 1 + (1-rq)*va
##############################################
liab<-Xc%*%be+rnorm(nindiv,mean=0,sd=sqrt(resliab)) # liability PROBIT 
#  (excluding mu)

liab<-liab+mu # include mu
xb<-Xc%*%be
var(xb) # GENETIC/GENOMIC VARIANCE AT THE LEVEL OF LIABILITY
p1 <-pnorm(mu+xb,mean=0,sd=sqrt(resliab)) # PROBIT MODEL Pr(Y=1|mu,x,b)
#####################################
## Prob(p1 > 0.5) ####
length(which(p1>0.5))/length(p1)
######################################
y.train <- rbinom(nindiv,1,p1) # GENERATE TRAINING BINARY RECORDS
y.val <- rbinom(nindiv,1,p1) # GENERATE VALIDATING BINARY RECORDS
######################################################
mean(y.train) # PROPORTION OF 1'S TRAIN
mean(y.val) # PROPORTION OF 1'S VAL
sum(varXc[IDq]*QTLeff^2) # GENETIC VARIANCE
##########################################################################
######## PROBABILITIES Pr(y=1|X,QTLeff) FOR THREE GENOTYPES AT 1 LOCUS ###
pnorm(mu+2*QTLeff)
pnorm(mu+QTLeff)
pnorm(mu)
############# END OF GENERATION OF DATA ####################################
############################################################################
# 2. GWAS USING A LOGISTIC REGRESSION WITH BONFERRONI CORRECTION (TRAIN DATA)
GWAS=matrix(nrow=ncol(Xc),ncol=4)
colnames(GWAS)=c('estimate','SE','t-value','p-value')
ptm <- proc.time()
for(i in 1:ncol(Xc)){
  fm=glm(y.train~Xc[,i],family = binomial(link = "probit"))
  GWAS[i,]= summary(fm)$coef[2,]
}
proc.time()-ptm
#############################################################################
GWASdetctBon<-which(-log10(GWAS[,4]) > -log10(0.05/nmark))
length(GWASdetctBon)  # DISCOVERY SET BONFERRONI
### COMPUTE PREDICTIONS BASED ON GWAS/GLM
fm2 <- glm(y.train~Xc[,GWASdetctBon],family = binomial(link = "probit"))
coef(fm2)
pred2b <- predict(fm2,type="response")
summary(pred2b)

ind <- which(pred2b==max(pred2b))
ind
### PREDICTION FOR y.val[ind] 
Xc[ind,GWASdetctBon]
genscore <- t(as.matrix(Xc[ind,GWASdetctBon]))%*%
  coef(fm2)[2:(length(GWASdetctBon)+1)]
pind<- pnorm(coef(fm2)[1]+genscore)
pind
pred2b[ind]
y.val[ind]
mean(y.val)
pred2b[ind]/mean(y.val)
###################################################################

plot(-log10(GWAS[,4]),type='o',ylab='-log10-pValue',xlab='Marker ID',
     cex=0.5,col=4,cex.lab=1.3)
abline(h=-log10(0.05/nmark),lty=2,col=2,lwd=2)
points(IDq,-log10(GWAS[IDq,4]),pch="*",cex=2,col="orange")

cat('Bonferoni',-log10(0.05/nmark),'\n')
GWASdetctBon<-which(-log10(GWAS[,4]) > -log10(0.05/nmark))
GWASdetctUnc<-which(-log10(GWAS[,4]) > -log10(0.05))
discovBonf <- length(GWASdetctBon) # DISCOVERY SET BONFERRONI
discovUnc <- length(GWASdetctUnc) # DISCOVERY SET UNCORRECTED
sort(GWASdetctBon)
sort(IDq) # LABELS FOR QTL MARKER LOCI


###########################################################################
true<-IDq # LABELS FOR QTL MARKERS
idmark<-seq(1:nmark)
fals<-idmark[-IDq]
discov<-GWASdetctBon
length(discov)
accepth0<-which(-log10(GWAS[,4]) < -log10(0.05/nmark))

tdisc<-length(intersect(true,discov))
fnegt<-length(setdiff(true,discov))
fdisc<-length(setdiff(discov,true))
trneg<-length(intersect(fals,accepth0))

fdisc # NUMBER OF FALSE DICOVERIES (BONFERRONI)
##############################################################################


################################################################################
# 3. CODE FOR THE GIBBS SAMPLER
# USES DATA y AND MARKER MATRIX Xc AS INPUT

V <- 1
# Initialise Vg, genomic variance, which initially is set equal to the 
#  additive genetic variance
Vg<-va
sumvarqtl<-sum(apply(as.matrix(Xc[,IDq]),2,var))
# start Gibbs chain
#marker effect variance is initially set initially equal to genomic variance divided by 
# the sum of the variance of the columns of Xc corresponding to QTLs
Vb<-Vg/sumvarqtl
maxlogods<- 0
# **************************************
nindivt <- nindiv

genomicvaluetrain<- rep(0,nindivt)

# Compute col vector xx with sum of squared terms involved in t(Xc)%*%Xc for the
#  training data
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
  (pnorm(0,mean=meanliab,sd=sd)*(1-y.train) + 
     (1-pnorm(0,mean=meanliab,sd=sd))*y.train)
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
  varxb <- var(Xc%*%b)
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
  genomicvaluetest<- Xc%*%b # Pred genomic values for the current iteration - 
  #  testing data
  # genomicvaluetrain<- Xc%*%b # Pred genomic values for the current iteration - 
  #  training data
  probval <- pnorm(u+genomicvaluetest) # vector of Prob(y=1|x,b) for 
  #  validating/testing data
  #probtrain <- pnorm(u+genomicvaluetrain) # vector of Prob(y=1|x,b) 
  #  for training data
  probtrain <- probval
  avrerrorvar <- mean(probtrain*(1-probtrain)) # Avr Var(y|b,x) - training data
  mse2briertest[i]<-mean((probval-y.val)^2) # MSE_v2 using validating data 
  #  based on probabilities
  mse2briertrain[i]<-mean((probval-y.train)^2) # MSE_v2 using validating data 
  #  based on probabilities
  
  y_predval <- as.numeric(ifelse(probval > 0.5, 1, 0))

  

  y_predtrain <- y_predval
  y_starval <- rbinom(length(probval),size=1,p=probval)
  y_startrain <- rbinom(length(probtrain),size=1,p=probtrain)
  mse2test[i] <- mean((y_predval - y.val) ^ 2) # TEST MSE using testing data 
  #  based on y_predval
  mse2train[i] <- mean((y_predtrain - y.train) ^ 2) # TTRAINING MSE using 
  #  training data based on y_predtrain
  mse3test[i] <- mean((y_starval - y.val) ^ 2) # TEST MSE using testing data 
  #  based on y_starval
  mse3train[i] <- mean((y_startrain - y.train) ^ 2) # TEST MSE using testing 
  #  data based on y_startrain

  resultprobtheta[i,] <- theta1
  result[i,]<- c(i,u,piqtl,nqtl/nmark,Vb,V,delta[IDq[1:3]],avrerrorvar)
  mse[i,] <- c(mse2train[i],mse2test[i],mse3train[i],mse3test[i],
               mse2briertest[i],mse2briertrain[i])
##############################################################################
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
  sumb <- sumb + b # SUM OF MARKER EFFECTS
  #############################################################################
  
  # ###########################################################################
  storeb[i,] <- b # STORE MARKER EFFECTS
  storevarxb[i,] <- varxb
  storeprobval[i,] <- probval # Store Prob(y=1|x,b) for validating/testing data

}
proc.time()-ptm
##############################################################################
postprob <- postprob/rep

plot(postprob,type='o',ylab='PostProb',xlab='Marker ID',cex=.5,col=4,cex.lab=1.3)
abline(h=0.8,lty=2,col=2,lwd=2)
points(IDq,postprob[IDq],pch="*",cex=2,col="orange")

selM <- which(postprob>0.4) # SELECTED MARKERS BASED ON PROBS > 0.4
colmeanprob <- apply(storeprobval[1:rep,],2,mean)
varprobval <- apply(storeprobval[1:rep,],2,var)
summary(colmeanprob)
summary(varprobval)
tail(sort(colmeanprob))
###### MEAN AND VARIANCE OF HIGHEST AND LOWEST PROBABILITES - ALL MARKERS
probyeq1max <- which(colmeanprob==max(colmeanprob))
colmeanprob[probyeq1max]
probyeq1min <- which(colmeanprob==min(colmeanprob))
colmeanprob[probyeq1min]
var(storeprobval[,which(colmeanprob==max(colmeanprob))])
var(storeprobval[,which(colmeanprob==min(colmeanprob))])
##########################################################################
##############################################################################
# COMPUTE ESTIMATED JOINT COUNTS
# TO OBTAIN ESTIMATED JOINT PROBABILTIES DIVIDE BY NINDIV=2500
# Pr(Y*>2*null,Yv=1)
p11 <- length(which(colmeanprob > (2*mean(y.val)) & y.val == 1)) 
# Pr(Y*>2*null,Yv=0)
p10 <- length(which(colmeanprob > (2*mean(y.val)) & y.val == 0))
# Pr(Y*<2*null,Yv=1)
p01 <- length(which(colmeanprob < (2*mean(y.val)) & y.val == 1))
# Pr(Y*<2*null,Yv=0)
p00 <- length(which(colmeanprob < (2*mean(y.val)) & y.val == 0)) 
# COMPUTE ESTIMATED CONDITIONAL PROBABILITIES
p11/length(which(y.val==1)) # Pr(Y*>2*null|Yv=1)
p10/length(which(y.val==0)) # Pr(y*>2*null|Yv=0)
p01/length(which(y.val==1)) # Pr(Y*<2*null|Yv==1)
p00/length(which(y.val==0)) # Pr(Y*<2*null|Y=0)
################################################################################
################################################################################
###########   SELECTED MARKERS BASED ON POSTPROB > 0.4    ######################
selM <- which(postprob>0.4) # LABLE FOR HIGHEST SCORING MARKERS (FOR PROB > 0.4)
## CONSTRUCT PREDICTORS BASED ON GENOTYPES OF HIGHEST SCORING MARKERS
for (i in 1:rep){
  predsel[i,] <- pnorm(result[i,2] + Xc[,selM]%*%storeb[i,selM])
}
# which(predsel == max(predsel))
# BRIER SCORE USING SELECTED MARKERS, ACOUNTING FOR VARIANCE OF THE PREDICTOR
for (i in 1:rep){
  msebriersel[i] <- mean((y.val - predsel[i,])^2)
}
msepredselBrier <- mean(msebriersel) # BRIER SCORE SELECTED MARKERS
# BRIER SCORE USING ALL MARKERS, ACOUNTING FOR VARIANCE OF THE PREDICTOR
mseprobvalBrier <- mean(mse[,5])
## MEAN AND VARIANCE OF POSTERIOR PREDICTIVE PREDICTIONS BASED ON Pr(Y*=1|X,Yt)
## USING SELECTED MARKERS:
meanpredsel <- apply(predsel,2,mean)
varpredsel <- apply(predsel,2,var)

summary(varpredsel)  # VAR OF Pr(Y*=1|X,Yt) USING SEL MARKERS
summary(meanpredsel)
tail(sort(meanpredsel))
summary(varprobval) # VAR OF Pr(Y*=1|X,Yt) USING ALL MARKERS
############################################################################

### THEORETICAL EXPECTATIONS OF THE BRIER SCORES ACCOUNTING FOR VAR OF PREDICTOR
#1. SELECTED MARKERS:
as <- mean(varpredsel)
bs <- mean((y.val-meanpredsel)^2)
expmseSEL <- as + bs #  theoretical expectation
expmseSEL
msepredselBrier # McMC estimate of theoretical expectation
#2. ALL MARKERS
aa <- mean(varprobval)
ba <- mean((y.val-colmeanprob)^2)
expmseALL <- aa + ba #  theoretical expectation
expmseALL
mean(mse[,5]) # McMC estimate of theoretical expectation
#############################################################################
############ GENOMIC VARIANCE USING ALL MARKERS
summary(storevarxb[50:rep])
quantile(storevarxb[50:rep],c(0.025,0.975))
summary(storevarxb[50:rep]/(1+storevarxb[50:rep]))
quantile(storevarxb[50:rep]/(1+storevarxb[50:rep]),c(0.025,0.975))
###############   BRIER SCORES SELECTED MARKERS - ALL DATA ###############
mean((y.val-meanpredsel)^2) # BRIER score using selected markers
mean((y.val-colmeanprob)^2) # BRIER score using all markers
mean((y.val-mean(y.val))^2) # BRIER score of null model
mean((y.val-pnorm(mean(result[,2])))^2) # BRIER SCORE NULL MODEL 
###################

### GENOMIC VARIANCE USING SELECTED MARKERS ##########################
avrb <- sumb/rep # MCMC E(b|y_t): Xb_hat
summary(avrb)
var(Xc[,selM]%*%avrb[selM])
for(i in 1:rep){
  varxbsel[i] <- var(Xc[,selM]%*%storeb[i,selM])
}
## SUMMARY GENOMIC VARIANCE SEL MARKERS:
summary(varxbsel[50:rep])
quantile(varxbsel[50:rep],c(0.025,0.975))
summary(varxbsel[50:rep]/(1+varxbsel[50:rep]))
quantile(varxbsel[50:rep]/(1+varxbsel[50:rep]),c(0.025,0.975))
##################################################################
## MEAN AND VARIANCE OF POSTERIOR PREDICTIVE PREDICTIONS BASED ON Pr(Y*=1|X,Yt)
## USING SELECTED MARKERS:
meanpredsel <- apply(predsel,2,mean)
varpredsel <- apply(predsel,2,var)
summary(varpredsel)  # VAR OF Pr(Y*=1|X,Yt) USING SEL MARKERS
summary(meanpredsel)
tail(sort(meanpredsel))
summary(varprobval) # VAR OF Pr(Y*=1|X,Yt) USING ALL MARKERS

########################################################
#### COMPUTE MEAN AND 95% CI FOR MAX AND MIN Pr(Y*=1|X,Yt) SELECTED MARKERS
probyeq1selmax <- which(meanpredsel==max(meanpredsel))
probyeq1selmin <- which(meanpredsel==min(meanpredsel))
length(probyeq1selmin)
minp <- probyeq1selmin[1]
meanpredsel[probyeq1selmax]
meanpredsel[probyeq1selmin[1]]
var(predsel[,which(meanpredsel==max(meanpredsel))])
var(predsel[,minp])
sd(predsel[,minp])
quantile(predsel[,probyeq1selmax],c(0.025,0.975))
meanpredsel[probyeq1selmin][1]
quantile(predsel[,probyeq1selmax],c(0.025,0.975))
quantile(predsel[,probyeq1selmin],c(0.025,0.975))
X[4,selM]
###########################################################################
####### LARGEST/SMALLEST 100 PROBABILITIES BASED ON SELECTED MARKERS
# 1. LARGST
sortpredsel <- sort(meanpredsel)
## LARGEST 100
large100 <- sortpredsel[2401:2500]
mean(large100)
min(large100)
idlarge100 <- which(meanpredsel >= min(large100))
mean(y.val[idlarge100])
mean(y.val)
# 2. SMALLEST
sortpredsel <- sort(meanpredsel)
## SMALLEST 100
smol100 <- sortpredsel[1:100]
mean(smol100)
min(smol100)
idsmol100 <- which(meanpredsel <= max(smol100))
mean(y.val[idsmol100])
mean(y.val)
###########################################################################
##############################################################################
################  logscores  #####################################
probnullmod <- pnorm(mean(result[50:rep,2]))
logliknull <- sum(y.val*log(probnullmod)+(1-y.val)*log((1-probnullmod)))
loglikall <- sum(y.val*log(colmeanprob)+(1-y.val)*log((1-colmeanprob)))
logliksel <- sum(y.val*log(meanpredsel)+(1-y.val)*log((1-meanpredsel)))
##############################################################################
#########################################################################
#### AS GOOD AS IT GETS #######################################################
genomicvaluetrue <- Xc%*%be 
probtrue <- pnorm(mu+genomicvaluetrue)
summary(probtrue)
probtrue[probyeq1max]
pnorm(mu + Xc[probyeq1max,IDq]%*%be[IDq])
y_predtrue <- as.numeric(ifelse(probtrue > 0.5, 1, 0))
mean((y_predtrue-y.val)^2)
##############################################################################

##############################################################################
# ******************************************
# COMPUTE BAYESIAN FALSE DISCOVERY RATE
#install.packages("dplyr", .libPaths()[1])
#library(dplyr)
ordpp<-order(-postprob)
sortpp<-postprob[ordpp]
localfdr <- 1-sortpp
#ordpp<-order(-postprobability)
#sortpp<-postprobability[ordpp]

for (i in 1:nmark){
  avr[i]<-mean(sortpp[1:i])
  fd[i]<-1-avr[i]# THIS COMPUTES THE (CUMULATIVE) EXPECTED PROPORTION OF BFDR 
  #                  FOR EACH DISCOVERY SET OF SIZE 1 TO NMARK
  fdloc[i] <- mean(localfdr[1:i]) # MEAN OF LOCAL FDR IN DISCOVERY SET OF SIZE i
  # THE SAME AS fd[]
  # BELOW: PROPORTION OF TRUE FALSE DISCOVERIES IN THE DISCOVERY SET OF SIZE i
  truefd[i]<-length(setdiff(ordpp[1:i],IDq))/i 
  
}

####################################################################



##############################################################################
# FUNCTION TO GENERATE HISTOGRAM OF BAYESIAN FDR USING GIBBS OUTPUT
# INPUT:
# 1. resultprobtheta: A FILE WITH DRAWS FROM (7.57)
# 2. fd: THE TRUE PROPORTION OF FALSE DISCOVERIES OBSERVED IN THE SAMPLE
# 3. prob: THE USER-USED THRESHOLD THAT DEFINES THE DISCOVERY SET
# 4. postprob: VECTOR OF POSTERIOR PROBABILITIES 
#   (OF BELONGING TO THE SLAB COMPONENT OF THE MIXURE) FOR EACH GENETIC MARKER
fdrhist <- function(resultprobtheta,truefd,prob){
  discset <- which(postprob > prob)
  fdisc <- apply(1-(resultprobtheta[,discset]),1,mean)
  # MEAN AND POSTERIOR INTERVAL FOR BAYES FDR:  
  avfdis <- mean(fdisc)
  quantilefdis <- quantile(fdisc,c(0.025,0.975))

 # HISTOGRAM OF McMC ESTIMATE OF POSTERIOR DISTRIBUTION OF BAYES FDR:  
  hist(fdisc,breaks=20,xlab='McMC-Bayes FDR',main=NULL, freq=FALSE,cex.lab=1.3)
  ### TRUE FDR: 
  abline(v=truefd[length(discset)],col="red",lwd=2) # TRUE FDR IN SAMPLE
 # dev.off()
  return <- c(avfdis,quantilefdis,length(discset),truefd[length(discset)]*
                length(discset))
}
out <- fdrhist(resultprobtheta,truefd,0.8)
out <- fdrhist(resultprobtheta,truefd,0.52)
out <- fdrhist(resultprobtheta,truefd,0.4)
out <- fdrhist(resultprobtheta[500:rep,],truefd,0.12)

out
##############################################################################
#############################################################################
#   4. FDR-BH (BenjaminiHochberg)
# USES OUTPUT FROM THE GWAS ANALYSIS AS INPUT (P-VALUES)
bID<-seq(1:nmark)
trueH1<-IDq
trueH0<-bID[-IDq]
# ALFA HERE IS THE q IN BOOK 
# (chosen value of expected proportion of false discoveries)
#alfa<-0.05
alfa<-0.15 
alfa <- 0.10
alfa <- 0.20
#alfa <- 0.5
adjpv<-vector()
qvalue<-vector()
qvbh<-vector()
pval<-GWAS[,4]
bhat<-GWAS[,1] 
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
sizediscovset
fdiscov<-length(setdiff(signif$bID,trueH1))
fdiscov # NO. TRUE FALSE DISCOVERIES
fdiscov/sizediscovset # TRUE OBSERVED PROPORTION OF FD IN SAMPLE 
# EXPECTED PROPORTION ALFA (CONCEPTUALLY AVERAGED OVER SAMPLES)
alfa
#################################################################


