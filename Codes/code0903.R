# CODE0903
rm(list=ls()) # CLEAR WORKSPACE
library(glmnet)
set.seed(30337)
va<-0.3
#va <- 1
p<-0.25
#p <- 0.5
mu <- qnorm(p)
ve <- 1
nindiv<-2000
nloci<-20
nmarker<-500
be<-matrix(data=0.0,nrow=nmarker,ncol=1) # parameter true model
IDq<-sample(1:nmarker,nloci,replace=F) # choose nloci as QTL
WT<-matrix(nrow=nindiv,ncol= nmarker,
           rbinom(n=nindiv*nmarker,size=2,p=.5))
XT<-matrix(data=NA,nrow=nindiv,ncol=nmarker) # NO INTERCEPT
XTi<-matrix(data=NA,nrow=nindiv,ncol=nmarker+1) # INTERCEPT
# CENTER MARKER MATRIX
cm<-colMeans(WT)
for (i in 1:nmarker)
{
  XT[,i]<-(WT[,i]-cm[i])
}
meanXT <- apply(XT,2,mean)
varXT<-apply(XT,2,var)
# Compute 2*Sum(p(1-p)); Sum over nqtl QTL:
sumvar<-sum(varXT[IDq]) 
QTLeff<-sqrt(va/sumvar) # QTL effect 
#    so that genetic variance is VA
be[IDq] <- QTLeff # QTL EFFECT
xb<-XT%*%be
p1 <- pnorm(mu+xb) # PROBIT MODEL
be[IDq]<-QTLeff # TRUE MARKER EFFECT = QTLeff; REST ARE ZERO
var(XT%*%be)
y<-rep(0,nindiv)
y <- rbinom(nindiv,1,p1)
mean(y)
################################################################
# FITTING FULL MODEL AND LASSO WITH GLMNET
set.seed(3337)
#library(glmnet)

train=sample(1:nrow(XTi),nrow(XTi)/2)
test=(-train)
y.test=y[test]
y.train<-y[train]
# ********** FIT GLMNET TO FIND THE BEST LAMBDA ***************
# STEP 1
cv.out=cv.glmnet(
  XT[train,],y[train],alpha=1,standardize=TRUE,family="binomial")
plot(cv.out)
bestlam=cv.out$lambda.min
length(which(as.vector(coef(cv.out,s=bestlam))!=0))

# STEP 2

fm.pred0=predict(
  cv.out,s=0,newx=XT[test,],family="binomial",type="class")
fm.pred=predict(
  cv.out,s=bestlam,newx=XT[test,],family="binomial",type="class")


mean((as.numeric(fm.pred)-y.test)^2) # VALIDATION MSE: LASSO

mean((as.numeric(fm.pred0)-y.test)^2) # VALIDATION MSE: FULL MODEL

# ERROR RATE OF NULL MODEL: y = mu + e
pnull<-mean(y.test)
pnull
ynull<-rep(0,length(y.test))
if(pnull > 0.5){ynull<-1}
mean((ynull-y.test)^2) # VALIDATION MSE BASED ON NULL MODEL
########################################################################
# GENERATING ROC CURVES
#install.packages("pROC", .libPaths()[1])
library(pROC)
set.seed(420)
par(pty="s") # GENERATES "PRETTY" FIGURES

fm.pred=predict(
  cv.out,s=bestlam,newx=XT[test,],family="binomial",type="response")
fm.pred0=predict(
  cv.out,s=0,newx=XT[test,],family="binomial",type="response")
pred<-as.numeric(fm.pred)
pred0<-as.numeric(fm.pred0)
roc.info<-roc(y.test, pred,plot=FALSE,legacy.axes=TRUE) # LASSO
# FULL MODEL:
roc.info0<-roc(y.test, pred0,plot=FALSE,legacy.axes=TRUE) 


s <-roc(y.test, pred,plot=FALSE,legacy.axes=TRUE,percent=TRUE,
        xlab="False Positive Percentage",ylab="True Positive Percentage",
        col="blue",print.auc=TRUE,cex.lab=1.3)

# CREATA DATA FRAME TO BE USED IN THE NEXT CODE
# THAT CONTAINS TPOS, FPOS FOR THE COMPLETE RANGE
# OF THRESHOLDS t
roc.df<-data.frame(tpos=roc.info$sensitivities*100,
                   fpos=(1-roc.info$specificities)*100,
                   thresholds=roc.info$thresholds)
roc.df0<-data.frame(tpos0=roc.info0$sensitivities*100,
                    fpos0=(1-roc.info0$specificities)*100,
                    thresholds0=roc.info0$thresholds)
head(roc.df)
f50<-roc.df[min(which(roc.df$thresholds > 0.49995)),]
f50_0<-roc.df0[min(which(roc.df0$thresholds > 0.49995)),]
par(pty="m")
f50
mean(y.test)
#########################################################################
# COMPUTE AUC BY NUMERICAL INTEGRATION
tpos <- roc.df$tpos/100
fpos <- roc.df$fpos/100
n<-length(tpos)
sum((tpos[-1]+tpos[-n])/2*(fpos[-n]-fpos[-1]))

# OR AS A FUNCTION
auc <- function(fpos,tpos){
  n <- length(tpos)
  abs(sum((tpos[-1]+tpos[-n])/2*(fpos[-n]-fpos[-1])))
}
with(roc.df,auc(fpos/100,tpos/100)) # AUC FOR LASSO
with(roc.df0,auc(
  fpos0/100,tpos0/100)) # AUC FULL MODEL: glmnet: LAMBDA=0
