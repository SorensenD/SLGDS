# CODE1305
# LIKELIHOOD PROBLEMS II QUESTION 4
# DATA BASED ON GENOMIC MODEL AND OBTAIN THE SVD OF WW'(1/m)
rm(list=ls()) # CLEAR WORKSPACE
set.seed(1327)
nindiv<-2000
nmark<-20000
nt<-nindiv*nmark
X<-matrix(nrow= nindiv,ncol= nmark,rbinom(n=nt,size=2,p=.5))
stdev <- matrix(data=NA,nrow= nmark,ncol=1)
W <- matrix(data=NA,nrow= nindiv,ncol=nmark)
U <- matrix(data=NA,nrow= nindiv,ncol= nindiv)
G<-matrix(data=NA,nrow= nindiv,ncol= nindiv)
cm <- colMeans(X)
# CREATE MATRIX OF STANDARDISED MARKER GENOTYPE CODES
for (i in 1:nmark)
{
  W[,i] <-( X[,i]-cm[i]) / sd(X[,i])
}
# COULD USE INSTEAD:
# W <- scale(X)
#qr(X)$rank
#qr(W)$rank
#nmark MARKER VALUES: REALISATIONS FROM N(0,I sqrt(0.001))
g <- (1/sqrt(nmark))*W%*%rnorm(nmark,mean=0,sd=sqrt(10))
# GENERATE nindiv PHENOTYPES WITH MEAN 0, VAR=10+15, 
# GENOMIC HERITABILITY=10/(10+15)=0.4
#PARAMETER k = Vg/Ve = 10/15 =0.67
y <- g+rnorm(nindiv,mean=0,sd=sqrt(15))
# GENOMIC RELATIONSHIP MATRIX G
G <- (1/nmark)*W%*%t(W)
# SVD OF G
EVD <- eigen(G)
names(EVD)
head(EVD$values)
U <- EVD$vector
val <- EVD$values
val[length(y)] <-0
D <- diag(val,nrow=nindiv)
ytilde <- t(U)%*%y
dim(ytilde)
#END OF GENERATION OF DATA
###################################################################
# CODE1305(cont)
# FUNCTIONS loglik AND logliktransf TO COMPUTE LOG-LIKELIHOODS
loglik <- function(data,par)
{
  ve <- par[1]
  k <- par[2]
  ll <- -0.5*(length(ytilde)*log(ve)+sum(log(val*k+1))+
                (1/ve)*sum(ytilde^2/(val*k+1)))
  return(-ll)
}
# FUNCTION logliktransf TO COMPUTE TRANSFORMED LOG-LIKELIHOOD
logliktransf <- function(data,par)
{
  nue <- par[1]
  nug <- par[2]
  lltransf<--0.5*(length(ytilde)*nue+sum(log(val*exp(nug-nue)+1))+
                    (1/exp(nue))*sum(ytilde^2/(val*exp(nug-nue)+1)))
  
  return(-lltransf)
}
##################################################################
# FUNCTION OTPIM TO COMPARE WITH RESULTS TO COME
result1 <-optim(par=c(5,0.5),loglik,data =ytilde,hessian=TRUE)
result1$par
# OBTAIN ASYMPTOTIC VARIANCES BY INVERSION OF THE -HESSIAN
solve(result1$hessian)
# USE OPTIM TO MAXIMIZE TRANSFORMED LOG-LIKELIHOOD
result2 <-
  optim(par=c(exp(5),exp(0.5)),logliktransf,data=ytilde,hessian=TRUE)
result2$par
solve(result2$hessian)
##################################################################
# NEWTON-RAPHSON COMPUTATIONS
nit <- 20
resultnr<-matrix(data=NA,nrow=nit,ncol=3)
llike<-matrix(data=NA,nrow=nit,ncol=1)
ve <- matrix(data=NA, nrow=nit+1,ncol=1)
k <- matrix(data=NA, nrow=nit+1,ncol=1)
ve[1] <- 7
k[1] <- 0.4
llike[1] <- -loglik(ytilde,c(ve[1],k[1]))
for (i in 1:nit)
{
  vc11 <- - 0.5*((2/ve[i]^3)*sum(ytilde^2/(1+k[i]*val))-
                   length(ytilde)/ve[i]^2)
  vc12 <- -0.5*(1/ve[i]^2)*sum(val*ytilde^2/(1+k[i]*val)^2)
  vc22 <- -0.5*((1/ve[i])*sum(2*val^2*ytilde^2/(1+k[i]*val)^3)-
                  sum(val^2/(1+k[i]*val)^2))
  vcmat <- matrix(c(vc11,vc12,vc12,vc22),nrow=2,ncol=2)
  vcmatinv <- solve(vcmat)
  fd1 <- -0.5*((length(ytilde)/ve[i])-
                 (1/ve[i]^2)*sum(ytilde^2/(1+k[i]*val)))
  fd2<--0.5*(sum(val/(1+k[i]*val))-
               (1/ve[i])*sum(val*ytilde^2/((1+k[i]*val)^2)))
  fd <- matrix(c(fd1,fd2),nrow=2,ncol=1)
  sol0 <- matrix(c(ve[i],k[i]),nrow=2,ncol=1)
  sol1 <- sol0 - (vcmatinv%*%fd)
  ve[i+1] <- sol1[1]
  k[i+1] <- sol1[2]
  llike[i+1] <- loglik(ytilde,c(ve[i+1],k[i+1]))
  resultnr[i,] <- c(sol1[1],sol1[2],llike[i])
}
#PRINT RESULTS
tail(resultnr)
-vcmatinv
################################################################
# CODE1305 (cont)
# EM COMPUTATIONS
emiter<-1000
vgem<-5
veem<-10
kem<-vgem/veem
resultem<-matrix(data=NA,nrow=emiter,ncol=4)
for (i in 1:emiter)
{
  expalfa<-kem*(val/(val*kem+1))*ytilde
  tol<-0.00001
  trdiv<-vgem*sum(1/(val[val>tol]*kem+1))
  trv<- vgem*sum(val/(val*kem+1))
  vgem<-((kem^2)*sum(ytilde^2*val/(val*kem+1)^2)+trdiv)/(nindiv-1)
  veem<-(sum(y^2)-2*sum(expalfa*ytilde)+sum(expalfa^2)+trv)/nindiv
  kem<-vgem/veem
  loglike<-loglik(ytilde,c(veem,kem))
  resultem[i,]<-c(vgem,veem,kem,-loglike)
}
tail(resultem)
