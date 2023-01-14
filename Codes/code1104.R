# CODE1104
rm(list=ls()) # CLEAR WORKSPACE
set.seed(195021)
# LOCAL LINEAR REGRESSION SETS p = 1
p <- 1
x<-seq(from=0, to=2*pi,by=0.2)
f0<-function(x){ 100+sin(2*x)+cos(x/2) }
R2<-2/3
y<-f0(x)+rnorm(n=length(x),sd=sqrt(var(f0(x))*(1-R2)/R2))

w <- matrix(data=NA,nrow=length(x),ncol=length(x))
X <- matrix(data=NA, nrow=length(y), ncol=(p+1))
W <- matrix(data=NA, nrow=length(y), ncol=length(y))

one <- rep(1,length(y))

X <- cbind(one,x)
Xt <- t(X)
k <- seq(1:length(x))

mhx <- function(k,w,one,X){
  W <- diag(w[k,])
  Xt <- t(X)
  solve(Xt%*%W%*%X,Xt%*%W%*%y)
}

mpred <- function(k,estx,x){
  estx[1,k]+estx[2,k]*x[k]
}

# CONSTRUCT DISTANCE MATRIX 
dst <- as.matrix(dist(x))
d <- as.matrix(dist(x))^2

# CONSTRUCT GAUSSIAN KERNEL
# CHOOSE h
h <- 0.2
Kh25 <- exp(-(1/(2*h^2))*d)
div <- apply(Kh25,1,sum)

# SCALE GAUSSIAN KERNEL: PLACE RESULT IN w
for( i in 1:nrow(Kh25)){
  w[i,] <- Kh25[i,]/div[i]
}

estx <- sapply(k,mhx,w,one,X)

fitx <- sapply(k,mpred,estx,x)
fitx[which.min(abs(x-5.12))] ## APROX PREDICTION FOR A NEW X=5.12
fitx[which.min(abs(x-5.0))] ## FIT FOR X=5.0
#setwd("C:/Users/au223137/Dropbox/Rsessions/MarkDown")
#pdf("C:/Users/au223137/Dropbox/Rsessions/MarkDown/Figures/
#    localpoly.pdf")

plot(x,y,main='Local linear regression')
lines(x,fitx,lty=1,col="blue")

h <- 0.5
Kh25 <- exp(-(1/(2*h^2))*d)
div <- apply(Kh25,1,sum)

# SCALE GAUSSIAN KERNEL: PLACE RESULT IN w
for( i in 1:nrow(Kh25)){
  w[i,] <- Kh25[i,]/div[i]
}

estx <- sapply(k,mhx,w,one,X)

fitx <- sapply(k,mpred,estx,x)
fitx[which.min(abs(x-5.12))] ## APROX PREDICTION FOR A NEW X=5.12
fitx[which.min(abs(x-5.0))] ## FIT FOR X=5.0

lines(x,fitx,lty=1,col="red")
legend(5, 103, legend=c("h=0.20", "h=0.50"),
       col=c("blue", "red"),lty=1:1, lwd=c(1.5,1.5), cex=0.8)
#dev.off()