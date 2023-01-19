# CODE1103
# GAUSSIAN KERNEL REGERSSION
rm(list=ls()) # CLEAR WORKSPACE
set.seed(195021)
# GENERATE DATA 
x<-seq(from=0, to=2*pi,by=0.2)
#x<-seq(from=0, to=2*pi,length.out=33)

f0<-function(x){ 100+sin(2*x)+cos(x/2) }
R2<-2/3
y<-f0(x)+rnorm(n=length(x),sd=sqrt(var(f0(x))*(1-R2)/R2))
# CHOOSE h
h <- 0.25
# CONSTRUCT DISTANCE MATRIX AND GAUSSIAN KERNEL
d <- as.matrix(dist(x))^2
Kh25 <- exp(-(1/(2*h^2))*d)
sc25 <- apply(Kh25,1,sum)
mhgaus25 <- Kh25%*%y/sc25
plot(y~x,main='Gaussian kernel regression')
lines(x,mhgaus25,col="red")
# CHOOSE h
h <- 0.10
Kh10 <- exp(-(1/(2*h^2))*d)
sc10 <- apply(Kh10,1,sum)
mhgaus10 <- Kh10%*%y/sc10
lines(x,mhgaus10,col="blue")
legend(5, 103, legend=c("h=0.25", "h=0.10"),
       col=c("red", "blue"),lty=1:1, lwd=c(1.5,1.5), cex=0.8)
