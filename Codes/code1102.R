# CODE1102
rm(list=ls()) # CLEAR WORKSPACE
set.seed(195021)
x<-seq(from=0, to=2*pi,by=0.2)
f0<-function(x){ 100+sin(2*x)+cos(x/2) }
R2<-2/3
y<-f0(x)+rnorm(n=length(x),sd=sqrt(var(f0(x))*(1-R2)/R2))
d <- as.matrix(dist(x))
h <- 0.2
d2 <- ifelse(d <= h,1,0)
div <- apply(d2,1,sum)
rx <- d2%*%y/div


plot(y~x,main='Uniform kernel using a grid of values of x')
lines(x,rx,col="blue")
h <- 0.8
d2 <- ifelse(d <= h,1,0)
div <- apply(d2,1,sum)
rx <- d2%*%y/div
lines(x,rx,col="red")
legend(5, 103, legend=c("h=0.2", "h=0.8"),
       col=c("blue","red"),lty=1:1, lwd=c(1.5,1.5), cex=0.8)
