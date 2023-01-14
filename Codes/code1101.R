# CODE1101
rm(list=ls()) # CLEAR WORKSPACE
set.seed(195021)
x<-seq(from=0, to=2*pi,by=0.2)
f0<-function(x){ 100+sin(2*x)+cos(x/2) }
R2<-2/3
y<-f0(x)+rnorm(n=length(x),sd=sqrt(var(f0(x))*(1-R2)/R2))
z1 <- cut(x,breaks=seq(from=min(x),to=max(x+.01),
                       length=7),right=F)
f1 <- lm(y~z1)
z2 <- cut(x,breaks=seq(from=min(x),to=max(x+.01),
                       length=20),right=F)
f2 <- lm(y~z2)
#setwd("C:/Users/au223137/Dropbox/Rsessions/MarkDown")
#pdf("C:/Users/au223137/Dropbox/Rsessions/MarkDown/Figures/
#binestimator.pdf")

plot(y~x,main='Binned estimator')
#lines(x=x,y=f0(x),col='red',lwd=2)
lines(x=x,y=predict(f1),col='blue',lwd=2)
lines(x=x,y=predict(f2),col='green',lwd=2)
legend(5, 103, legend=c("6 bins", "19 bins"),
       col=c("blue", "green"), lty=1:1, cex=0.8)
#dev.off()
