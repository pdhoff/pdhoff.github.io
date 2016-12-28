#### ---- 
alpha<-.05
source("functions.R")



#### ----
pdf("fig_fabz.pdf",height=2.5,width=7.5,family="Times")
par(mar=c(3,3,0.1,1),mgp=c(1.75,.75,0))
layout(matrix(c( 1,1,1,2,2,2,3,3,4),3,3))


## -- setups for plots
mu<-0
t2s<-c(1/4,1,4)
thetas<-seq(-3,3,length=200) 


## -- compute w functions
W<-NULL
for(theta in thetas) 
{
  w<-NULL 
  for(k in 1:3){ w<-c(w,wfabz(theta,mu,t2s[k],alpha) ) }
  W<-rbind(W,c(w))
}

## -- plot w functions
plot(range(thetas),range(W),
     type="n",xlab=expression(theta),ylab=expression(w(theta)))
for(k in 1:3){ lines(thetas,W[,k] ,col=gray(k/4),lwd=2) }
legend(-3,1,lwd=2,col=gray(3:1/4),
       legend=c(expression(tau^2==4),expression(tau^2==1),
                expression(tau^2==1/4)),bty="n")


## -- compute FAB CI 
plot(range(thetas),range(thetas),type="n",xlab="y",ylab="CI endpoints")
for(k in 1:3)
{
  YUL<-NULL
  for(y in thetas){ YUL<-rbind(YUL, c(fabz_ci(y,mu,t2s[k],alpha))) }
  matplot(thetas,YUL,type="l",lty=1,col=gray(k/4),add=TRUE,lwd=2)
}

abline(0,1,col="gray")
abline(h=0,col="gray")
abline(v=0)
abline( qnorm(1-alpha/2),1,lty=2 )
abline( qnorm(alpha/2),1,lty=2 )


## -- compute expected widths 
RISK<-NULL
for(theta in thetas)
{
  ew<-NULL ;for(k in 1:3){ ew<-c(ew,fabz_ew(theta,t2s[k],alpha=alpha) ) }
  RISK<-rbind(RISK,ew)
}

## -- plot risk 
plot(range(thetas),range(RISK),type="n",xlab="",xaxt="n",ylab="expected width")
for(k in 1:3){ lines(thetas,RISK[,k] ,col=gray(k/4),lwd=2)  }
abline(h=2*qnorm(1-alpha/2),lty=2)

## -- plot priors
plot(range(thetas),c(0,.8),
     type="n",xlab=expression(theta),ylab=expression(pi(theta)))
for(k in 1:3){ lines(thetas,dnorm(thetas,0,sqrt(t2s[k])),col=gray(k/4),lwd=2)}



#### ---- done plotting
dev.off()



## -- relative width and expected width 
RRISK<-NULL
t2<-.25

diff(fabz_ci(0,0,t2,alpha))
diff(fabz_ci(0,0,t2,alpha))/(2*qnorm(1-alpha/2))

fabz_ew(0,t2,alpha)/(2*qnorm(1-alpha/2))

diff(fabz_ci(0,0,t2,.5))/(2*qnorm(.75))
fabz_ew(0,t2,.5)/(2*qnorm(.75))




for(alpha in seq(.01,.5,by=.02))
{
  RRISK<-rbind(RRISK,
               c(alpha,diff(fabz_ci(0,0,t2,alpha))/(2*qnorm(1-alpha/2)),
                            fabz_ew(0,t2,alpha)/(2*qnorm(1-alpha/2)) ))
}




