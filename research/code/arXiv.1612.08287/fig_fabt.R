#### ---- 
alpha<-.05
source("functions.R")



#### ----
pdf("fig_fabt.pdf",height=2.5,width=7.5,family="Times")
par(mar=c(3,3,0.1,1),mgp=c(1.75,.75,0))
layout(matrix(c( 1,1,1,2,2,2,3,3,4),3,3))


## -- setups for plots
mu<-0
t2s<-c(.25,1,4)
thetas<-seq(-3,3,length=200) 
n<-10 ; nu0<-2 ; s20<-n 

## -- compute w functions
W<-NULL
for(theta in thetas) 
{
  w<-NULL 
  for(k in 1:3){ w<-c(w,wfabt(theta,mu,t2s[k],nu0,s20,n,alpha) ) }
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
s<-sqrt(n) 
plot(range(thetas),range(thetas),type="n",
     xlab=expression(bar(y)),ylab="CI endpoints")
for(k in 1:3)
{
  YUL<-NULL
  for(y in thetas){YUL<-rbind(YUL, c(fabt_ci(y,s,n,mu,t2s[k],nu0,s20,alpha)))}
  matplot(thetas,YUL,type="l",lty=1,col=gray(k/4),add=TRUE,lwd=2)
}

abline(0,1,col="gray")
abline(h=0,col="gray")
abline(v=0)
abline(  qt(1-alpha/2,n-1)*s/sqrt(n),1,lty=2 )
abline( qt(alpha/2,n-1)*s/sqrt(n),1,lty=2 )


## -- compute expected widths 
if( !is.element("ERISK", system("ls",intern=TRUE)) )
{
  s2<-n
  ERISK<-NULL
  thetass<-seq(0,3,length=10) 
  for(theta in thetass)
  {
    ew<-NULL 
    for(k in 1:3){ ew<-c(ew,fabt_ew(theta,s2,n,mu,t2s[k],nu0,s20,alpha=alpha))}
    ERISK<-rbind(ERISK,ew)
    cat(theta,ERISK[nrow(ERISK),],"\n") 
  } 
  save(ERISK,file=paste0("ERISK",100*(1-alpha))) 
} else{ thetass<-seq(0,3,length=10) ; load("ERISK") }

## polynomial apporoximation  
ERISK<-rbind(ERISK[nrow(ERISK):1,],ERISK) 
thetass<-c( -rev(thetass),thetass) 
RISK<-NULL
for(k in 1:3)
{ 
  fit<-lm( ERISK[,k]~ thetass + I(thetass^2) + I(thetass^3))$coef 
  RISK<-cbind(RISK, cbind(1,thetas,thetas^2,thetas^3)%*%fit ) 
}

## -- plot risk 
plot(range(thetas),range(RISK),type="n",xlab="",xaxt="n",ylab="expected width")
for(k in 1:3){ lines(thetas,RISK[,k] ,col=gray(k/4),lwd=2)  }
abline(h=2*qt(1-alpha/2,n-1),lty=2)

## -- plot priors
plot(range(thetas),
      c(0,.8),type="n",xlab=expression(theta),ylab=expression(pi(theta)))
for(k in 1:3){ lines(thetas,dnorm(thetas,0,sqrt(t2s[k])),col=gray(k/4),lwd=2)}


#### ---- done plotting
dev.off()


