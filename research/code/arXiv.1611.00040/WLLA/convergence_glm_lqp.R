## --  functions and data
source("functions_glm_lqp.R")
source("functions_lla.R") 

K<-4
p<-100
n<-round(p*1.5) 
ns<-100
yS<-lambdaS<-NULL
XS<-array(dim=c(n,p,ns))
betaT<-betaS<-matrix(nrow=ns,ncol=p) 

for(s in 1:ns)
{
  set.seed(s) 
  beta0<-rnorm(p)*rbinom(p,1,.5)*.5
  X<-matrix(rnorm(n*p),n,p) 
  y<-rbinom(n,1,exp( X%*%beta0 )/(1+ exp(  X%*%beta0 )  ) )

  beta_hat<-apply(X,2,function(x){ glm(y~ -1+x,family="binomial")$coef } )
  t2<-var(beta_hat)
  qpow<-2/K
  lambda<-2*(  gamma(3/qpow)/( gamma(1/qpow)*t2 ) )^(qpow/2)

  yS<-cbind(yS,y)
  lambdaS<-c(lambdaS,lambda)
  XS[,,s]<-X
  betaS[s,]<-beta_hat
  betaT[s,]<-beta0
}
 

## --  evaluate convergence 
NIT<-OBJ<-BD<-SP<-SD<-NULL
for(s in 1:ns)
{
  fith<-hpp_glm(yS[,s],XS[,,s],lambda=lambdaS[s],K=K)
  fitl<-lqa_glm(yS[,s],XS[,,s],lambda=lambdaS[s],q=2/K)
  fits<-lla_glm(yS[,s],XS[,,s],lambda=lambdaS[s],q=2/K)

  ## -- number of iterations and objective at convergence
  objh<-fith$obj
  objl<-fitl$obj 
  objs<-fits$obj 
  NIT<-rbind(NIT, c(length(objh),length(objl),length(objs) ) )
  OBJ<-rbind( OBJ,c( min(objh), min(objl),min(objs) ) )

  ## -- checking for differences in beta estimates 
  bh<-fith$beta
  bl<-fitl$beta
  BD<-c(BD,  mean(abs(bh-bl))/mean(abs(betaT[s,]))  )

  ## -- checkng sparsity pattern 
  sbh<-1*(abs(bh)>1e-6)
  sbl<-1*(abs(bl)>1e-6)
  sparsity<-c( sum(sbh*sbl),sum(sbh*(1-sbl)),sum((1-sbh)*sbl),
               sum((1-sbh)*(1-sbl)) )
  mbh_blz<-mean( abs(bh[ sbh==1 & sbl==0 ] ) )
  mbl_bhz<-mean( abs(bl[ sbh==0 & sbl==1 ] ) )
  SD<-rbind(SD,c(mbh_blz,mbl_bhz))
  SP<-rbind(SP,sparsity)

  cat(s,"\n")
}


## -- summary 
range( (OBJ[,2]-OBJ[,1])/abs(OBJ[,1]) )
range( (OBJ[,3]-OBJ[,1])/abs(OBJ[,1]) )

range(BD)
apply(NIT,2,median)
mean( NIT[,2]/NIT[,1] ) 

mean( SP[,2]+SP[,3] > 0 ) 
mean( OBJ[ SP[,2]+SP[,3] > 0 ,2 ] > OBJ[ SP[,2]+SP[,3] > 0 ,1 ] )
mean( (SP[,2]+SP[,3])/p ) 
mean( SD,na.rm=TRUE)


## -- convergence trajectory
OBJH<-OBJL<-OBJS<-NULL
iter<-50
for(s in 1:ns)
{
  objh<-hpp_glm(yS[,s],XS[,,s],lambda=lambdaS[s],itmin=iter,itmax=iter,K=K)$obj
  objl<-lqa_glm(yS[,s],XS[,,s],lambda=lambdaS[s],itmin=iter,itmax=iter,q=2/K)$obj
  objs<-lla_glm(yS[,s],XS[,,s],lambda=lambdaS[s],itmin=iter,itmax=iter,q=2/K)$obj

  OBJH<-rbind( OBJH, objh[-1] )
  OBJL<-rbind( OBJL, objl[-1] )
  OBJS<-rbind( OBJS, objs[-1] )
}


## -- compute average progress towards minimum 
mnobj<-pmin( apply(OBJH,1,min,na.rm=TRUE) ,  apply(OBJL,1,min,na.rm=TRUE) ,
             apply(OBJS,1,min,na.rm=TRUE) )
mxobj<-pmax( apply(OBJH,1,max,na.rm=TRUE) ,  apply(OBJL,1,max,na.rm=TRUE) ,
             apply(OBJS,1,max,na.rm=TRUE) )
aph<-apply(sweep(sweep(-OBJH,1,mxobj,"+") ,1, mxobj-mnobj,"/" ), 2,mean)[-iter]
apl<-apply(sweep(sweep(-OBJL,1,mxobj,"+") ,1, mxobj-mnobj,"/" ), 2,mean)[-iter]
aps<-apply(sweep(sweep(-OBJS,1,mxobj,"+") ,1, mxobj-mnobj,"/" ), 2,mean)[-iter]



## -- plots
pdf(paste0("glm_K",K,"p_results.pdf"),family="Times",height=3.5,width=8)
par(mfrow=c(1,2),mar=c(3,3,1,1),mgp=c(1.75,.75,0) )

plot(log10(1-aph),type="l",
     xlab="iteration",
     ylab=expression(log[10] (1-w)),col="black",
     ylim=range(c(log10(1-aph)), log10(1-apl) )  )
lines( log10(1-apl),lty=2,col="black" )
lines( log10(1-aps),lty=1,col="gray" )


legend(30,0,lty=c(1,2,1),col=c("black","black","gray"),
      legend=c("HPP","LQA","LLA"),bty="n" )


boxplot(log10(NIT[,1:2]%*%diag(c(K,1))),names=c("HPP","LQA"),col="gray",
        ylab=expression(log[10] ~ "IRLS optimizations")) 

dev.off() 


## -- timings
timeh<-system.time(for(s in 1:ns){hpp_glm(yS[,s],XS[,,s],lambdaS[s],K=K)}) 
timel<-system.time(for(s in 1:ns){lqa_glm(yS[,s],XS[,,s],lambdaS[s],q=2/K)})
times<-system.time(for(s in 1:ns){lla_glm(yS[,s],XS[,,s],lambdaS[s],q=2/K)})


## -- one inner step per iteration 
timeh1<-system.time(for(s in 1:ns){hpp_glm(yS[,s],XS[,,s],lambdaS[s],K=K,iitmax=1)})
timel1<-system.time(for(s in 1:ns){lqa_glm(yS[,s],XS[,,s],lambdaS[s],q=2/K,iitmax=1)})
times1<-system.time(for(s in 1:ns){lla_glm(yS[,s],XS[,,s],lambdaS[s],q=2/K,iitmax=1)}) 


## -- check convergence of this approach
NIT1<-OBJ1<-NULL 
for(s in 1:ns)
{  
  fith<-hpp_glm(yS[,s],XS[,,s],lambda=lambdaS[s],K=K,iitmax=1)
  fitl<-lqa_glm(yS[,s],XS[,,s],lambda=lambdaS[s],q=2/K,iitmax=1)
  fits<-lla_glm(yS[,s],XS[,,s],lambda=lambdaS[s],q=2/K,iitmax=1)

  ## -- number of iterations and objective at convergence
  objh<-fith$obj
  objl<-fitl$obj
  objs<-fits$obj
  NIT1<-rbind(NIT1, c(length(objh),length(objl),length(objs) ) )
  OBJ1<-rbind( OBJ1,c( min(objh), min(objl),min(objs) ) ) 
  cat(s,"\n") 
}

## -- 
plot(OBJ[,1],OBJ1[,1]) ; abline(0,1) 
points(OBJ[,1],OBJ1[,2],col="green") 
points(OBJ[,1],OBJ1[,3],col="red") 

apply( (OBJ1-OBJ)/(OBJ) ,2,mean)*100




