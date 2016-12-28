#### ---- functions and data
source("functions.R")
radon<-dget("http://www.stat.washington.edu/~pdhoff/courses/560/Data/mn.radon.data")

#### ---- data description
dim(radon) 
g<-radon[,1]
y<-radon[,2]
ng<-table(g)
ybg<-tapply(y,g,mean) 
sdg<-tapply(y,g,sd) 


#### ---- simulation study 
thetag<-tapply(y,g,mean)
sg<-tapply(y,g,sd) ; sg[is.na(sg)]<-0

tmp<-rep(0,length(ng)) ; tmp[ng==1]<-NA 
bcover<-ucover<-fcover<-bwidth<-uwidth<-fwidth<-tmp 

nsim<-2500
alpha<-.05
for(s in 1:nsim)
{ 
  
  ## -- simulate data
  set.seed(s+nsim*(sgroup-1))
  ysim<-NULL
  for(j in 1:length(ng)){  ysim<-c(ysim,rnorm(ng[j],thetag[j],sg[j])) } 
  ybgsim<-tapply(ysim,g,mean)
  sdgsim<-tapply(ysim,g,sd)


  ## -- construct intervals
  BCI<-UCI<-FCI<-NULL  
  for(j in 1:length(ng))
  { 

    ## -- data from this group
    bci<-uci<-fci<-c(-Inf,Inf)
    if(ng[j]>1)
    {

      ## -- eBayes estimates from other groups 
      gamj<-ebayes_fabt(ysim[g!=j],g[g!=j])
 
      ## -- FAB-t interval 
      fci<-fabt_ci(ybgsim[j],sdgsim[j],ng[j],
                    gamj$mu,gamj$t2,gamj$nu0,gamj$s20,alpha=alpha)
 
      ## -- UMA-t interval
      uci<- ybgsim[j]+c(-1,1)*sdgsim[j]*qt(1-alpha/2,ng[j]-1)/sqrt(ng[j]) 

      ## -- Bayes credible intervals
      ethj<-(gamj$mu/gamj$t2+ybgsim[j]*ng[j]/sdgsim[j]^2)/
            (1/gamj$t2+ng[j]/sdgsim[j]^2)
      psdj<-1/sqrt(1/gamj$t2 + ng[j]/sdgsim[j]^2 )   
      bci<-ethj+c(-1,1)*psdj*qt(1-alpha/2,ng[j]-1)
      
    }
    BCI<-rbind(BCI,bci)  
    UCI<-rbind(UCI,uci)
    FCI<-rbind(FCI,fci) 
  }

  
  fcover<-fcover+1*( FCI[,1] < thetag & thetag < FCI[,2] )    
  ucover<-ucover+1*( UCI[,1] < thetag & thetag < UCI[,2] )    
  bcover<-bcover+1*( BCI[,1] < thetag & thetag < BCI[,2] )

  fwidth<-fwidth + (FCI[,2]-FCI[,1]) 
  uwidth<-uwidth + (UCI[,2]-UCI[,1]) 
  bwidth<-bwidth + (BCI[,2]-BCI[,1])

  cat(s,"\n") 
  if(s%%10==0){ save(fcover,ucover,bcover,fwidth,uwidth,bwidth,s,sgroup,
                     file=paste0("simres",sgroup,"_alpha",alpha))  }


  ## -- simulation study plot  
  par(mfrow=c(1,2),mar=c(3,3,1,1),mgp=c(1.75,.75,0))

  pbcover<-bcover/s 
  plot(thetag,pbcover,xlab=expression(theta),ylab="credible interval coverage",pch=16,
       col="gray")
  segments(thetag,pbcover+1.96*sqrt(pbcover*(1-pbcover)/s),
         thetag,pbcover-1.96*sqrt(pbcover*(1-pbcover)/s) ) 
  abline(h=1-alpha,lty=2)

  plot(thetag,fwidth/uwidth,ylim=range(c(fwidth/uwidth,bwidth/uwidth),na.rm=T),pch=16,
       xlab=expression(theta),ylab="width ratio")
  points(thetag,bwidth/uwidth,col="gray",pch=16)
  segments(thetag,fwidth/uwidth,thetag,bwidth/uwidth)
  abline(h=1,lty=2)
  legend(2.25,1.2,legend=c("FAB","Bayes"),pch=c(16,16),col=c("black","gray")  )
}



