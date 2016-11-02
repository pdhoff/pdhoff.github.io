#### ---- function for plotting voxel slices 
svplot<-function(Y,zlevels=c(1,2,3),cut=0)
{
  x1<-1:dim(Y)[1]
  x2<-1:dim(Y)[2]
  x1c<-c(matrix( x1,nrow=dim(Y)[1],ncol=dim(Y)[2] ))
  x2c<-c(t(matrix( x2,nrow=dim(Y)[2],ncol=dim(Y)[1] )) )

  for(z in zlevels)
  {
    vpsig<-which( !is.na(Y[,,z]) & Y[,,z]>  cut )
    vnsig<-which( !is.na(Y[,,z]) & Y[,,z]< -cut )

    plot( DTIdata[DTIdata[,3]==z,1:2],col="pink",
           pch=16,xaxt="n",yaxt="n",xlab="",ylab="")
    points(  x1c[vpsig],x2c[vpsig] ,col="darkgreen",pch=16)
    points(  x1c[vnsig],x2c[vnsig] ,col="blue",pch=16)
  }
}



#### ---- load and tensorize data
load(url("http://statweb.stanford.edu/~ckirby/brad/LSI/datasets-and-programs/data/DTIdata.RData")) 
#load("DTIdata.RData")
DTIdata[,1]<-match(DTIdata[,1],sort(unique(DTIdata[,1])) )
DTIdata[,2]<-match(DTIdata[,2],sort(unique(DTIdata[,2])) )
DTIdata[,3]<-match(DTIdata[,3],sort(unique(DTIdata[,3])) )
m<-apply(DTIdata[,1:3],2,max)
Y<-array(dim=m)
Y[ DTIdata[,1:3] ] <- DTIdata[,4] 
IDX<-which(!is.na(Y),arr.ind=TRUE) 



#### ---- empirical Bayes for unstructured shrinkage
lambda0<- 2*sqrt(2)/sqrt(mean(Y^2,na.rm=TRUE)-1)
Theta0<-array(dim=dim(Y))
Theta0[!is.na(Y) & abs(Y)<lambda0/2 ] <-  0
Theta0[!is.na(Y) &  Y> lambda0/2] <-  Y[ !is.na(Y) &  Y>lambda0/2 ] - lambda0/2
Theta0[!is.na(Y) & -Y> lambda0/2] <-  Y[ !is.na(Y) & -Y>lambda0/2]+lambda0/2
sum( Theta0!=0 ,na.rm=TRUE ) 
svplot(Theta0,z=12:14)



#### ---- empirical Bayes for spatial autocorrelation

## -- use a subset of the data to obtain eBayes estimate of spat corr
set.seed(1) 
sIDX<-IDX[ sample(1:nrow(IDX),1500) , ]

## -- neighborhood matrix 
G<-matrix(0,nrow(sIDX),nrow(sIDX))
for(i in 1:nrow(sIDX))
{
  neighbors<-which(apply(abs(sweep(sIDX,2,sIDX[i,],"-")),1,sum)==1)
  G[i,neighbors]<-G[neighbors,i]<-1/length(neighbors)
}
neighborblock<-rbind(diag(3),-diag(3))

## -- find eBayes correlation estimate
Theta<- Y/(1+1/(mean(Y^2,na.rm=TRUE)-1) ) 
v<-sign( Theta[sIDX] ) * sqrt( abs(Theta[sIDX]) )   

llCAR<-function(theta)
{
  rho<-1/(1+exp(-theta[1])) ; t2<-exp(theta[2]) ; p<-nrow(G)  
  Psi<-(diag(nrow(G))-rho*G)/t2 
  ePsi<-eigen(Psi)  
  m2ll<-Inf
  if(all(ePsi$val>0))
  { 
    m2ll<- (t(v)%*%Psi%*%v) - sum(log(ePsi$val)) 
  }
  cat(rho,t2,m2ll,"\n")
  m2ll
}

fit<-optim( c(0,-2), llCAR )
rho<- 1/(1+exp(-fit$par[1])) 
icvar<-1/exp(fit$par[2])



#### ---- BCD algorithm
Theta<-Y
U<- sqrt(abs(Theta))
V<-Theta/U
del<-1 
nit<-0 
cut<-min( abs(Theta0[Theta0!=0]) ,na.rm=TRUE)/10
PRG<-NULL

while(del>1e-10)
{

  for(vox in 1:nrow(IDX))
  { 
    i<-IDX[vox,1] ; j<-IDX[vox,2] ; k<-IDX[vox,3] 

    ## -- get neighborhood
    neighbors<-sweep(neighborblock,2,c(i,j,k),"+"  )
    neighbors<-neighbors[ apply( neighbors>0,1,all) &
                          neighbors[,1]<=m[1] &
                          neighbors[,2]<=m[2] &
                          neighbors[,3]<=m[3]  , ] 

    ## -- update U 
    q<-V[i,j,k]^2 + icvar ; l<-Y[i,j,k]*V[i,j,k]
    l<-l + icvar*rho*sum(U[neighbors],na.rm=TRUE)/max(1,nrow(neighbors)) 
    U[i,j,k] <- l/q 

    ## -- update V
    q<-U[i,j,k]^2 + icvar ; l<-Y[i,j,k]*U[i,j,k]
    l<-l + icvar*rho*sum(V[neighbors],na.rm=TRUE)/max(1,nrow(neighbors)) 
    V[i,j,k] <- l/q

  }

  UV<-U*V 
  del<-mean( (UV - Theta)^2 ,na.rm=TRUE)/mean(Theta^2,na.rm=TRUE)  
  Theta<-UV 
  nit<-nit+1 
  PRG<-rbind(PRG, c(del,sum(abs(Theta)>cut,na.rm=TRUE) ) ) 
  cat(nit,del,sum(abs(Theta)>cut,na.rm=TRUE),"\n")  
  par(mfrow=c(1,3))
  plot(log(PRG[,1]),type="l") ; plot(PRG[,2],type="l") 
  if(nit>1) {  plot( log(1-diff(PRG[,2])),type="l")  } 
  if(nit%%10==0){ svplot(Theta,cut=cut,z=12:14) }
}



#### ---- plot results 
pdf("brainraw.pdf",height=2,width=6)
par(mfrow=c(1,3),mar=c(1,1,1,1),mgp=c(0,0,0))


x1c<-c(matrix( 1:dim(Y)[1],nrow=dim(Y)[1],ncol=dim(Y)[2] ))
x2c<-c(t(matrix( 1:dim(Y)[2],nrow=dim(Y)[2],ncol=dim(Y)[1] )) )
for(z in 12:14)
{
  vp<-which( !is.na(Y[,,z]) & Y[,,z]>  qnorm(.9) )
  vn<-which( !is.na(Y[,,z]) & Y[,,z]<  qnorm(.1) )
  vpsig<-which( !is.na(Y[,,z]) & Y[,,z]>  qnorm(.975) )
  vnsig<-which( !is.na(Y[,,z]) & Y[,,z]<  qnorm(.025) )

  plot(DTIdata[DTIdata[,3]==z,1:2],col="pink",
       pch=16,xaxt="n",yaxt="n",xlab="",ylab="")
  points(x1c[vp],x2c[vp] ,col="green",pch=16)
  points(x1c[vn],x2c[vn] ,col="lightblue",pch=16)
  points(x1c[vpsig],x2c[vpsig] ,col="darkgreen",pch=16)
  points(x1c[vnsig],x2c[vnsig] ,col="blue",pch=16)
}
dev.off() 

pdf("brainshrink.pdf",height=4,width=6)
par(mfrow=c(2,3),mar=c(1,1,1,1),mgp=c(0,0,0))
svplot(Theta,cut=cut*10,z=12:14)
svplot(Theta0,z=12:14)
dev.off() 



