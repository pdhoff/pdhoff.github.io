## ----readData,echo=FALSE------------------------------------------------------

## Load the following objects:
# X : vectorized source EEMs 
# Y : vectorized mixed EEMs 
# B0 : mixing proportions (abundances) 
# idx : indicator matrix of meaningful pixels 
# EEM2vec : vectorization function for an EEM
# vec2EEM : inverse of EEM2vec 
load("DOMmix.RData") 


p<-dim(X)[[1]]
s<-dim(X)[[2]]
n<-dim(X)[[3]] 
m<-dim(Y)[[2]]


## ----procVar,fig.height=2.75,fig.width=5,out.width='100%',fig.align='center',echo=FALSE----

## Plots of across-replicate variation 

# sample average EEMs
MX<-apply(X,c(1,2),mean) 
MY<-apply(Y,c(1,2),mean) 

# plot replicate EEMs versus average EEMs
par(mfrow=c(2,3),mar=c(3,3,1,1),mgp=c(1.75,.75,0)) 
for(k in 1:2){   
  for(i in 1:n){ 
    plot(MY[,k],Y[,k,i],xlab="averaged EEM",ylab="replicate EEM") 
    abline(0,1,col="green" ) }
} 


## ----procVarEst,echo=FALSE,size="normalsize",results="asis"-------------------

## Estimate procedural effects and variation 

# divide EEMs by sample average to estimate multiplicative effect 
RX<-sweep(X,c(1,2),MX,"/") ; RX[is.na(RX)]<-0 ; AX<-apply(RX,c(2,3),mean) 
RY<-sweep(Y,c(1,2),MY,"/") ; RY[is.na(RY)]<-0 ; AY<-apply(RY,c(2,3),mean) 

# create a summary table 
RAX<-NULL
for(k in 1:s){ 
  Xk<-X[,k,] 
  rho<-cor(X[,k,],MX[,k]) 
  RAX<-rbind(RAX,c(rho,AX[k,],sd(AX[k,]) ) )
}

RAY<-NULL
for(k in 1:m){ 
  Yk<-Y[,k,] 
  rho<-cor(Y[,k,],MY[,k]) 
  RAY<-rbind(RAY,c(rho,AY[k,],sd(AY[k,]) ) )
}

RA<-rbind(RAX,RAY) 
colnames(RA)<-c("$\\hat \\rho_1$","$\\hat \\rho_2$","$\\hat\\rho_3$","$\\hat a_1$","$\\hat a_2$","$\\hat a_3$","$\\hat \\sigma_a$") 
rownames(RA)<-c("$s_1$","$s_2$","$s_3$",
                "$m_1$","$m_2$","$m_3$","$m_4$","$m_5$","$m_6$","$m_7$") 

knitr::kable(round(RA[,1:7],3),
     align="|ccc|ccc|c",format="latex",escape=FALSE,linesep=c(rep("",2),
     "\\hline",rep("",8)) ) 

# pooled sample estimate of procedural variance
s2a<-mean( c(apply(AX,1,var),apply(AY,1,var)) )  


## ----fitRes,fig.height=2.25,fig.width=5,out.width='100%',fig.align='center',echo=FALSE----
## Fitted versus residual plots

# center replications with EEM-level multiplicative effect 
CX<-sweep(X,c(2,3),AX,"/") 
CY<-sweep(Y,c(2,3),AY,"/") 

# plot sample average EEMs versus  additive residuals 
par(mfrow=c(1,3),mar=c(3,3,1,1),mgp=c(1.75,.6,0)) 
for(k in 1:s){ 
  plot(range(MX[,k]),range( sweep(CX[,k,],c(1),MX[,k],"-")),type="n",
    xlab=expression(hat(mu)), ylab=expression( y/hat(a)-hat(mu))) 
  for(i in 1:n){ points(MX[,k] , CX[,k,i] - MX[,k] ) } 
  abline(h=0,col="green")  
}


## ----mvRel,fig.height=2.25,fig.width=5,out.width='100%',fig.align='center',echo=FALSE----

## Mean-variance relationship 

# pixel-specific across-replicate standard deviation 
SX<-apply(X,c(1,2),sd) 
SY<-apply(Y,c(1,2),sd) 

# plots of mean versus standard deviation
par(mfrow=c(1,3),mar=c(3,3,1,1),mgp=c(1.75,.75,0)) 
for(k in 1:3){   
  plot(MX[,k],SX[,k],xlab="replication average",
                     ylab="replication standard deviation") 
  lines(lowess(MX[,k],SX[,k]),col="green" )
}



## ----expVar,echo=FALSE,size="normalsize",results="asis"-----------------------

## Estimate measurement variation

# sweep out mean and procedural effects to get measurement residuals
EX<-sweep( sweep(X,c(2,3),AX,"/") ,c(1,2),MX,"/") ; EX[is.na(EX)]<-1
EY<-sweep( sweep(Y,c(2,3),AY,"/") ,c(1,2),MY,"/") ; EY[is.na(EY)]<-1

# measurement variance estimates 
s2eX<-apply(EX,2,sd)^2
s2eY<-apply(EY,2,sd)^2 

# SNR 
snrX<-(1/( (s2a+1)*(s2eX+1) -1 ))
snrY<-(1/( (s2a+1)*(s2eY+1) -1 ))

# table
STAB<-round(cbind( c(s2eX,s2eY),c(snrX,snrY) ),3 )
rownames(STAB)<-c("$s_1$","$s_2$","$s_3$",
                "$m_1$","$m_2$","$m_3$","$m_4$","$m_5$","$m_6$","$m_7$") 
colnames(STAB)<-c("$\\hat \\sigma_e$","$\\hat{\\text{SNR}}$") 
knitr::kable( STAB,  align="|ccc|",format="latex",escape=FALSE)


## ----htest,fig.height=5.5,fig.width=5,out.width='100%',fig.align='center',echo=FALSE----

## variance of sample averages using sample=specific estimates of s2e
VXB<-MX^2*( (s2a+1)*(matrix(s2eX,p,s,byrow=TRUE)+1) -1 )/n 
VYB<-MY^2*( (s2a+1)*(matrix(s2eY,p,m,byrow=TRUE)+1) -1 )/n 

Z<-P<-U<-NULL
FDR<-.05
par(mfrow=c(4,2),mar=c(3,3,1,1),mgp=c(1.75,.75,0)) 
for(k in 1:m){ 
 
  ## z-statistic 
  zn<- MY[,k] - MX%*%B0[k,] 
  zd<-sqrt(VYB[,k] + apply(sweep(VXB,2,B0[k,]^2,"*"),1,sum )  )
  z<-zn/zd  

  ## BH procedure
  pv<-2*(1-pnorm(abs(z))) 
  opv<-sort(pv) 
  rc<-max( which( opv< FDR*(1:length(pv))/length(pv) )) 
  id<- which( pv <= opv[rc] ) 
  u<-sign(z)*(pv<= opv[rc] ) 

  ## plot significant deviations in EEMs
  uEEM<-vec2EEM(u,idx) 
  image(1:ncol(uEEM),1:nrow(uEEM),t(uEEM) , xaxt="n",yaxt="n",xlab="",ylab="",
        col=c("blue","gray","yellow"))  
  axis(1,at=1:ncol(uEEM),labels=colnames(uEEM),tick=FALSE,las=1) 
  axis(2,at=1:nrow(uEEM),labels=rownames(uEEM),tick=FALSE,las=1)
  mtext(paste0("m",k),side=3,cex=.7) 
  P<-cbind(P,pv) ; Z<-cbind(Z,z) ; U<-cbind(U,u) 
}


## fluorescence/p-value relationship
plot(MY,log(P),col=c("blue","gray","yellow")[U+2] ,
     ylab="log p-value",xlab=expression(bar(y)),pch="." )



## ----estVar,echo=FALSE,fig.height=3.75,fig.width=5.5--------------------------
## All combinations of mixtures and endmembers 
BHAT<-array(dim=c(m,3,81))
for(k in 1:m){ 
  l<-0
  for(iy in 1:n){ for(ix1 in 1:n){ for(ix2 in 1:n){ for(ix3 in 1:n) {
  l<-l+1
  yl<-Y[,k,iy] 
  Xl<-cbind( X[,1,ix1], X[,2,ix2],X[,3,ix3] ) 
  BHAT[k,,l]<- nnls::nnls(Xl,yl)$x  
  }}}}
} 

O<-matrix(rep(c(-1,0,1),times=rep(m,s))/100,m,s) 
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0)) 
plot(range(B0+O),range(BHAT),type="n",
      xlab="true abundances",ylab="estimated abundances") ; abline(0,1) 

for(i in 1:dim(BHAT)[3]){ 
points((B0+O)[,2],BHAT[,2,i],col="green",pch=2)  
points((B0+O)[,3],BHAT[,3,i],col="lightblue",pch=3) 
points((B0+O)[,1],BHAT[,1,i],col="black",pch=1) 
}

segments((B0+O)[,2],apply(BHAT[,2,],1,min),  
         (B0+O)[,2],apply(BHAT[,2,],1,max),col="green") 
segments((B0+O)[,3],apply(BHAT[,3,],1,min),  
         (B0+O)[,3],apply(BHAT[,3,],1,max),col="lightblue") 
segments((B0+O)[,1],apply(BHAT[,1,],1,min),  
         (B0+O)[,1],apply(BHAT[,1,],1,max),col="black") 

legend("topleft",pch=c(1,2,3),col=c("black","green","lightblue"),
    legend=c("groundwater","streamwater","wastewater") )


## ----estVarTab,echo=FALSE,size="normalsize",results="asis"--------------------
MB<- apply(BHAT,c(1,2),mean)
SB<- apply(BHAT,c(1,2),sd)

TB<-cbind(B0,MB,SB)[,c(1,4,7,2,5,8,3,6,9) ]
rownames(TB)<-c("$m_1$","$m_2$","$m_3$","$m_4$","$m_5$","$m_6$","$m_7$")  
colnames(TB)<-c("$b_1$","mean($\\hat b_1$)","sd($\\hat b_1$)",
  "$b_2$","mean($\\hat b_2$)","sd($\\hat b_2$)",
"$b_3$","mean($\\hat b_3$)","sd($\\hat b_3$)" )

library(kableExtra) 
kable(round(TB,3),
   align="ccc|ccc|ccc",format="latex",escape=FALSE,
   linesep=c(rep("",9) ))  %>%
  add_header_above(c(" "=1, "groundwater" = 3, "streamwater" = 3, "wastewater" = 3))


## ----echo=FALSE---------------------------------------------------------------
plotEEM<-function(x){   
  fields::image.plot( 
    as.numeric(dimnames(idx)[[2]]),as.numeric(dimnames(idx)[[1]]),
    t(vec2EEM(x,idx)),
    col=topo.colors(10),  xlab="Ex. (nm)", ylab="Ex. (nm)",   
    legend.lab="QSU",legend.cex=.8)
}


## ----echo=FALSE,fig.height=4,fig.width=5.5,fig.align='center',out.width='5in'----
enames<-c(expression(s[1]),expression(s[2]),expression(s[3]))
enames<-c(expression(m[1]),expression(m[2]),expression(m[3]))
par(mfrow=c(2,2),mar=c(3,4,1,1.5),mgp=c(1.75,.75,0),oma=c(0,1,0,2)) 
for(k in 1:3){ plotEEM(Y[,k,1]) ; mtext(enames[k],side=3,cex=.8) }


## ----echo=FALSE,fig.height=8,fig.width=5.5,fig.align='center',out.width='5in'----
enames<-c(expression(m[1]),expression(m[2]),expression(m[3]),expression(m[4]),
          expression(m[5]),expression(m[6]),expression(m[7]))
par(mfrow=c(4,2),mar=c(3,4,1,1.5),mgp=c(1.75,.75,0),oma=c(0,1,0,2))
for(k in 1:7){ plotEEM(Y[,k,1]) ; mtext(enames[k],side=3,cex=.8) }

