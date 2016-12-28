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


#### ---- EDA
anova(lm(y~as.factor(g)) )        

library(lme4)
summary(lm(y~as.factor(g)))
s2<-summary(lm(y~as.factor(g)))$sigma^2
fit_lme<-lmer( y ~ 1 + (1|g) )

## --  Levene's test
d<-abs(y-tapply(y,g,median)[g])
anova( lm(d~as.factor(g)) )


#### ---- EFAB t intervals for heteroscedastic groups 
alpha<-.05
UCI<-FCI<-NULL  
for(j in 1:length(ng))
{ 

  uci<-fci<-c(-Inf,Inf)
  if(ng[j]>1)
  {
   
    ## -- eBayes estimates from other groups 
    gamj<-ebayes_fabt(y[g!=j],g[g!=j]) 

    ## -- FAB-t interval 
    fci<-fabt_ci(ybg[j],sdg[j],ng[j],gamj$mu,gamj$t2,gamj$nu0,gamj$s20,alpha=alpha)
 
    ## -- UMA-t interval
    uci<- ybg[j]+c(-1,1)*sdg[j]*qt(1-alpha/2,ng[j]-1)/sqrt(ng[j]) 
  }

  UCI<-rbind(UCI,uci)
  FCI<-rbind(FCI,fci) 
  cat(j,"\n") 

}

## -- plot
pdf("radon.pdf",height=4,width=9,family="Times")  
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0)) 
bg<-which(ng>=3)
ybarg<-tapply(y,g,mean)
plot(range(ybarg[bg]),range(c(UCI[bg,],FCI[bg,])),type="n", 
      xlab=expression(bar(y)),ylab=expression(theta)) 

segments(ybarg[bg],UCI[bg,1],ybarg[bg],UCI[bg,2],col="lightgray",lwd=4) 
segments(ybarg[bg],FCI[bg,1],ybarg[bg],FCI[bg,2]) 

abline(v=summary(fit_lme)$coef[1],col="gray",lty=2 )  
abline(h=summary(fit_lme)$coef[1],col="gray",lty=2 )

abline(0,1,col="gray",lty=2)
dev.off() 

WR<-( UCI[,2]-UCI[,1] )/( FCI[,2]-FCI[,1] ) 
sum(WR>1,na.rm=T)
sum(!is.na(WR)) 
mean(WR,na.rm=TRUE) 
tapply(y,g,mean)[ which( WR<1 ) ]
ng[ which( WR<1 ) ] 


#### ---- simulation results
COVER<-WIDTH<-matrix(0,ncol=3,nrow=85) 
nsim<-0
for(sgroup in 1:4)
{
  load(paste0("simres",sgroup,"_alpha",alpha)) 
  COVER[,1]<-COVER[,1]+ucover 
  COVER[,2]<-COVER[,2]+fcover 
  COVER[,3]<-COVER[,3]+bcover  

  WIDTH[,1]<-WIDTH[,1]+uwidth 
  WIDTH[,2]<-WIDTH[,2]+fwidth
  WIDTH[,3]<-WIDTH[,3]+bwidth

  nsim<-nsim+s  
} 
COVER<-COVER/nsim ; WIDTH<-WIDTH/nsim 


#### ---- plots
pdf("radon_sim.pdf",height=3.5,width=8.5,family="Times")
par(mfrow=c(1,2),mar=c(3,3,1,1),mgp=c(1.75,.75,0))

## -- expected widths given theta
plot(ybg,WIDTH[,2]/WIDTH[,1],
     ylim=range(c(WIDTH[,2]/WIDTH[,1],WIDTH[,3]/WIDTH[,1]),na.rm=T),
     pch=16,xlab=expression(theta),ylab="width ratio")
points(ybg,WIDTH[,3]/WIDTH[,1],col="gray",pch=16)
segments(ybg,WIDTH[,2]/WIDTH[,1],ybg,WIDTH[,3]/WIDTH[,1])
abline(h=1,lty=2)
legend(2.15,1.25,legend=c("FAB","Bayes"),pch=c(16,16),col=c("black","gray")  )
abline(v=mean(ybg),lty=2,col="gray")

## -- coverage of Bayes intervals
plot(ybg,COVER[,3],xlab=expression(theta),ylab="posterior interval coverage",pch=16,
     col="gray")
segments(ybg,COVER[,3]+1.96*sqrt(COVER[,3]*(1-COVER[,3])/nsim),
       ybg,COVER[,3]-1.96*sqrt(COVER[,3]*(1-COVER[,3])/nsim),col="gray" )
abline(h=1-alpha,lty=2)
abline(v=mean(ybg),lty=2,col="gray")

dev.off()

sum( WIDTH[,2]>WIDTH[,1],na.rm=TRUE )

apply(WIDTH,2,mean,na.rm=TRUE) 
apply(COVER,2,mean,na.rm=TRUE)

range(COVER[,3],na.rm=T )

