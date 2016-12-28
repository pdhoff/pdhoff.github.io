#### ----
wfabz<-function(theta,mu,t2,alpha=.05) 
{
  igfun<-function(x,alpha)
  {
    gwmx <-function(w){ qnorm(alpha*w) - qnorm(alpha*(1-w)) - x }
    uniroot(gwmx,interval=c(0,1),maxiter=2000,tol=.Machine$double.eps^0.5)$root
  }
  igfun( 2*(theta-mu)/t2,alpha)
}



#### ----
fabz_ci<-function(y,mu,t2,alpha=.05) 
{ 
  tmp<-function(th){ (y+qnorm(1-alpha + pnorm(y-th)) + 2*(mu-th)/t2 ) - th }  
  if(tmp( y+qnorm(1-alpha)+1e-12) < 0 ){ thU<-y+qnorm(1-alpha) } else
  { thU<-uniroot(tmp,interval=c(y+qnorm(1-alpha)+1e-12,y+10))$root }

  tmp<-function(th){ (y+qnorm(alpha - pnorm(th-y)) + 2*(mu-th)/t2 ) - th }   
  if(tmp( y+qnorm(alpha)-1e-12 ) > 0 ){ thL<-y+qnorm(alpha) }  else 
   { thL<-uniroot(tmp,interval=c(y-10,y+qnorm(alpha)-1e-12))$root }

  c(thL,thU) 
}



#### ---- function for computing expected width at mu=0
fabz_ew<-function(theta,t2,alpha=.05)
{ 

  hx<-function(xs,alpha=.05)
  {
    hxs<-NULL
    for(x in xs)
    {
    ih_mx <-function(w){ qnorm(alpha*w) - qnorm(alpha*(1-w)) - x }
    hxs<-c(hxs,uniroot(ih_mx,interval=c(0,1),maxiter=2000,
               tol = .Machine$double.eps^0.5)$root  )
    }
    hxs
  }

  H1<-function(x,theta,t2,alpha=.05)
  {
    h1<-NULL
    for(xi in x)
    {
      h1<-c(h1,.5*t2*pnorm( theta - .5*t2*xi + qnorm( 1-alpha*hx(xi,alpha) ) ) )
    }
    h1
  } 
 
  H2<-function(x,theta,t2,alpha=.05)
  {
    h2<-NULL
    for(xi in x)
    {
      h2<-c(h2,.5*t2*pnorm( theta - .5*t2*xi - qnorm( 1-alpha*hx(-xi,alpha) ) ) )
    }
    h2
  }

  H<-function(x){ H1(x,theta,t2,alpha) - H2(x,theta,t2,alpha) }
  integrate(H,lower= -Inf,upper=Inf,
            subdivisions = 800L)$value
}



#### ---- functions for t-interval
paccept<-function(w,theta,mu,t2,nu0,s20,n,alpha)
{
  a<-nu0/2 ; b<-nu0*s20/2

  f<-function(s2)
  {
    c<-sqrt(s2/n)/sqrt(s2/n+t2)
    ncp<- c*(mu - theta)/(sqrt(s2/n))
    suppressWarnings(
    pacc<- pt( c*qt(1-alpha*(1-w),n-1),n-1,ncp ) -
           pt( c*qt(alpha*w,n-1),n-1,ncp )
                    )

    a<-nu0/2 ; b<-nu0*s20/2
    ds2<-exp( a*log(b) - lgamma(a) - (a+1)*log(s2) - b/s2 )

    pacc*ds2
  }

  int<-try(integrate(f,lower=0,upper=Inf)$val,silent=TRUE)  
  if(is.character(int)){int<-1}  
  int    
}


wfabt<-function(theta,mu,t2,nu0,s20,n,alpha)
{
  optimize( paccept, lower=0,upper=1,
            theta=theta,mu=mu,t2=t2,nu0=nu0,s20=s20,n=n,alpha=alpha)$min
}


fabt_ci<-function(ybar,s,n,mu,t2,nu0,s20,alpha=.05)
{

  ubroot<-function(theta)
  {   
    w<-wfabt(theta,mu,t2,nu0,s20,n,alpha)
    ybar + s*qt(1-alpha*w,n-1)/sqrt(n) - theta
  }  
  a<-b<-ybar + .99*(s/sqrt(n))*qt(1-alpha,n-1)  
  while(ubroot(b)>0){ b<- b + (s/sqrt(n))*qnorm(1-alpha)*n/(n+4) } 
  thetaU<-uniroot(ubroot,c(a,b))$root

  lbroot<-function(theta)
  { 
    w<-wfabt(theta,mu,t2,nu0,s20,n,alpha)
    ybar + s*qt(alpha*(1-w),n-1)/sqrt(n) - theta 
  } 
  a<-b<-ybar + .99*(s/sqrt(n))*qt(alpha,n-1) 
  while(lbroot(a)<0){ a<- a + (s/sqrt(n))*qnorm(alpha)*n/(n+4) } 
  thetaL<-uniroot(lbroot,c(a,b))$root  

  c(thetaL,thetaU)
}



fabt_ew<-function(theta,s2,n,mu,t2,nu0,s20,alpha=.05) 
{
  dys2fabtw<-function(ys)
  {
    diff( fabt_ci(ys[1],ys[2],n,mu,t2,nu0,s20,alpha) )* 
    dnorm(ys[1],theta,sqrt(s2/n))* 
    dchisq((n-1)*ys[2]^2/s2,n-1)*(n-1)/s2 
  }

  adaptIntegrate(dys2fabtw,lowerLimit=c(-Inf,0),upperLimit=c(Inf,Inf))

adaptIntegrate(dys2fabtw,lowerLimit=c(0,0),upperLimit=c(3,3))

}


fabt_ew<-function(theta,s2,n,mu,t2,nu0,s20,alpha=.05,nmc=1000) 
{
   ybar<-rnorm(nmc,theta,sqrt(s2/n)) 
   s<-sqrt( s2*rchisq(nmc,n-1)/(n-1) ) 
   ciw<-NULL
   for(i in 1:nmc)    
   {
     ciw<-c(ciw,diff(fabt_ci(ybar[i],s[i],n,mu,t2,nu0,s20,alpha)))
   }
   mean(ciw)  
}



#### ---- compute ebayes estimate of (mu,t2,nu0,s20) 
ebayes_fabt<-function(y,g)
{
  n<-c(table(g))
  ss<-(n-1)*c(tapply(y,g,var))
  ybar<-tapply(y,g,mean)

  n<-n[!is.na(ss)]
  ss<-ss[!is.na(ss)]
  ybar<-ybar[!is.na(ss)]

  ## -- get mmle for prior over s2
  mllab<-function(lab)
  {
    a<-exp(lab[1]) ; b<-exp(lab[2])
    ap<-a+(n-1)/2
    bp<-b+ss/2
    -sum( ( a*log(b)-lgamma(a) ) - ( ap*log(bp)-lgamma(ap) ) )
  }
  l<-(n-1)/ss ; lab0<-log(c( mean(l)^2/var(l) ,mean(l)/var(l) ))
  ab<-exp(optim(lab0, mllab)$par)

  ## -- get posterior modes for each s2 
  a<-ab[1] ; b<-ab[2]
  ap<-a+(n-1)/2
  bp<-b+ss/2
  s2<-bp/(ap+1)

  ## -- mmle of mu, t2 with plugin s2
  mllmut2<-function(mut2)
  {
    -sum(dnorm(ybar,mut2[1],sqrt(s2/n+exp(mut2[2])),log=TRUE) )
  }
  mut20<-c(mean(ybar),log(var(ybar)) )
  mut2<-optim(mut20,mllmut2)$par
  mut2[2]<-exp(mut2[2])

  list(mu=mut2[1],t2=mut2[2],nu0=ab[1]*2,s20=ab[2]/ab[1])
}











