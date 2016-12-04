#### ---- objective function
lqobjective_glm<-function(beta,y,X,lambda,q=1,family="binomial")
{
  eta<-X%*%beta    
  if(family=="binomial"){ A<-A_bin } 
  if(family=="poisson"){ A<-A_pois }
  -2*sum( eta*y- A(eta) ) + lambda*sum( abs(beta)^q )
}



#### ---- LQA 
lqa_glm<-function(y,X,lambda,q=1,family="binomial",add_int=FALSE, 
                  beta=NULL,tol=1e-6,itmin=1,itmax=Inf,iitmax=Inf,epsilon=1e-12)  
{ 
  if(family=="binomial"){ dA<-dA_bin ; ddA<-ddA_bin }
  if(family=="poisson"){ dA<-dA_pois ; ddA<-ddA_pois } 

  if(add_int)
  { 
    X<-cbind(1,X) 
    if(is.null(beta))
    {
      beta<-apply(X,2,function(x){ glm(y~x,family=family)$coef } )
      beta<-c(mean(beta[1,]),beta[2,] )
    }    
    ridge<-lambda*c(0,rep(1,ncol(X)-1 ))
  }

  if(!add_int)
  { 
    if(is.null(beta))
    {
      beta<-apply(X,2,function(x){ glm(y~ -1+x,family=family)$coef } )
    }
    ridge<-lambda*rep(1,ncol(X))
  }


  del<-1 ; obj<-lqobjective_glm(beta,y,X,lambda,q) 
  while( (length(obj)<=itmin) | ( del>tol & length(obj)<=itmax) )
  {
    W<-.5*q*diag( ridge*( c( abs(beta)^(q-1) )/c(abs(beta)+epsilon)  )  ) 
    beta<-irls(y,X,beta=beta,family=family,ridge=W,itmax=iitmax)
    obj<-c(obj, lqobjective_glm(beta,y,X,lambda,q) )  
    del<-( obj[length(obj)-1] - obj[length(obj)] )/abs( obj[length(obj)] )
  }
  list(beta=beta,obj=obj) 
}



#### ---- fractionally penalized glm using HPP 
hpp_glm<-function(y,X,lambda,K=2,family="binomial",add_int=FALSE, 
                  beta=NULL,tol=1e-6,itmin=1,itmax=Inf,iitmax=Inf) 
{
  if(family=="binomial"){ dA<-dA_bin ; ddA<-ddA_bin }
  if(family=="poisson"){ dA<-dA_pois ; ddA<-ddA_pois }

  if(add_int)
  {
    X<-cbind(1,X)
    if(is.null(beta))
    {
      beta<-apply(X,2,function(x){ glm(y~x,family=family)$coef } )
      beta<-c(mean(beta[1,]),beta[2,] )
    }
    ridge<-diag(lambda*c(0,rep(1,ncol(X)-1 )))/K
  }

  if(!add_int)
  {
    if(is.null(beta))
    {
      beta<-apply(X,2,function(x){ glm(y~ -1+x,family=family)$coef } )
    }
    ridge<-diag(lambda*rep(1,ncol(X)))/K 
  }

  U<-matrix(abs(beta)^(1/K),ncol(X),K)
  U[,1]<- U[,1]*sign(beta)

  del<-1 ; obj<-lqobjective_glm(beta,y,X,lambda,2/K) 
  while( (length(obj)<=itmin) | ( del>tol & length(obj)<=itmax) )
  {
    for(k in 1:K)
    {
      v<-apply(U[,-k,drop=FALSE],1,prod)
      Xv<-sweep(X,2,v,"*") 
      U[,k]<-irls(y,Xv,beta=U[,k],family=family,ridge=ridge,itmax=iitmax)
    }
    beta<-apply(U,1,prod) 
    obj<-c(obj, lqobjective_glm(beta,y,X,lambda,2/K) )    
    del<-( obj[length(obj)-1] - obj[length(obj)] )/abs( obj[length(obj)] )
  }
  list(beta=beta,obj=obj)
}



#### ---- 
irls<-function(y,X,family="binomial",ridge=diag(0,ncol(X)),
               tol=1e-6,beta=beta,itmax=Inf)
{
  if(family=="binomial"){ dA<-dA_bin ; ddA<-ddA_bin }
  if(family=="poisson"){ dA<-dA_pois ; ddA<-ddA_pois }

  del<-1 ; it<-0
  while( del>tol & it<itmax  )
  {
    eta<-c(X%*%beta)
    dl<- apply( X*(y - dA(eta)) ,2,sum)  - ridge%*%beta
    ddl<- -crossprod(X, X*ddA(eta))    - ridge

    R<-chol( -ddl )
    beta_new<-beta + backsolve( R , forwardsolve(t(R), dl )  )

    del<-mean( (beta-beta_new)^2 )/mean(beta^2)
    beta<-beta_new 
    it<-it+1 
  }
  beta
}



#### ---- link function and derivatives 
A_pois<-function(eta){ exp(eta) }
dA_pois<-function(eta){ exp(eta) }
ddA_pois<-function(eta){ exp(eta) }

A_bin<-function(eta){ log(1+exp(eta))  }
dA_bin<-function(eta){ 1/(1+exp(-eta))  }
ddA_bin<-function(eta){ 1/( (1+exp(-eta))*(1+exp(eta)) ) }



