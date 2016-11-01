#### ---- check KKT conditions for L1 penalized regression estimates
KKT<-function(beta,Q,l,lambda,tol=1e-6)
{
  ssg<-2*(l - Q%*%beta)/lambda
  ( sum( (ssg-sign(beta))[abs(beta)>tol ]^2 ) < tol ) &
  ( all( abs(ssg[abs(beta)<tol]) <1 ) )
}



#### ---- objective function 
l1objective<-function(beta,Q,l,lambda)
{
  t(beta)%*%( Q%*%beta -  2*l) + lambda*sum(abs(beta))
}



#### ---- LQA 
lqa_lm<-function(y,X,lambda,
                 beta=NULL,tol=1e-6,itmin=1,itmax=Inf,epsilon=1e-12)
{
  Q<-crossprod(X)
  l<-crossprod(X,y)
  if(is.null(beta) & nrow(X)>=ncol(X)){ beta<-chol2inv(chol(Q))%*%l }
  if(is.null(beta) & nrow(X)<ncol(X))
  { 
    beta<-apply(X,2,function(x){ lm(y~ -1+ x)$coef } )
  } 

  obj<-l1objective(beta,Q,l,lambda) 
  while( ( length(obj)<=itmin ) |
         ( !KKT(beta,Q,l,lambda,tol=tol) & length(obj)<=itmax) )
  {
    W<-.5*lambda*diag( c(abs(beta)+epsilon)^(-1) ) 
    R<-chol(  Q + W  )
    beta<-backsolve( R , forwardsolve(t(R), l )  ) 

    obj<-c(obj,l1objective(beta,Q,l,lambda)) 
  }
  list(beta=beta,obj=obj)
}



#### ---- HPP
hpp_lm<-function(y,X,lambda,
                 beta=NULL,tol=1e-6,itmin=1,itmax=Inf)
{
  Q<-crossprod(X)
  l<-crossprod(X,y) 
  if(is.null(beta) & nrow(X)>=ncol(X)){ beta<-chol2inv(chol(Q))%*%l }
  if(is.null(beta) & nrow(X)<ncol(X))
  { 
    beta<-apply(X,2,function(x){ lm(y~ -1+ x)$coef } )
  } 

  v<-beta/sqrt(abs(beta)) ; u<-beta/v
  S<-lambda*diag(length(l))/2 

  obj<-l1objective(beta,Q,l,lambda)
  while( ( length(obj)<=itmin ) | 
         ( !KKT(u*v,Q,l,lambda,tol=tol) & length(obj)<=itmax) )
  {
    R<-chol( ( Q * v%*%t(v) ) + S )
    u<-backsolve( R , forwardsolve(t(R), l*v )  )

    R<-chol( ( Q * u%*%t(u) ) + S )
    v<-backsolve( R , forwardsolve(t(R), l*u )  )

    obj<-c(obj,l1objective(u*v,Q,l,lambda)) 
  }
  list(beta=u*v,obj=obj) 
}


