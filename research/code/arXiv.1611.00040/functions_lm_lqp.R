#### ---- objective function 
lqobjective<-function(beta,Q,l,lambda,q)
{
  t(beta)%*%( Q%*%beta -  2*l) + lambda*sum(abs(beta)^q)
}



#### ---- LQA 
lqa_lm<-function(y,X,lambda,q=1,
                 beta=NULL,tol=1e-6,itmin=1,itmax=Inf,epsilon=1e-12)
{
  Q<-crossprod(X)
  l<-crossprod(X,y) 
  if(is.null(beta) & nrow(X)>=ncol(X)){ beta<-chol2inv(chol(Q))%*%l }
  if(is.null(beta) & nrow(X)<ncol(X))
  { 
    beta<-apply(X,2,function(x){ lm(y~ -1+ x)$coef } )
  }

  obj<-lqobjective(beta,Q,l,lambda,q) ; del<-1
  while( ( length(obj)<=itmin ) | ( del>tol & length(obj)<=itmax) )
  {
    W<-.5*lambda*q*diag( c(abs(beta)^(q-1))/c(abs(beta)+epsilon) )
    R<-chol(  Q + W  )
    beta<-backsolve( R , forwardsolve(t(R), l )  )

    obj<-c(obj,lqobjective(beta,Q,l,lambda,q))
    del<-( obj[length(obj)-1] - obj[length(obj)] )/abs( obj[length(obj)] )

  }
  list(beta=beta,obj=obj)
}



#### ---- HPP
hpp_lm<-function(y,X,lambda,K=2,
                 beta=NULL,tol=1e-6,itmin=1,itmax=Inf)
{
  Q<-crossprod(X)
  l<-crossprod(X,y) 
  if(is.null(beta) & nrow(X)>=ncol(X)){ beta<-chol2inv(chol(Q))%*%l }
  if(is.null(beta) & nrow(X)<ncol(X))
  { 
    beta<-apply(X,2,function(x){ lm(y~ -1+ x)$coef } )
  }

  U<-matrix(abs(beta)^(1/K),ncol(X),K)
  U[,1]<- U[,1]*sign(beta)
  lIpK<-diag(length(l))*lambda/K

  obj<-lqobjective(beta,Q,l,lambda,2/K) ; del<-1
  while( ( length(obj)<=itmin ) | ( del>tol & length(obj)<=itmax) ) 
  {

    for(k in 1:K)
    {
      v<-apply(U[,-k,drop=FALSE],1,prod)
      R<-chol( ( Q * v%*%t(v) ) + lIpK )
      U[,k]<-backsolve( R , forwardsolve(t(R), l*v )  )
    } 

    beta<-apply(U,1,prod)
    obj<-c(obj,lqobjective(beta,Q,l,lambda,2/K))
    del<-( obj[length(obj)-1] - obj[length(obj)] )/abs( obj[length(obj)] )

  }
  list(beta=beta,obj=obj)
}



