#### ---- 
lla_lm<-function(y,X,lambda,q=1,
                  beta=NULL,tol=1e-6,itmin=1,itmax=Inf,iitmax)
{
  Q<-crossprod(X)
  l<-crossprod(X,y)
  if(is.null(beta) & nrow(X)>=ncol(X)){ beta<-chol2inv(chol(Q))%*%l }
  if(is.null(beta) & nrow(X)<ncol(X))
  {
    beta<-apply(X,2,function(x){ lm(y~ -1+ x)$coef } )
  }

  nzbeta<-1:length(beta) 
  obj<-lqobjective(beta,Q,l,lambda,q) ; del<-1
  while( ( length(obj)<=itmin ) | ( del>tol & length(obj)<=itmax) )
  {
    w<- q*(abs(beta[nzbeta])^(q-1))
    pfit<-penalized::penalized(y,X[,nzbeta],unpenalized=~0,lambda1=.5*lambda*w,
                              trace=FALSE,maxiter=iitmax)
    beta[nzbeta]<-penalized::coefficients(pfit,"all")
    nzbeta<-which(beta!=0)
 
    obj<-c(obj,lqobjective(beta,Q,l,lambda,q))
    del<-del<-( obj[length(obj)-1] - obj[length(obj)] )/abs( obj[length(obj)] )
  }
  list(beta=beta,obj=obj)
}



#### ---- 
lla_glm<-function(y,X,lambda,q=1,family="binomial",add_int=FALSE,
                  beta=NULL,tol=1e-6,itmin=1,itmax=Inf,iitmax=Inf)
{

 if(is.null(beta))
 {
   beta<-apply(X,2,function(x){ glm(y~ -1+x,family=family)$coef } )
 }

  nzbeta<-1:length(beta)
  del<-1 ; obj<-lqobjective_glm(beta,y,X,lambda,q)
  while((length(obj)<=itmin) | ( del>tol & length(obj)<=itmax) )
  {
    w<- q*(abs(beta[nzbeta])^(q-1)) 
    pfit<-penalized::penalized(y,X[,nzbeta],unpenalized=~0,lambda1=.5*lambda*w,
                              trace=FALSE,maxiter=iitmax)
    beta[nzbeta]<-penalized::coefficients(pfit,"all") 
    nzbeta<-which(beta!=0) 

    obj<-c(obj, lqobjective_glm(beta,y,X,lambda,q) )
    del<-( obj[length(obj)-1] - obj[length(obj)] )/abs( obj[length(obj)] )
  }
  list(beta=beta,obj=obj)
}


