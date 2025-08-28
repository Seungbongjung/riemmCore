########### Libraries

library(RSpectra)
library(covKCD)
library(pracma)
library(expm)
library(MASS)

########### Utensils

array.cov=function(dat){
  
  dat.dim=dim(dat); n=dat.dim[1]; p1=dat.dim[2]; p2=dat.dim[3]; p=p1*p2
  
  S=matrix(0,p,p)
  for(i in 1:n){
    S=S+as.numeric(dat[i,,])%*%t(as.numeric(dat[i,,]))/n
  }
  return(S)
}

skew=function(A){
  return((A-t(A))/2)
}

sym=function(A){
  return((A+t(A))/2)
}

sym.root=function(Sigma){
  
  Sigma.eig=eigen(Sigma)
  Q=Sigma.eig$vectors; D=Sigma.eig$values
  return(Q%*%diag(sqrt(D))%*%t(Q))
  
}

sym.inv.root=function(Sigma){
  
  Sigma.eig=eigen(Sigma)
  Q=Sigma.eig$vectors; D=Sigma.eig$values
  return(Q%*%diag(1/sqrt(D))%*%t(Q))
  
}

core.geod=function(A,V){
  
}

ai.geod=function(Omega,V){
  
  Omega.root=sym.root(Omega) 
  Omega.inv.root=sym.inv.root(Omega)
  return(Omega.root%*%as.matrix(expm::expm(Omega.inv.root%*%V%*%Omega.inv.root))%*%Omega.root)
  
}

lc.geod=function(L,V){
  
  L.diag=diag(L)
  L.lower=L-diag(L.diag)
  V.diag=diag(V)
  V.lower=V-diag(V.diag)
  return(L.lower+V.lower+diag(L.diag)%*%diag(exp(V.diag/L.diag)))
  
}

optim.lambda=function(S,K1.inv,K2.inv,nu,A){
  
  K.inv=kronecker(K2.inv,K1.inv)/nu
  S.whit=K.inv%*%S%*%t(K.inv)
  
  p=ncol(S); r=ncol(A)
  
  A.svd=svd(A)
  d=A.svd$d; u=A.svd$u
  
  obj=function(lambda){
    A.sq_inv=(1-lambda)/lambda*u%*%diag(d^2/((1-lambda)*d^2+lambda))%*%t(u)
    return(sum(diag(S.whit))/lambda-sum(diag(S.whit%*%A.sq_inv))
           +(p-r)*log(lambda)+sum(log((1-lambda)*d^2+lambda)))
  }
  
  lambda.optim=optimize(obj,interval=c(1e-3,1))
  
  return(list(lambda=lambda.optim$optimum,optim=(lambda.optim$value+2*p*log(nu))))
  
  
}

optim.K1=function(S,K1,K1.inv,K2.inv,lambda,A,metric){
  
  p1=ncol(K1); p=ncol(S); r=ncol(A)
  
  K.inv=kronecker(K2.inv,K1.inv)
  S.whit=K.inv%*%S%*%t(K.inv)
  
  A.svd=svd(A)
  u=A.svd$u; d=A.svd$d
  A.sq_inv=(1-lambda)/lambda*u%*%diag(d^2/((1-lambda)*d^2+lambda))%*%t(u)
  
  if(metric=="AI"){
    
    K1.euclid=matrix(0,p1,p1)
    
    K1.hess=matrix(0,p1^2,choose(p1+1,2))
    
    
  }else{
    
  }
  
  return(K1)
}

optim.K2=function(S,K1.inv,K2,K2.inv,lambda,A,metric){
  
  p2=ncol(K2); p=ncol(S); r=ncol(A)
  
  K.inv=kronecker(K2.inv,K1.inv)
  S.whit=K.inv%*%S%*%t(K.inv)
  
  A.svd=svd(A)
  u=A.svd$u; d=A.svd$d
  A.sq_inv=(1-lambda)/lambda*u%*%diag(d^2/((1-lambda)*d^2+lambda))%*%t(u)
  
  if(metric=="AI"){
    
  }else{
    
  }
  
  return(K2)
  
}

optim.nu=function(S,K1.inv,K2.inv,lambda,A){
  
  p=ncol(S)
  
  K.inv=kronecker(K2.inv,K1.inv)
  S.whit=K.inv%*%S%*%t(K.inv)
  
  A.svd=svd(A)
  u=A.svd$u; d=A.svd$d
  A.sq_inv=(1-lambda)/lambda*u%*%diag(d^2/((1-lambda)*d^2+lambda))%*%t(u)
  
  return(sqrt(sum(diag(S.whit)/lambda-diag(S.whit%*%A.sq_inv))/p))
  
}

optim.A=function(S,K1.inv,K2.inv,nu,lambda){
  
  p=nrow(A); r=ncol(A)
  
  K.inv=kronecker(K2.inv,K1.inv)/nu
  S.whit=K.inv%*%S%*%t(K.inv)
  
  A.euclid=matrix(0,p,r)
  
  for(i in 1:p){
    for(j in 1:r){
      basis=matirx(0,p,r)
      basis[i,j]=1
      
    }
  }
  
  A.hess=matrix(0,p*r,p*r)
  
  
  
}


obj=function(S,K1.inv,K2.inv,nu,lambda,A){
  
  p=ncol(S); r=ncol(A)
  
  K.inv=kronecker(K2.inv,K1.inv)/nu
  S.whit=K.inv%*%S%*%t(K.inv)
  
  A.svd=svd(A)
  u=A.svd$u; d=A.svd$d
  A.sq_inv=(1-lambda)/lambda*u%*%diag(d^2/((1-lambda)*d^2+lambda))%*%t(u)
  
  lik=sum(diag(S.whit))/(lambda)-sum(S.whit*A.sq_inv)
  +2*p*log(nu)+sum(log((1-lambda)*d^2+lambda))+(p-r)*log(lambda)
  return(lik)
  
}

optim.lik=function(dat,r,metric="AI",maxiter=100,eps=5e-04,center=TRUE){
  
  ##### Arguments
  ### dat : data array of dimension (n x p1 x p2) (sample size x row dimension x column dimension)
  ### r : partial isotropy rank 
  ### metric : metric on SPD cone, should be either "AI" (affine-invariant) or "LC" (Log-Cholesky)
  ### maxiter : maximum number of iteration 
  ### eps : tolerance parameter for relative convergence criterion
  ### center : whether center the data (TRUE or FALSE)
  
  n=dim(dat)[1]; p1=dim(dat)[2]; p2=dim(dat)[3]
  
  if(r<=(p1/p2+p2/p1) && n<=(p1/p2+p2/p1)){
    stop("Both r and n should be larger than r>p1/p2+p2/p1.")
  }
  if(metric!="AI" && metric!="LC"){
    stop("Metric should either AI or LC.")
  }
  
  ##### Initialization 
  
  if(center==TRUE){
    dat=apply(dat,MARGIN=c(2,3),FUN=scale,scale=FALSE)   
  }
  S=array.cov(dat)
  S.kcd=covKCD(S,p1,p2)
  
  K1=S.kcd$K1; K2=S.kcd$K2; C=S.kcd$C 
  C.eig=eigs_sym(C,r) 
  lambda=(p-sum(C.eig$values))/(p-r)
  C.r=C.eig$vectors%*%diag(C.eig$values)%*%t(C.eig$vectors)
  C.r=covKCD(C.r,p1,p2)$C
  
  C.r.eig=eigs_sym(C.r,r)
  
  if(metric=="AI"){
    
    K1.eig=eigen(K1); K2.eig=eigen(K2)
    K1=K1.eig$vectors%*%diag(sqrt(K1.eig$values))%*%t(K1.eig$vectors)
    K2=K2.eig$vectors%*%diag(sqrt(K2.eig$values))%*%t(K2.eig$vectors)
    
    nu1=exp(mean(log(sqrt(K1.eig$values))))
    nu2=exp(mean(log(sqrt(K2.eig$values))))
    K1=K1/nu1; K2=K2/nu2
    nu=nu1*nu2
    
  }else{
    
    
    K1.root=sym.root(K1)
    K2.root=sym.root(K2)
    
    K1=t(chol(K1)); K2=t(chol(K2))
    nu1=exp(mean(log(diag(K1)))); nu2=exp(mean(log(diag(K2))))
    K1=K1/nu1; K2=K2/nu2
    nu=nu1*nu2
    
    O2
    O1
    O=kronecker(O2,O1)
    
    
    
  }
  
  
  return(list(K1=K1,K2=K2,nu=nu,A=A,lambda=lambda))
  
  
  
}

