########### Libraries

library(RSpectra)
library(covKCD)
library(pracma)
library(matrixcalc)
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

vec2array=function(v,p1,p2,r){

  p=p1*p2
  A=array(0,dim=c(r,p1,p2))
  for(i in 1:r){
    A[i,,]=matrix(v[((i-1)*p+1):(i*p)],nrow=p1,ncol=p2)
  }
  
  return(A)

}

vec2mat=function(v,r){
  
  p=length(v)/r 
  
}

array2mat=function(A){
  
  A.dim=dim(A); r=A.dim[1]; p1=A.dim[2]; p2=A.dim[3]; p=p1*p2
  
  
}

core.tangent=function(A){
  
  A.dim=dim(A); r=A.dim[1]; p1=A.dim[2]; p2=A.dim[3]; p=p1*p2
  
  commute1=commutation.matrix(p1,p1)
  commute2=commutation.matrix(p2,p2)

  a=NULL
  for(i in 1:r){
    a=c(a,as.numeric(A[i,,]))
  }
  
  A.R=matrix(aperm(A,perm=c(2,3,1)),nrow=p1,ncol=p2*r)
  A.C=matrix(aperm(A,perm=c(3,2,1)),nrow=p2,ncol=p1*r)
  
  J1=kronecker(A.R,diag(p1)); J1=(diag(p1^2)+commute1)%*%J1/p-2/(p*p1)*as.numeric(diag(p1))%*%t(a)
  J2=NULL
  for(i in 1:r){
    J2=cbind(J2,kronecker(diag(p2),t(A[i,,])))
  }
  
  J2=(diag(p2^2)+commute2)%*%J2/p-2/(p*p2)*as.numeric(diag(p2))%*%t(a)
  J3=2*t(a)

  return(null(rbind(J1,J2,J3)))
}

A=matrix(rnorm(40),ncol=4)
S=A%*%t(A); S.kcd=covKCD(S,5,2)
K1.inv=sym.inv.root(S.kcd$K1)
K2.inv=sym.inv.root(S.kcd$K2)
A.array=array(0,dim=c(4,5,2))
for(i in 1:4){
  A.array[i,,]=K1.inv%*%matrix(A[,i],5,2)%*%K2.inv
}
S=matrix(0,10,4)
for(i in 1:4){
  S[,i]=as.numeric(A.array[i,,])
}
N=core.tangent(A.array)
b=matrix(rnorm(40))
V=N%*%t(N)%*%b
V=matrix(V,ncol=4)
M=S+V
M=M%*%t(M)
M.kcd=covKCD(M,5,2)
covKCD(S%*%t(S),5,2)
covKCD(M%*%t(M),5,2)

core.hessian=function(A){
  
}

core.retract=function(A,euclid.deriv,metric){
  
  A.dim=dim(A); r=dim(A)[1]; p1=dim(A)[2]; p2=dim(A)[3]
  
  A.mat=array2mat(A)
  
  v=-ginv(M%*%t(N))%*%N%*%t(N)%*%euclid.deriv

  V.mat=vec2mat(v,r)
  M=A.mat+V.mat; M=M%*%t(M)
  M.kcd=covKCD(M,p1,p2)
  
  if(metric=="AI"){
    K1.inv=sym.inv.root(M.kcd$K1); K2.inv=sym.inv.root(M.kcd$K2)
  }else{
    K1.chol=t(chol(M.kcd$K1)); K2.chol=t(chol(M.kcd$K2))
    K1.inv=backsolve(K1.chol,diag(ncol(p1)),upper.tri=FALSE)
    K2.inv=backsolve(K2.chol,diag(ncol(p2)),upper.tri=FALSE)
  }
  
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
    R=t(u)%*%S.whit%*%u
    return(sum(diag(S.whit))/lambda-(1-lambda)/lambda*sum((diag(R)*d^2)/((1-lambda)*d^2+lambda))
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

optim.A=function(S,K1.inv,K2.inv,A,nu,lambda){
  
  p=nrow(A); r=ncol(A)
  p1=ncol(K1.inv); p2=ncol(K2.inv)
  
  K.inv=kronecker(K2.inv,K1.inv)/nu
  S.whit=K.inv%*%S%*%t(K.inv)
  
  A.svd=svd(A)
  u=A.svd$u; d=A.svd$d
  A.sq_inv=diag(p)/lambda-(1-lambda)/lambda*u%*%diag(d^2/((1-lambda)*d^2+lambda))%*%t(u)
  A.euclid=-2*(1-lambda)*A.sq_inv%*%S.whit%*%t(A.sq_inv)%*%A+2*(1-lambda)%*%A.sq_inv%*%A
  
  anc.mat1=A.sq_inv%*%S.whit%*%A.sq_inv
  anc.mat2=t(A)%*%A.sq_inv%*%A
  anc.mat3=A.sq_inv%*%A
  
  N=core.tangent(A)
  A.hess=2*(1-lambda)*kronecker(diag(r),A.sq_inv)-2*(1-lambda)*kronecker(diag(r),anc.mat1)
  -2*(1-lambda)^2*kronecker(anc.mat2,A.sq_inv)
  
  
  A.hess=A.hess%*%N
  
  
  
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

