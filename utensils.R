################ Library

library(Matrix)
library(parallel)
library(future.apply)
library(matrixcalc)
library(covKCD)
library(MASS)
library(pracma)
library(expm)
library(RSpectra)
source("lsmr.R")

######## Utensils 

sym.root=function(Sigma){
  
  ###### Compute a symmetric square root of Sigma
  
  Sigma.eig=eigen(Sigma); Q=Sigma.eig$vectors
  return(Q%*%diag(sqrt(Sigma.eig$values))%*%t(Q))
}

sym.inv.root=function(Sigma){
  
  ###### Compute a symmetric inverse square root of Sigma
  
  Sigma.eig=eigen(Sigma); Q=Sigma.eig$vectors
  return(Q%*%diag(sqrt(1/Sigma.eig$values))%*%t(Q))
}

sym.mat=function(V){
  return((V+t(V))/2)
}

det.proj=function(Sigma,Sigma.inv,V){

  p=ncol(Sigma)
  return(V-sum(diag(Sigma.inv%*%V))*Sigma/p)
  
}

det.chol.proj=function(L,L.inv,V){
  
  p=ncol(L)
  return(V-sum(diag(L.inv%*%V))*diag(diag(L))/p)

}

sq.mat=function(mat){
  return(mat%*%t(mat))
}

array.sq=function(dat){
  
  dim=dim(dat)
  return(array(t(apply(dat,MARGIN=1,FUN=sq.mat)),dim=c(dim[1],dim[2],dim[2]))) 
  
}

array.cov=function(dat){
  
  dat.dim=dim(dat); n=dat.dim[1]; p1=dat.dim[2]; p2=dat.dim[3]; p=p1*p2
  
  S=matrix(0,p,p)
  for(i in 1:n){
    S=S+as.numeric(dat[i,,])%*%t(as.numeric(dat[i,,]))/n
  }
  return(S)
}

spd.geodesic=function(Sigma,V){
  
 ###### Compute an exponential map emanating from Sigma in the direction of V under affine-invariant metric
  
  Sigma.root=sym.root(Sigma); Sigma.inv.root=sym.inv.root(Sigma)
  return(Sigma.root%*%as.matrix(expm::expm(Sigma.inv.root%*%V%*%Sigma.inv.root))%*%Sigma.root)
}

chol.geodesic=function(L,V){
  
  ###### Compute an exponential map emanating from L in the direction of V under Cholesky metric
  
  L.diag=diag(L); V.diag=diag(V)
  L.lower=L; L.lower[upper.tri(L.lower,diag=TRUE)]=0
  V.lower=V; V.lower[upper.tri(V.lower,diag=TRUE)]=0
  return(L.lower+V.lower+diag(L.diag*exp(V.diag/L.diag)))
}

mat2array=function(mat,p1){
  
  n=nrow(mat); p=ncol(mat); p2=p/p1
  mat.array=array(0,dim=c(n,p1,p2))
  for(i in 1:n){
    mat.array[i,,]=matrix(mat[i,],p1,p2)
  }
  return(mat.array)
}

pt.mat=function(mat,MARGIN,p1,p2){
  
  ###### Compute the partial trace operators
  ###### Argument 
  ## mat: p1p2 x p1p2 matrix
  ## MARGIN: 1 (first type, p1 x p1), 2 (second type, p2 x p2)
  ## p1: row dimension
  ## p2: column dimension 
  
  p=ncol(mat)
  if(p!=(p1*p2)){
    stop(message("The row and column dimensions are incorrectly specified."))
  }
  if((MARGIN!=1) && (MARGIN!=2)){
    stop(message("MARGIN should be either 1 or 2."))
  }
  
  if(MARGIN==1){
    temp=matrix(0,p1,p1)
    if(isSymmetric(mat)==1){
      for(i in 1:p2){
        temp=temp+mat[((i-1)*p1+1):(i*p1),((i-1)*p1+1):(i*p1)]
      }
    }else{
      mat=(mat+t(mat))/2
      for(i in 1:p2){
        temp=temp+mat[((i-1)*p1+1):(i*p1),((i-1)*p1+1):(i*p1)]
      }
    }
  }else{
    temp=matrix(0,p2,p2)
    if(isSymmetric(mat)==1){
      for(i in 1:p2){
        for(j in i:p2){
          temp[i,j]=temp[j,i]=sum(diag(mat[((i-1)*p1+1):(i*p1),((j-1)*p1+1):(j*p1)]))
        }
      }
    }else{
      mat=(mat+t(mat))/2
      for(i in 1:p2){
        for(j in i:p2){
          temp[i,j]=temp[j,i]=sum(diag(mat[((i-1)*p1+1):(i*p1),((j-1)*p1+1):(j*p1)]))
        }
      }
    }
  }
  return(temp)
}


core.tangent.proj=function(mat,p1,p2){
  
  ###### Compute the orthogonal projection of symmetric matrix onto the tangent space of core covariance manifold
  ###### Argument 
  ## mat: p1p2 x p1p2 symmetric matrix
  ## p1: row dimension
  ## p2: column dimension 

  p=ncol(mat)
  if(isSymmetric(mat)!=1){
    stop(message("A matrix should be symmetric."))
  }
  if(p!=(p1*p2)){
    stop(message("The row and column dimensions are incorrectly specified."))
  }

  pt1=pt.mat(mat,1,p1,p2)
  pt2=pt.mat(mat,2,p1,p2)
  return(mat-kronecker(diag(p2),pt1)/p2-kronecker(pt2,diag(p1))/p1+diag(p)*mean(diag(mat)))   
  
}

core.tangent=function(p1,p2){
  
  ###### Randomly generate a tangent vector of core covariance manifold.
  ###### Argument 
  ## p1: row dimension
  ## p2: column dimension 
  
  p=p1*p2
  Y=matrix(rnorm(p*p),ncol=p); Y=(Y+t(Y))/2
  return(core.tangent.proj(Y,p1,p2))

}

sym2chol=function(V){
  
  ###### Lower-Triangularization of Symmetric matrix 
  
  if(isSymmetric(V)!=1){
    stop(message("The input should be a symmetric matrix"))
  }
  
  V.diag=diag(V)
  V.lower=V; V.lower[upper.tri(V)]=0
  return(V.lower-diag(V.diag)/2)
  
}

mat2lower=function(V){
  
  ##### Lower-Trianugularization of a square matrix
  V.lower=V; V.lower[upper.tri(V)]=0
  return(V.lower)
  
}

J.mat=function(A,p1){
  
  p=nrow(A); p2=p/p1
  r=ncol(A); a=as.numeric(A)
  
  K1=commutation.matrix(p1,p1)
  K2=commutation.matrix(p2,p2)
  J1=NULL
  for(i in 1:r){
    J1=cbind(J1,kronecker(matrix(A[,i],p1,p2),diag(p1)))
  }
  J1=(diag(p1^2)+K1)%*%J1/p-2/(p*p1)*as.numeric(diag(p1))%*%t(a)
  J2=NULL
  for(i in 1:r){
    J2=cbind(J2,kronecker(diag(p2),t(matrix(A[,i],p1,p2))))
  }
  J2=(diag(p2^2)+K2)%*%J2/p-2/(p*p2)*as.numeric(diag(p2))%*%t(a)
  return(rbind(J1,J2,2*t(a)))
}

diff.h=function(K1,K2,U1,U2,root="sym"){
  
  ###### Compute differential of h 
  ###### Argument 
  ## K1: p1 x p1 row covariance 
  ## K2: p2 x p2 column covariance
  ## U1: p1 x p1 symmetric
  ## U2: p2 x p2 symmetric 
  ## root: "sym" (Symmetric) or "chol" (Cholesky) for separable component
  
  if(root=="sym"){
    
    p1=ncol(K1); p2=ncol(K2)
    K1.eig=eigen(K1); K2.eig=eigen(K2)
    Q1=K1.eig$vectors; Q2=K2.eig$vectors
    lambda1=K1.eig$values; lambda2=K2.eig$values 
    lambda1.sqrt=sqrt(lambda1); lambda2.sqrt=sqrt(lambda2)
    Q=kronecker(Q2,Q1)
    Lambda=kronecker(rep(1,p2)%*%t(lambda2.sqrt),rep(1,p1)%*%t(lambda1.sqrt)); Lambda=Lambda+t(Lambda)
    temp=kronecker(diag(lambda2),t(Q1)%*%U1%*%Q1)+kronecker(t(Q2)%*%U2%*%Q2,diag(lambda1))
    temp=temp*(1/Lambda)
    return(Q%*%temp%*%t(Q))
    
    
  }else{
    
    L1=t(chol(K1)); L2=t(chol(K2))
    L1.inv=solve(L1); L2.inv=solve(L2)
    p1=ncol(L1); p2=ncol(L2)
    L=kronecker(L2,L1)
    temp1=sym2chol(L1.inv%*%U1%*%t(L1.inv))
    temp2=sym2chol(L2.inv%*%U2%*%t(L2.inv))
    return(L%*%(kronecker(diag(p2),temp1)+kronecker(temp2,diag(p1))))
    
  }
  
}

diff.g=function(K1,K2,C,U1,U2,W,root="sym"){
  
  ###### Compute differential of g
  ###### Argument 
  ## K1: p1 x p1 row covariance 
  ## K2: p2 x p2 column covariance
  ## C: p x p core covariance (p=p1p2)
  ## U1: p1 x p1 symmetric
  ## U2: p2 x p2 symmetric 
  ## W: tangent vector of core covariance manifold 
  ## root: "sym" (Symmetric) or "chol" (Cholesky) for separable component
  
  if(root=="sym"){
    
    H=diff.h(K1,K2,U1,U2,root=root)
    S1=sym.root(K1); S2=sym.root(K2); S=kronecker(S2,S1)
    return(S%*%W%*%t(S)+H%*%C%*%t(S)+S%*%C%*%t(H))
  
    
  }else{
    
    H=diff.h(K1,K2,U1,U2,root=root)
    L1=t(chol(K1)); L2=t(chol(K2)); L=kronecker(L2,L1)
    return(L%*%W%*%t(L)+H%*%C%*%t(L)+L%*%C%*%t(H))
    
  }
  
}

hess.k=function(K1,K2,C,V,p1,p2){
  
  ###### Compute the Riemannian Hessian operator of 
  ###### Argument 
  ## Sigma: p x p covariance matrix (p=p1p2)
  ## V: p x p symmetric matrix
  ## p1: row dimension
  ## p2: column dimension
  ## sep: whether Sigma is separable (TRUE) or not (FALSE) 
  
  K1.inv.root=sym.inv.root(K1); K1.root=sym.root(K1)
  K2.inv.root=sym.inv.root(K2); K2.root=sym.root(K2)
  K.inv.root=kronecker(K2.inv.root,K1.inv.root)
  V.whit=K.inv.root%*%V%*%t(K.inv.root)
  p=p1*p2

  V1=pt.mat(V.whit,1,p1,p2); V1=K1.root%*%(V1-sum(diag(V.whit))/p1*diag(p1))%*%t(K1.root)/p2
  V2=pt.mat(V.whit,2,p1,p2); V2=K2.root%*%V2%*%t(K2.root)/p1

  f=function(u){
    
    u1=u[1:(choose(p1+1,2)-1)]
    u2=u[(choose(p1+1,2)):(choose(p1+1,2)+choose(p2+1,2)-1)]
  
    U1=matrix(0,p1,p1); U1[lower.tri(U1)]=u1[p1:(choose(p1+1,2)-1)]
    U1=U1+t(U1)
    diag(U1)[1:(p1-1)]=u1[1:(p1-1)]
    diag(U1)[p1]=-sum(u1[1:(p1-1)])
    
    U2=matrix(0,p2,p2); U2[lower.tri(U2,diag=TRUE)]=u2
    U2=U2+t(U2)-diag(diag(U2))
    U2.whit=K2.inv.root%*%U2%*%t(K2.inv.root)
    
    M1=matrix(0,p1,p1); M2=matrix(0,p2,p2)
    
    for(i in 1:p2){
      for(j in 1:p2){
        M1=M1+U2.whit[i,j]*C[((j-1)*p1+1):(j*p1),((i-1)*p1+1):(i*p1)]
        M2[i,j]=sum(diag(C[((i-1)*p1+1):(i*p1),((j-1)*p1+1):(j*p1)]%*%U1))
      }
    }
    U1=K1.root%*%U1%*%t(K1.root)
    R1=U1+K1.root%*%M1%*%t(K1.root)/p2-sum(diag(U2.whit))*K1/p2
    R2=U2+K2.root%*%M2%*%t(K2.root)/p1
    
    return(sqrt(sum((K1.inv.root%*%(R1-V1)%*%t(K1.inv.root))^2))+sqrt(sum((K2.inv.root%*%(R2-V2)%*%t(K2.inv.root))^2)))
  }
 
   u=optim(par=rep(1,(choose(p1+1,2)+choose(p2+1,2)-1)),fn=f,method="BFGS",hessian=TRUE)$par
   u1=u[1:(choose(p1+1,2)-1)]
   u2=u[(choose(p1+1,2)):(choose(p1+1,2)+choose(p2+1,2)-1)]
   U1=matrix(0,p1,p1); U1[lower.tri(U1)]=u1[p1:(choose(p1+1,2)-1)]
   U1=U1+t(U1)
   diag(U1)[1:(p1-1)]=u1[1:(p1-1)]
   diag(U1)[p1]=-sum(u1[1:(p1-1)])
   U1=K1.root%*%U1%*%t(K1.root)
   
   U2=matrix(0,p2,p2); U2[lower.tri(U2,diag=TRUE)]=u2
   U2=U2+t(U2)-diag(diag(U2))
  
   return(list(U1=U1,U2=U2))
  
}

diff.k=function(Sigma,V,p1,p2,sep=FALSE){
  
  ###### Compute differential of k
  ###### Argument 
  ## Sigma: p x p covariance matrix (p=p1p2)
  ## V: p x p symmetric matrix
  ## p1: row dimension
  ## p2: column dimension
  ## sep: whether Sigma is separable (TRUE) or not (FALSE) 

  
  if(sep==TRUE){
    
    Sigma.KCD=covKCD(Sigma,p1,p2)
    K1=Sigma.KCD$K1; K2=Sigma.KCD$K2
    K1.inv.root=sym.inv.root(K1)
    K2.inv.root=sym.inv.root(K2)
    K.inv.root=kronecker(K2.inv.root,K1.inv.root)
    V.whit=K.inv.root%*%V%*%t(K.inv.root)
    
    K1.root=sym.root(K1); K2.root=sym.root(K2)
    K.root=kronecker(K2.root,K1.root)
    p=p1*p2
    V1=pt.mat(V.whit,1,p1,p2)
    V2=pt.mat(V.whit,2,p1,p2)
    return(K.root%*%(kronecker(V2,diag(p1))/p1+kronecker(diag(p2),V1)/p2-diag(p)*sum(diag(V.whit))/p)%*%t(K.root))
    
  }else{
    
    Sigma.KCD=covKCD(Sigma,p1,p2)
    K1=Sigma.KCD$K1; K2=Sigma.KCD$K2
    C=Sigma.KCD$C
    
    U=hess.k(K1,K2,C,V,p1,p2)
    U1=U[[1]]; U2=U[[2]]
    return(kronecker(U2,K1)+kronecker(K2,U1))
  
  }
  
}

diff.c=function(Sigma,V,p1,p2,root="sym",sep=FALSE){
  
  ###### Compute differential of c
  ###### Argument 
  ## Sigma: p x p covariance matrix (p=p1p2)
  ## V: p x p symmetric matrix
  ## p1: row dimension
  ## p2: column dimension
  ## root: "sym" (Symmetric) or "chol" (Cholesky) for separable component 
  ## sep: whether Sigma is separable (TRUE) or not (FALSE) 
  
  
  if(sep==TRUE){
    
    Sigma.KCD=covKCD(Sigma,p1,p2)
    K1=Sigma.KCD$K1; K2=Sigma.KCD$K2
    
    K1.inv.root=sym.inv.root(K1); K1.root=sym.root(K1)
    K2.inv.root=sym.inv.root(K2); K2.root=sym.root(K2)
    K.inv.root=kronecker(K2.inv.root,K1.inv.root)
    V.whit=K.inv.root%*%V%*%t(K.inv.root)
    V1=pt.mat(V.whit,1,p1,p2); V2=pt.mat(V.whit,2,p1,p2)
    U1=K1.root%*%(V1-sum(diag(V.whit))/p1*diag(p1))%*%t(K1.root)/p2
    U2=K2.root%*%V2%*%t(K2.root)/p1
    
    if(root=="sym"){
      
      K1.eig=eigen(K1); K2.eig=eigen(K2)
      Q1=K1.eig$vectors; Q2=K2.eig$vectors
      lambda1=K1.eig$values; lambda2=K2.eig$values 
      lambda1.sqrt=sqrt(lambda1); lambda2.sqrt=sqrt(lambda2)
      Q=kronecker(Q2,Q1)
      Lambda=kronecker(rep(1,p2)%*%t(lambda2.sqrt),rep(1,p1)%*%t(lambda1.sqrt)); Lambda=Lambda+t(Lambda)
      temp=kronecker(diag(lambda2),t(Q1)%*%U1%*%Q1)+kronecker(t(Q2)%*%U2%*%Q2,diag(lambda1))
      temp=temp*(1/Lambda)
      temp=Q%*%kronecker(diag(1/lambda2.sqrt),diag(1/lambda1.sqrt))%*%temp%*%t(Q)
        
      temp=V.whit-temp-t(temp)
      return(temp)
      
    }else{
      
      L1=t(chol(K1)); L1.inv=solve(L1)
      L2=t(chol(K2)); L2.inv=solve(L2)
      L.inv=kronecker(L2.inv,L1.inv)
      
      temp=kronecker(diag(p2),sym2chol(L1.inv%*%U1%*%t(L1.inv)))+kronecker(sym2chol(L2.inv%*%U2%*%t(L2.inv)),diag(p1))
      temp=L.inv%*%V%*%t(L.inv)-temp-t(temp)
      return(temp)
      
    }
   
  }else{
    
    Sigma.KCD=covKCD(Sigma,p1,p2)
    K1=Sigma.KCD$K1; K2=Sigma.KCD$K2
    C=Sigma.KCD$C
    
    U=hess.k(K1,K2,C,V,p1,p2)
    U1=U[[1]]; U2=U[[2]]
    H=diff.h(K1,K2,U1,U2,root=root)
    
    if(root=="sym"){
      H.inv=kronecker(sym.inv.root(K2),sym.inv.root(K1))
    }else{
      H.inv=kronecker(solve(t(chol(K2))),solve(t(chol(K1))))
    }
    temp=H.inv%*%H
    return(H.inv%*%V%*%t(H.inv)-temp%*%C-t(C)%*%t(temp))
    
  }

}

diff.f=function(Sigma,V,p1,p2,root="sym",sep=FALSE){
  
  K=diff.k(Sigma,V,p1,p2,root=root,sep=sep)
  C=diff.c(Sigma,V,p1,p2,root=root,sep=sep)
  return(list(diff.k=K,diff.c=C))
  
}

