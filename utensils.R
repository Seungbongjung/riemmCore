################ Library

library(Matrix)
library(parallel)
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

chol.exp=function(L,V){
  
  ###### Compute an exponential map emanating from L in the direction of V under Cholesky metric
  
  L.diag=diag(L); V.diag=diag(V)
  L.lower=L; L.lower[upper.tri(L.lower,diag=TRUE)]=0
  V.lower=V; V.lower[upper.tri(V.lower,diag=TRUE)]=0
  return(L.lower+V.lower+L.diag%*%exp(V.diag/L.diag))
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

det.proj=function(Sigma,V,metric="AI"){
  
  ###### Compute the orthogonal projection of symmetric matrix onto the tangent space of submanifold of PD cone with unit determinant
  ###### Argument
  ## Sigma: a covariance matrix
  ## V: a asymmetric matrix (tangent vector)
  ## metric: should be either affine-invariant (AI) or log-Cholesky (LC)
  
  if(metric!="AI" && metric!="LC"){
    stop(message("metric should be either affine-invariant (AI) or log-Cholesky (LC)."))
  }
  
  p=ncol(Sigma)
  if(metric=="AI"){
    return(V-sum(diag(solve(Sigma)%*%V))*Sigma/p) 
  }else{
    return()
  }
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
  return(mat-kronecker(diag(p2),pt1)/p2-kronecker(pt2,diag(p1))/p1+diag(p)*sum(diag(mat))/p)   
  
}

core.tangent.proj.ai=function(V,p1,p2){
  
}

core.tangent.proj.lc=function(V,p1,p2){
  
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

sym.basis=function(p,Sigma=0,unit.det=FALSE){
  
  ####### Create a basis 
  ####### Argument 
  
  
}

hess.k=function(Sigma,p1,p2){
  
}

diff.k=function(Sigma,K1=0,K2=0,p1,p2,sep=FALSE){
  
}

diff.c=function(Sigma,p1,p2,root="sym",sep=FALSE){
  
}

diff.f=function(Sigma,V,p1,p2,root="sym",sep=FALSE){
  
}

diff.g=function(K1,K2,C,V,W,root="sym"){
  
  if(root=="sym"){
    
  }else{
    
  }
  
}



