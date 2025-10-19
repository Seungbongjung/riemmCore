################## Utensils

source("utensils.R")

optim.K1=function(dat,K1,K1.inv,K2.inv,A,lambda){
  
  p1=nrow(K1.inv); p2=nrow(K2.inv); p=p1*p2 
  n=dim(dat)[1]
  r=ncol(A)
  
  A.svd=svd(A,r)
  U=A.svd$u
  d=A.svd$d
  alpha=d^2/((1-lambda)*d^2+lambda)
  
  anc.array=apply(aperm(dat,c(2,3,1)),3,FUN=function(X){K1.inv%*%X%*%t(K2.inv)})
  anc.array=aperm(array(anc.array,dim=c(p1,p2,n)),c(3,1,2))
  anc.array.sq=array.sq(anc.array)
  
  euclid.deriv=matrix(0,p1,p1)
  for(i in 1:n){
    temp=anc.array.sq[i,,]
    
    euclid.deriv=-2/(n*lambda)*sym.mat(temp%*%K1.inv)
    for(j in 1:r){
      temp2=anc.array[i,,]%*%t(matrix(U[,j],p1,p2))
      euclid.deriv=euclid.deriv+2*(1-lambda)/(n*lambda)*alpha[j]*sum(diag(temp2))*sym.mat(temp2%*%K1.inv)
    }
  }
  
  grad=K1%*%euclid.deriv%*%K1
  grad=det.proj(K1,K1.inv,grad)
  grad=grad[lower.tri(grad,diag=TRUE)]
  
  coeff=NULL
  for(j in 1:p1){
    for(i in j:p1){
      basis=matrix(0,p1,p1)
      if(i==j){
        basis[i,j]=1
      }else{
        basis[i,j]=basis[j,i]=1
      }
      basis=det.proj(K1,K1.inv,basis)
      hess=matrix(0,p1,p1)
      
      for(k in 1:n){
        temp=anc.array.sq[k,,]
        hess=hess+2/(n*lambda)*(sym.mat(K1.inv%*%basis%*%temp%*%K1.inv)+sym.mat(temp%*%t(basis)%*%t(K1.inv)%*%K1.inv)+
                                  sym.mat(temp%*%K1.inv%*%basis%*%K1.inv))
        for(l in 1:r){
          temp2=anc.array[k,,]%*%t(matrix(U[,l],p1,p2))
          hess=hess-2*(1-lambda)/(n*lambda)*alpha[l]*sum(diag(K1.inv%*%basis%*%temp2))*sym.mat(temp2%*%K1.inv)-2*(1-lambda)/(lambda*n)*alpha[l]*sum(diag(temp2))*(sym.mat(K1.inv%*%basis%*%temp2%*%K1.inv)+sym.mat(temp2%*%K1.inv%*%basis%*%K1.inv))
        }
      }
      temp=basis%*%euclid.deriv%*%K1/2+K1%*%euclid.deriv%*%basis/2+K1%*%hess%*%K1
      temp=det.proj(K1,K1.inv,temp)
      temp=temp[lower.tri(temp,diag=TRUE)]
      coeff=cbind(coeff,temp)
    }
  }
  
  v=-lsmr(coeff,grad)$x
  V=matrix(0,p1,p1)
  V[lower.tri(V,diag=TRUE)]=v; V=V+t(V)-diag(diag(V))
  K1=spd.geodesic(K1,det.proj(K1,K1.inv,V))
  
  return(K1)
  
}

optim.K2=function(dat,K1.inv,K2,K2.inv,A,lambda){
  
  p1=nrow(K1.inv); p2=nrow(K2.inv); p=p1*p2 
  n=dim(dat)[1]
  r=ncol(A)
  
  A.svd=svd(A,r)
  U=A.svd$u
  d=A.svd$d
  alpha=d^2/((1-lambda)*d^2+lambda)
  
  anc.array=apply(aperm(dat,c(2,3,1)),3,FUN=function(X){K1.inv%*%X%*%t(K2.inv)})
  anc.array=array(anc.array,dim=c(p1,p2,n))
  anc.array=aperm(array(apply(anc.array,3,FUN=t),dim=c(p2,p1,n)),c(3,1,2))
  anc.array.sq=array.sq(anc.array)
  
  euclid.deriv=matrix(0,p2,p2)
  for(i in 1:n){
    temp=anc.array.sq[i,,]
    euclid.deriv=-2/(n*lambda)*sym.mat(temp%*%K2.inv)
    for(j in 1:r){
      temp2=anc.array[i,,]%*%matrix(U[,j],p1,p2)
      euclid.deriv=euclid.deriv+2*(1-lambda)/(n*lambda)*alpha[j]*sum(diag(temp2))*sym.mat(temp2%*%K2.inv)
    }
  }
  grad=K2%*%euclid.deriv%*%K2
  grad=det.proj(K2,K2.inv,grad)
  grad=grad[lower.tri(grad,diag=TRUE)]
  
  
  coeff=NULL
  for(j in 1:p2){
    for(i in j:p2){
      basis=matrix(0,p2,p2)
      if(i==j){
        basis[i,j]=1
      }else{
        basis[i,j]=basis[j,i]=1
      }
      basis=det.proj(K2,K2.inv,basis)
      hess=matrix(0,p2,p2)
      
      for(k in 1:n){
        temp=anc.array.sq[k,,]
        hess=hess+2/(n*lambda)*(sym.mat(K2.inv%*%basis%*%temp%*%K2.inv)+sym.mat(temp%*%t(basis)%*%t(K2.inv)%*%K2.inv)+
                                  sym.mat(temp%*%K2.inv%*%basis%*%K2.inv))
        for(l in 1:r){
          temp2=anc.array[k,,]%*%matrix(U[,l],p1,p2)
          hess=hess-2*(1-lambda)/(lambda*n)*alpha[l]*sum(diag(temp2))*(sym.mat(K2.inv%*%basis%*%temp2%*%K2.inv)+sym.mat(temp2%*%K2.inv%*%basis%*%K2.inv))-
            2*(1-lambda)/(n*lambda)*alpha[l]*sum(diag(K2.inv%*%basis%*%temp2))*sym.mat(temp2%*%K2.inv)
        }
      }
      temp=basis%*%euclid.deriv%*%K2/2+K2%*%euclid.deriv%*%basis/2+K2%*%hess%*%K2
      temp=det.proj(K2,K2.inv,temp)
      temp=temp[lower.tri(temp,diag=TRUE)]
      coeff=cbind(coeff,temp)
    }
  }
  
  v=-lsmr(coeff,grad)$x
  V=matrix(0,p2,p2)
  V[lower.tri(V,diag=TRUE)]=v; V=V+t(V)-diag(diag(V))
  K2=spd.geodesic(K2,det.proj(K2,K2.inv,V))
  
  return(K2)
  
}


optim.A=function(S,K1.inv,K2.inv,nu,A,lambda){
  
  K.inv=kronecker(K2.inv,K1.inv)/nu
  
  p1=ncol(K1.inv); p2=ncol(K2.inv); p=p1*p2
  r=ncol(A)
  
  A.svd=svd(A,r)
  U=A.svd$u
  d=A.svd$d
  alpha=d^2/((1-lambda)*d^2+lambda)
  
  S.whit=K.inv%*%S%*%t(K.inv)
  C.inv=diag(p)/lambda-(1-lambda)/lambda*U%*%diag(alpha)%*%t(U)
  
  euclid.deriv=-2*(1-lambda)*C.inv%*%S.whit%*%C.inv%*%A+2*(1-lambda)*C.inv%*%A
  J=J.mat(A,p1)
  J.ginv=pinv(J)
  proj.J=J.ginv%*%J
  proj.null=diag(p*r)-proj.J
  euclid.deriv=as.numeric(euclid.deriv)
  
  grad=euclid.deriv-proj.J%*%euclid.deriv
  
  coeff=NULL
  anc.mat1=C.inv%*%S.whit%*%C.inv
  
  euclid.hess.ftn=function(i){
    basis=matrix(proj.null[,i],p,r)
    anc.mat2=C.inv%*%(A%*%t(basis)+basis%*%t(A))
    euclid.hess=-2*(1-lambda)*anc.mat1%*%basis+2*(1-lambda)*C.inv%*%basis+
      2*(1-lambda)^2*C.inv%*%(A%*%t(basis)+basis%*%t(A))%*%anc.mat1%*%A+
      2*(1-lambda)^2*C.inv%*%S.whit%*%anc.mat2%*%C.inv%*%A-
      2*(1-lambda)^2*anc.mat2%*%C.inv%*%A
    J.basis=J.mat(basis,p1)
    temp1=as.numeric(euclid.hess)
    temp2=as.numeric((J.ginv%*%J.basis%*%proj.null+proj.null%*%t(J.basis)%*%t(J.ginv))%*%euclid.deriv)
    return(temp1-temp2)
  }
  
  coeff=mclapply(X=1:(p*r),FUN=euclid.hess.ftn,mc.cores=max(1,detectCores()-1))
  coeff=proj.null%*%do.call(cbind,coeff)
  
  v=-lsmr(coeff,grad)$x
  v=v-proj.J%*%v
  V=matrix(v,p,r)
  M=A%*%t(V)+V%*%t(A)
  M=A%*%t(A)+M
  M=eigs_sym(M,r)
  A=M$vectors%*%diag(sqrt(M$values+1e-03))
  A.kcd=covKCD(A%*%t(A),p1,p2)
  A=kronecker(sym.inv.root(A.kcd$K2),sym.inv.root(A.kcd$K1))%*%A
  
  return(A)
  
}

optim.nu=function(S,K1.inv,K2.inv,A,lambda){
  
  K.inv=kronecker(K2.inv,K1.inv)
  
  p1=ncol(K1.inv); p2=ncol(K2.inv); p=p1*p2
  r=ncol(A)
  
  A.svd=svd(A)
  U=A.svd$u
  d=A.svd$d
  alpha=d^2/((1-lambda)*d^2+lambda)
  
  S.whit=K.inv%*%S%*%t(K.inv); R=t(U)%*%S.whit%*%U
  
  nu=sqrt((sum(diag(S.whit))/lambda-(1-lambda)/lambda*sum(diag(R)*alpha))/p)
  return(nu)
  
}

optim.lambda=function(S,K1.inv,K2.inv,nu,A){
  
  K.inv=kronecker(K2.inv,K1.inv)/nu
  S.whit=K.inv%*%S%*%t(K.inv)
  p=nrow(A); r=ncol(A)
  A.svd=svd(A)
  U=A.svd$u
  d=A.svd$d
  R=t(U)%*%S.whit%*%U
  
  ell=function(lambda){
    alpha=d^2/((1-lambda)*d^2+lambda)
    temp=sum(diag(S.whit))/lambda-(1-lambda)/lambda*sum(diag(R)*alpha)+sum(log((1-lambda)*d^2+lambda))+(p-r)*log(lambda)
    return(temp)
  }
  lambda.optim=optimize(ell,interval=c(1e-3,1))
  return(lambda.optim$minimum)
  
}

obj=function(S,K1.inv,K2.inv,nu,A,lambda){
  
  p1=ncol(K1.inv)
  K.inv=kronecker(K2.inv,K1.inv)/nu
  S.whit=K.inv%*%S%*%t(K.inv)
  p=nrow(A); r=ncol(A)
  A.svd=svd(A,r)
  U=A.svd$u
  d=A.svd$d
  C=U%*%diag(d^2/((1-lambda)*d^2+lambda))%*%t(U)
  
  val=sum(diag(S.whit))/lambda-(1-lambda)/lambda*sum(diag(S.whit%*%C))+sum(log((1-lambda)*d^2+lambda))+(p-r)*log(lambda)+2*p*log(nu)
  return(val)
  
}

optim.lik=function(dat,r,eps=5e-4,max.iter=100,center=TRUE){
  
  n=dim(dat)[1]; p1=dim(dat)[2]; p2=dim(dat)[3]; p=p1*p2

  if(center==TRUE){
    dat=apply(dat,MARGIN=c(2,3),FUN=scale,scale=FALSE)   
  }
  
  S=array.cov(dat)
  S.KCD=covKCD(S,p1,p2)
  K1=S.KCD$K1; K2=S.KCD$K2
  
  K1=sym.root(K1); K2=sym.root(K2)
  c=det(K1)^(1/p1); d=det(K2)^(1/p2)
  K1=K1/c; K2=K2/d
  nu=c*d
  
  K1.inv=solve(K1); K2.inv=solve(K2)
  
  C=S.KCD$C; C.eig=eigs_sym(C,r)
  lambda=(p-sum(C.eig$values))/(p-r) 
  
  C.hat=C.eig$vectors%*%diag(C.eig$values)%*%t(C.eig$vectors)
  C.hat=covKCD(C.hat,p1,p2)$C 
  A=eigs_sym(C.hat,r)
  A=A$vectors%*%diag(sqrt(A$values))
  
  lik1=obj(S,K1.inv,K2.inv,nu,A,lambda)
  lik0=lik1/2
  
  lik=c(lik1)
  iter=0
  
  K1.temp=K1; K1.inv.temp=K1.inv; K2.temp=K2; K2.inv.temp=K2.inv; nu.temp=nu
  A.temp=A; lambda.temp=lambda
  
  start=Sys.time()
  while(abs(lik1-lik0)/abs(lik1)>eps && iter<max.iter){
    
    lik0=lik1
    lik
    K1.temp=optim.K1(dat,K1.temp,K1.inv.temp,K2.inv.temp,A.temp,lambda.temp); K1.inv.temp=solve(K1.temp) 
    K2.temp=optim.K2(dat,K1.inv.temp,K2.temp,K2.inv.temp,A.temp,lambda.temp); K2.inv.temp=solve(K2.temp)
    nu.temp=optim.nu(S,K1.inv.temp,K2.inv.temp,A.temp,lambda.temp)
    A.temp=optim.A(S,K1.inv.temp,K2.inv.temp,nu.temp,A.temp,lambda.temp)
    lambda.temp=optim.lambda(S,K1.inv.temp,K2.inv.temp,nu.temp,A.temp)
    lik1=obj(S,K1.inv.temp,K2.inv.temp,nu.temp,A.temp,lambda.temp)
    
    if(lik1>lik0){
      rel.conv=0
      break
    }else{
      rel.conv=1
      iter=iter+1
      lik=c(lik,lik1)
      K1=K1.temp; K2=K2.temp; A=A.temp; nu=nu.temp; lambda=lambda.temp
    }
  }
  end=Sys.time()
  
  if(iter>=max.iter){
    convergence=FALSE
  }else{
    convergence=TRUE
  }
  
  return(list(K1=K1,K2=K2,nu=nu,A=A,lambda=lambda,lik=lik,rel.conv=rel.conv,convergence=convergence,iter=iter,time=as.numeric(difftime(end,start,units="mins"))))
  
}
