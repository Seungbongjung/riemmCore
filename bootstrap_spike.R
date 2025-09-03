################## Utensils

library(RSpectra)
library(covKCD)
library(RMT4DS)

#### Estimation of the number of spikes using Passemier and Yao's method

thres.C=function(p,n,iter,seed,alpha){
  
  set.seed(seed)
  thres=NULL
  
  for(i in 1:iter){
    
    X=matrix(rnorm(n*p),ncol=p)
    S=t(X)%*%X/n
    S.eig=eigs_sym(S,2)$values
    
    thres=c(thres,S.eig[1]-S.eig[2])
    
  }
  return(quantile(thres,alpha))
}

dat.center=function(dat){
  return(apply(dat,MARGIN=c(2,3),FUN=scale,scale=FALSE))
}

C=thres.C(320,400,1000,1,0.05)


num.spike1=function(dat,C){
  
  # dat : centered data array
  
  dat.dim=dim(dat)
  n=dat.dim[1]; p1=dat.dim[2]; p2=dat.dim[3]
  p=p1*p2
  
  S=matrix(0,p,p)
  for(i in 1:n){
    S=as.numeric(dat[,,i])%*%t(as.numeric(dat[,,i]))/n
    S=covKCD(S,p1,p2)$C
  }
  
  if(p<=n){
    S.eig=eigen(S)$values
  }else{
    S.eig=eigs_sym(S,n)$values
  }
  l=length(S.eig)
  
  k=1; terminate=0
  spike1=S.eig[1]-S.eig[2]; spike2=S.eig[2]-S.eig[3]
  while(spike1>=C || spike2>=C){
    if(k==(l-1)){
      terminate=1
      break
    }
    k=k+1
    spike1=spike2; spike2=S.eig[k+1]-S.eig[k+2]
    
  }
  temp=c(k,terminate); names(temp)=c("Num.Spike","Termination")
  return(temp)
}

#### Estimation of the number of spikes using Jiang's method

sym.inv.root=function(Sigma){
  
  Sigma.eig=eigen(Sigma)
  Q=Sigma.eig$vectors
  D=Sigma.eig$values
  return(Q%*%diag(1/sqrt(D))%*%t(Q))
  
}

num.spike2=function(dat,k.lower,k.upper){
  
  # dat : centered data array
  
  dat.dim=dim(dat)
  n=dat.dim[1]; p1=dat.dim[2]; p2=dat.dim[3]
  p=p1*p2
  
  dat.mat=matrix(dat,nrow=n)
  c=p/n
  
  S=t(dat.mat)%*%dat.mat/n
  S.KCD=covKCD(S,p1,p2)
  S.core=S.KCD$C 
  S.K1=S.KCD$K1; S.K2=S.KCD$K2
  S.K1.root=sym.inv.root(S.K1)
  S.K2.root=sym.inv.root(S.K2)
  
  for(i in 1:n){
    dat[i,,]=S.K1.root%*%dat[i,,]%*%t(S.K2.root)
  }
  dat.mat=matrix(dat,nrow=n)
  
  test.stat=NULL
  p.val=NULL
  est.spike=list()
  j=0
  S.core.eig=eigen(S.core)$values
  
  for(i in k.lower:k.upper){
    
    j=j+1
    
    alpha=MPEst(dat.mat,n_spike=i)$d[1:i]
    m=as.numeric(table(alpha))
    l1=length(m)
    m=m[l1:1]
    alpha=unique(alpha)
    
    Sigma.hat=(p-sum(S.core.eig[1:i]))/(p-i)
    
    Sigma.c=0
    for(k1 in 1:l1){
      Sigma.c=Sigma.c+alpha[k1]*c*m[k1]/(alpha[k1]/Sigma.hat-1)
    }
    Sigma.c=Sigma.hat+Sigma.c/(p-i)
    
    nu=2*c
    temp=1/sqrt(nu)*(p/Sigma.c-sum(S.core.eig[1:i])/Sigma.c-(p-i)+sum(m*alpha/Sigma.c/(alpha/Sigma.c-1))*c)
    test.stat=c(test.stat,temp)
    p.val=c(p.val,1-pnorm(abs(temp)))
    est.spike[[j]]=list(spike=alpha,non.spike=Sigma.c,m=m)
  }
  return(list(test.stat=test.stat,p.val=p.val,est.spike=est.spike))
}

#### Implementation of Li, Han, and Yao's Test 

spike.test=function(dat,n.spike,level){
 
  dat.dim=dim(dat)
  n=dat.dim[1]; p1=dat.dim[2]; p2=dat.dim[3]
  p=p1*p2
  
  dat.mat=matrix(dat,nrow=n)
  y=p/n
  c=1+sqrt(y)
  
  S=t(dat.mat)%*%dat.mat/n
  S.core=covKCD(S,p1,p2)$C 
  
  S.core.eig=eigen(S.core)$values
  S.core.spike=S.core.eig[1:n.spike]
  test.stat=S.core.spike[n.spike]/((p-sum(S.core.spike))/(p-n.spike))
  
  alpha=NULL
  for(i in 1:n.spike){
   temp=-(1-y)/S.core.spike[i]+1/n*sum(1/(S.core.eig[-(1:n.spike)]-S.core.spike[i])) 
   alpha=c(alpha,-1/temp)
  }
  t=alpha/(p-sum(alpha))*(p-n.spike)

  z=qnorm(1-level)
  q=function(x){
    anc1=1/(1-sum(c(y/(1-1/t[1:(n.spike-1)]),y/(1-1/x)))/(p-n.spike))
    anc2=x+y/(1-1/x)
    anc3=1-y/(x-1)^2
    sigma=2*x^2*anc3*anc1^2-4*y*x^2/((p-i)*(x-1)^2)*anc2*anc1^3+2*y*n/(p-n.spike)^2*anc2^2*anc1^4+
      x^2*anc3^2/n*anc1^2*(4*y*x/(3*(x-1)^3)-4*y*x/(3*(x-1)^3*anc3^3)+2*y^2*x^2/(3*(x-1)^6*anc3^4)+
                             2*y*x^2/(x-1)^4+4*y^2*x^2/(3*(x-1)^6*anc3))
    
    obj=(x+y/(1-1/x))/(1-sum(c(y/(1-1/t[1:(n.spike-1)]),y/(1-1/x)))/(p-n.spike))+sigma*z/sqrt(n)
    return(obj)
  }
  q.optim=optimize(q,interval=c(c,t[(n.spike-1)]-0.001))$objective
  if(test.stat<q.optim){
    reject=1
    return(reject)
  }else{
    reject=0
    return(list(reject=reject,alpha=alpha))
  }
  
}

spike.pval=function(dat,n.spike,boot.size,seed){
 
  set.seed(seed)

  dat.dim=dim(dat)
  n=dat.dim[1]; p1=dat.dim[2]; p2=dat.dim[3]
  p=p1*p2
  
  dat.mat=matrix(dat,nrow=n)
  S=t(dat.mat)%*%dat.mat/n
  S.core=covKCD(S,p1,p2)$C 
  
  S.core.spike=eigs_sym(S.core,n.spike)$values
  test.stat=S.core.spike[n.spike]/((p-sum(S.core.spike))/(p-n.spike))
  
  test.stat.boot=NULL
  
  for(i in 1:boot.size){
    
    samp=sample(1:n,n,replace=TRUE)
    temp=dat.mat[samp,]
    S.temp=t(temp)%*%temp/n
    S.core.temp=covKCD(S.temp,p1,p2)$C 
    
    temp.spike=eigs_sym(S.core.temp,n.spike)$values
    test.stat.boot=c(test.stat.boot,temp.spike[n.spike]/((p-sum(temp.spike))/(p-n.spike)))
    
    if(i%%10==0){
      print(paste(i,"th iteration",sep=""))
    }
  }
  
  return((sum(test.stat.boot<=test.stat)+1)/(boot.size+1))
  
}

