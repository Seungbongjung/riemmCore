##################### Utensils

source("./riemmCore.R")

para.generate=function(r,p1,p2,seed=100){
  
  ###### Generate parameters for partial isotropy core covariance
  ###### Argument 
  ## r: partial isotropy rank
  ## p1: row dimension
  ## p2: column dimension
  ## seed: seed number
  
  set.seed(seed)
  
  X1=matrix(rnorm(p1^2*2),ncol=2*p1); K1=X1%*%t(X1)/(2*p1)
  X2=matrix(rnorm(p2^2*2),ncol=2*p2); K2=X2%*%t(X2)/(2*p2)
  
  A=matrix(rnorm(p1*p2*r),ncol=r); A.kcd=covKCD(A%*%t(A),p1,p2)
  A=kronecker(sym.inv.root(A.kcd$K2),sym.inv.root(A.kcd$K1))%*%A
  
  return(list(K1=K1,K2=K2,A=A))
}

dat.generate=function(n,para,lambda){

  ####### Generate data for specified covariance model
  ####### Argument
  ## n: number of samples 
  ## para: parameter for the specified model
  
  ####### Output
  ## n x p1 x p2 data array

  K1=para$K1; K2=para$K2; A=para$A
  K.root=kronecker(sym.root(K2),sym.root(K1))
  p=nrow(A)
  Sigma=K.root%*%((1-lambda)*A%*%t(A)+lambda*diag(p))%*%t(K.root)
  p1=ncol(K1); p2=ncol(K2)
    
  dat=matrix(rnorm(n*p),ncol=p)%*%sym.root(Sigma)
  dat.array=array(0,dim=c(n,p1,p2))
  for(i in 1:n){
    dat.array[i,,]=matrix(dat[i,],p1,p2)
  }
  return(dat.array)
  
}

pi.core.cov=function(est){
  
  K1=est$K1; K2=est$K2; nu=est$nu; A=est$A; lambda=est$lambda
  K.est=kronecker(K2,K1)*nu
  C.est=(1-lambda)*A%*%t(A)+lambda*diag(nrow(A))
  return(list(K.est=K.est%*%t(K.est),C.est=C.est,Sigma.est=K.est%*%C.est%*%t(K.est)))
}


optim.iter=function(n,para,lambda,iter,seed,print.iter,eps=0.001,max.iter=100){
  
  K1=para$K1; K2=para$K2; A=para$A
  K.root=kronecker(sym.root(K2),sym.root(K1))
  O.true=kronecker(solve(t(chol(K2)))%*%sym.root(K2),solve(t(chol(K1)))%*%sym.root(K1))
  p=nrow(A); r=ncol(A)
  
  p1=ncol(K1); p2=ncol(K2)

  K.true=K.root%*%t(K.root)
  C.true=(1-lambda)*A%*%t(A)+lambda*diag(p)
  Sigma=K.root%*%C.true%*%t(K.root)
  
  K.true.val=eigs_sym(K.true,1)$values
  C.true.val=eigs_sym(C.true,1)$values
  Sigma.val=eigs_sym(Sigma,1)$values
  
  kro.mle=matrix(0,nrow=iter,ncol=2)
  pt.est=matrix(0,nrow=iter,ncol=2)
  cse=matrix(0,nrow=iter,ncol=2)
  pi.core.init.ai=matrix(0,nrow=iter,ncol=3)
  pi.core.init.chol=matrix(0,nrow=iter,ncol=2)
  pi.core.ai=list()
  pi.core.chol=list()

  colnames(pi.core.init.ai)=c("K","C","Sigma")
  colnames(pi.core.init.chol)=c("C","Sigma")
  colnames(kro.mle)=c("C","Sigma")
  colnames(cse)=c("C","Sigma")
  
  set.seed(seed)
  
  for(i in 1:iter){
    dat=dat.generate(n,para,lambda)
    S=mcov(dat)
    S.kcd=covKCD(S,p1,p2)
   
    K.init=S.kcd$K; C.init=S.kcd$C
    K.init.root=kronecker(sym.root(S.kcd$K2),sym.root(S.kcd$K1))
    K.init.chol=kronecker(t(chol(S.kcd$K2)),t(chol(S.kcd$K1)))
    O=kronecker(solve(t(chol(S.kcd$K2)))%*%sym.root(S.kcd$K2),solve(t(chol(S.kcd$K1)))%*%sym.root(S.kcd$K1))
    C.rank_r=eigs_sym(C.init,r)
    lambda.init=(p-sum(C.rank_r$values))/(p-r)
    C.rank_r=(C.rank_r$vectors)%*%diag(C.rank_r$values)%*%t(C.rank_r$vectors)
    C.rank_r=covKCD(C.rank_r,p1,p2)$C
    C.rank_r=(1-lambda.init)*C.rank_r+lambda.init*diag(p)
    pi.core.init.ai[i,]=c(svds(K.init-K.true,1)$d,svds(C.rank_r-C.true,1)$d,svds(K.init.root%*%C.rank_r%*%K.init.root-Sigma,1)$d)/c(K.true.val,C.true.val,Sigma.val)
 
    C.rank_r=eigs_sym(O%*%C.init%*%t(O),r)
    lambda.init=(p-sum(C.rank_r$values))/(p-r)
    C.rank_r=(C.rank_r$vectors)%*%diag(C.rank_r$values)%*%t(C.rank_r$vectors)
    C.rank_r=covKCD(C.rank_r,p1,p2)$C
    C.rank_r=(1-lambda.init)*C.rank_r+lambda.init*diag(p)
    pi.core.init.chol[i,]=c(svds(C.rank_r-O.true%*%C.true%*%t(O.true),1)$d,svds(K.init.chol%*%C.rank_r%*%t(K.init.chol)-Sigma,1)$d)/c(C.true.val,Sigma.val)

    kro.mle[i,]=c(svds(diag(p)-C.true,1)$d,svds(K.init-Sigma,1)$d)/c(C.true.val,Sigma.val)
    w=attributes(covCSE(dat,n=n,p1=p1,p2=p2))$w
    cse.est=K.init.root%*%((1-w)*C.init+w*diag(p))%*%t(K.init.root)
    cse[i,]=c(svds((1-w)*C.init+w*diag(p)-C.true,1)$d,svds(cse.est-Sigma,1)$d)/c(C.true.val,Sigma.val)

    pi.ai.est=optim.lik(dat,r,eps=eps,max.iter=max.iter,center=TRUE,metric="AI")
    pi.chol.est=optim.lik(dat,r,eps=eps,max.iter=max.iter,center=TRUE,metric="Chol")
    
    pi.ai.cov.est=pi.core.cov(pi.ai.est)
    pi.ai.cov.eval=c(svds(pi.ai.cov.est$K.est-K.true,1)$d,svds(pi.ai.cov.est$C.est-C.true,1)$d,svds(pi.ai.cov.est$Sigma.est-Sigma,1)$d)/c(K.true.val,C.true.val,Sigma.val)
    names(pi.ai.cov.eval)=c("K","C","Sigma")
    pi.core.ai[[i]]=list(eval=pi.ai.cov.eval,lik=pi.ai.est$lik,rel.conv=pi.ai.est$rel.conv,convergence=pi.ai.est$convergence,iter=pi.ai.est$iter,
                         time=pi.ai.est$time)

    pi.chol.cov.est=pi.core.cov(pi.chol.est)
    pi.chol.cov.eval=c(svds(pi.chol.cov.est$K.est-K.true,1)$d,svds(pi.chol.cov.est$C.est-O.true%*%C.true%*%t(O.true),1)$d,svds(pi.chol.cov.est$Sigma.est-Sigma,1)$d)/c(K.true.val,C.true.val,Sigma.val)
    names(pi.chol.cov.eval)=c("K","C","Sigma")
    pi.core.chol[[i]]=list(eval=pi.chol.cov.eval,lik=pi.chol.est$lik,rel.conv=pi.chol.est$rel.conv,convergence=pi.chol.est$convergence,iter=pi.chol.est$iter,
                         time=pi.chol.est$time)

    if(i%%print.iter==0){
      print(paste(i,"th iteration done",sep=""))
    }
  }
  return(list(pi.core.init.ai=pi.core.init.ai,pi.core.init.chol=pi.core.init.chol,kro.mle=kro.mle,cse=cse,pi.core.ai=pi.core.ai,pi.core.chol=pi.core.chol))
}

##################### Simulation

######### r=4
para.16.10.4=para.generate(4,16,10,100)
save(para.16.10.4,file="./numerical.simul/para.16.10.4.rds")
lambda=0.4

######### n=p/8 (20)
pi.core.fit.16.10.4.20=optim.iter(20,para.16.10.4,lambda,100,100,1)
save(pi.core.fit.16.10.4.20,file="./numerical.simul/pi.core.fit.16.10.4.20.rds")

######### n=p/4 (40)
pi.core.fit.16.10.4.40=optim.iter(40,para.16.10.4,lambda,100,100,1)
save(pi.core.fit.16.10.4.40,file="./numerical.simul/pi.core.fit.16.10.4.40.rds")  

######### n=p/2 (80)
pi.core.fit.16.10.4.80=optim.iter(80,para.16.10.4,lambda,100,100,1)
save(pi.core.fit.16.10.4.80,file="./numerical.simul/pi.core.fit.16.10.4.80.rds")

######### n=p (160)
pi.core.fit.16.10.4.160=optim.iter(160,para.16.10.4,lambda,100,100,1)
save(pi.core.fit.16.10.4.160,file="./numerical.simul/pi.core.fit.16.10.4.160.rds")

######### n=2p (320)
pi.core.fit.16.10.4.320=optim.iter(320,para.16.10.4,lambda,100,100,1)
save(pi.core.fit.16.10.4.320,file="./numerical.simul/pi.core.fit.16.10.4.320.rds")

######### r=8
para.16.10.8=para.generate(8,16,10,100)
save(para.16.10.8,file="./numerical.simul/para.16.10.8.rds")
lambda=0.4

######### n=p/8 (20)
pi.core.fit.16.10.8.20=optim.iter(20,para.16.10.8,lambda,100,100,1)
save(pi.core.fit.16.10.8.20,file="./numerical.simul/pi.core.fit.16.10.8.20.rds")

######### n=p/4 (40)
pi.core.fit.16.10.8.40=optim.iter(40,para.16.10.8,lambda,100,100,1)
save(pi.core.fit.16.10.8.40,file="./numerical.simul/pi.core.fit.16.10.8.40.rds")

######### n=p/2 (80)
pi.core.fit.16.10.8.80=optim.iter(80,para.16.10.8,lambda,100,100,1)
save(pi.core.fit.16.10.8.80,file="./numerical.simul/pi.core.fit.16.10.8.80.rds")

######### n=p (160)
pi.core.fit.16.10.8.160=optim.iter(160,para.16.10.8,lambda,100,100,1)
save(pi.core.fit.16.10.8.160,file="./numerical.simul/pi.core.fit.16.10.8.160.rds")

######### n=2p (320)
pi.core.fit.16.10.8.320=optim.iter(320,para.16.10.8,lambda,100,100,1)
save(pi.core.fit.16.10.8.320,file="./numerical.simul/pi.core.fit.16.10.8.320.rds")
