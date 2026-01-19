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
  cse=matrix(0,nrow=iter,ncol=3)
  pi.core.init.ai=matrix(0,nrow=iter,ncol=3)
  pi.core.init.chol=matrix(0,nrow=iter,ncol=2)
  pi.core.ai=list()
  pi.core.chol=list()

  colnames(pi.core.init.ai)=c("K","C","Sigma")
  colnames(pi.core.init.chol)=c("C","Sigma")
  colnames(kro.mle)=c("C","Sigma")
  colnames(cse)=c("C","Sigma","w")
  
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
    cse[i,1:2]=c(svds((1-w)*C.init+w*diag(p)-C.true,1)$d,svds(cse.est-Sigma,1)$d)/c(C.true.val,Sigma.val)
    cse[i,3]=w
    
    pi.ai.est=optim.lik(dat,r,eps=eps,max.iter=max.iter,center=TRUE,metric="AI")
    pi.chol.est=optim.lik(dat,r,eps=eps,max.iter=max.iter,center=TRUE,metric="Chol")
    
    pi.ai.cov.est=pi.core.cov(pi.ai.est)
    pi.ai.cov.eval=c(svds(pi.ai.cov.est$K.est-K.true,1)$d,svds(pi.ai.cov.est$C.est-C.true,1)$d,svds(pi.ai.cov.est$Sigma.est-Sigma,1)$d)/c(K.true.val,C.true.val,Sigma.val)
    names(pi.ai.cov.eval)=c("K","C","Sigma")
    pi.core.ai[[i]]=list(eval=pi.ai.cov.eval,lik=pi.ai.est$lik,lambda=pi.ai.est$lambda,rel.conv=pi.ai.est$rel.conv,convergence=pi.ai.est$convergence,iter=pi.ai.est$iter,
                         time=pi.ai.est$time)

    pi.chol.cov.est=pi.core.cov(pi.chol.est)
    pi.chol.cov.eval=c(svds(pi.chol.cov.est$K.est-K.true,1)$d,svds(pi.chol.cov.est$C.est-O.true%*%C.true%*%t(O.true),1)$d,svds(pi.chol.cov.est$Sigma.est-Sigma,1)$d)/c(K.true.val,C.true.val,Sigma.val)
    names(pi.chol.cov.eval)=c("K","C","Sigma")
    pi.core.chol[[i]]=list(eval=pi.chol.cov.eval,lik=pi.chol.est$lik,lambda=pi.chol.est$lambda,rel.conv=pi.chol.est$rel.conv,convergence=pi.chol.est$convergence,iter=pi.chol.est$iter,
                         time=pi.chol.est$time)

    if(i%%print.iter==0){
      print(paste(i,"th iteration done",sep=""))
    }
  }
  return(list(pi.core.init.ai=pi.core.init.ai,pi.core.init.chol=pi.core.init.chol,kro.mle=kro.mle,cse=cse,pi.core.ai=pi.core.ai,pi.core.chol=pi.core.chol))
}

##################### Simulation

######### (p1,p2)=(16,12)
######### r=3
para.16.12.3=para.generate(3,16,12,100)
save(para.16.12.3,file="./numerical.simul_picore/para.16.12.3.rds")

######### lambda=0.2

lambda=0.2

######### n=p/8 (24)
pi.core.fit.16.12.3.24.lambda1=optim.iter(24,para.16.12.3,lambda,100,100,1)
save(pi.core.fit.16.12.3.24.lambda1,file="./numerical.simul_picore/pi.core.fit.16.12.3.24.lambda1.rds")

######### n=p/4 (48)
pi.core.fit.16.12.3.48.lambda1=optim.iter(48,para.16.12.3,lambda,100,100,1)
save(pi.core.fit.16.12.3.48.lambda1,file="./numerical.simul_picore/pi.core.fit.16.12.3.48.lambda1.rds")  

######### n=p/2 (96)
pi.core.fit.16.12.3.96.lambda1=optim.iter(96,para.16.12.3,lambda,100,100,1)
save(pi.core.fit.16.12.3.96.lambda1,file="./numerical.simul_picore/pi.core.fit.16.12.3.96.lambda1.rds")

######### n=p (192)
pi.core.fit.16.12.3.192.lambda1=optim.iter(192,para.16.12.3,lambda,100,100,1)
save(pi.core.fit.16.12.3.192.lambda1,file="./numerical.simul_picore/pi.core.fit.16.12.3.192.lambda1.rds")

######### n=2p (384)
pi.core.fit.16.12.3.384.lambda1=optim.iter(384,para.16.12.3,lambda,100,100,1)
save(pi.core.fit.16.12.3.384.lambda1,file="./numerical.simul_picore/pi.core.fit.16.12.3.384.lambda1.rds")

######### lambda=0.4

lambda=0.4

######### n=p/8 (24)
pi.core.fit.16.12.3.24.lambda2=optim.iter(24,para.16.12.3,lambda,100,100,1)
save(pi.core.fit.16.12.3.24.lambda2,file="./numerical.simul_picore/pi.core.fit.16.12.3.24.lambda2.rds")

######### n=p/4 (48)
pi.core.fit.16.12.3.48.lambda2=optim.iter(48,para.16.12.3,lambda,100,100,1)
save(pi.core.fit.16.12.3.48.lambda2,file="./numerical.simul_picore/pi.core.fit.16.12.3.48.lambda2.rds")  

######### n=p/2 (96)
pi.core.fit.16.12.3.96.lambda2=optim.iter(96,para.16.12.3,lambda,100,100,1)
save(pi.core.fit.16.12.3.96.lambda2,file="./numerical.simul_picore/pi.core.fit.16.12.3.96.lambda2.rds")

######### n=p (192)
pi.core.fit.16.12.3.192.lambda2=optim.iter(192,para.16.12.3,lambda,100,100,1)
save(pi.core.fit.16.12.3.192.lambda2,file="./numerical.simul_picore/pi.core.fit.16.12.3.192.lambda2.rds")

######### n=2p (384)
pi.core.fit.16.12.3.384.lambda2=optim.iter(384,para.16.12.3,lambda,100,100,1)
save(pi.core.fit.16.12.3.384.lambda2,file="./numerical.simul_picore/pi.core.fit.16.12.3.384.lambda2.rds")

######### lambda=0.6

lambda=0.6

######### n=p/8 (24)
pi.core.fit.16.12.3.24.lambda3=optim.iter(24,para.16.12.3,lambda,100,100,1)
save(pi.core.fit.16.12.3.24.lambda3,file="./numerical.simul_picore/pi.core.fit.16.12.3.24.lambda3.rds")

######### n=p/4 (48)
pi.core.fit.16.12.3.48.lambda3=optim.iter(48,para.16.12.3,lambda,100,100,1)
save(pi.core.fit.16.12.3.48.lambda3,file="./numerical.simul_picore/pi.core.fit.16.12.3.48.lambda3.rds")  

######### n=p/2 (96)
pi.core.fit.16.12.3.96.lambda3=optim.iter(96,para.16.12.3,lambda,100,100,1)
save(pi.core.fit.16.12.3.96.lambda3,file="./numerical.simul_picore/pi.core.fit.16.12.3.96.lambda3.rds")

######### n=p (192)
pi.core.fit.16.12.3.192.lambda3=optim.iter(192,para.16.12.3,lambda,100,100,1)
save(pi.core.fit.16.12.3.192.lambda3,file="./numerical.simul_picore/pi.core.fit.16.12.3.192.lambda3.rds")

######### n=2p (384)
pi.core.fit.16.12.3.384.lambda3=optim.iter(384,para.16.12.3,lambda,100,100,1)
save(pi.core.fit.16.12.3.384.lambda3,file="./numerical.simul_picore/pi.core.fit.16.12.3.384.lambda3.rds")

######### lambda=0.8

lambda=0.8

######### n=p/8 (24)
pi.core.fit.16.12.3.24.lambda4=optim.iter(24,para.16.12.3,lambda,100,100,1)
save(pi.core.fit.16.12.3.24.lambda4,file="./numerical.simul_picore/pi.core.fit.16.12.3.24.lambda4.rds")

######### n=p/4 (48)
pi.core.fit.16.12.3.48.lambda4=optim.iter(48,para.16.12.3,lambda,100,100,1)
save(pi.core.fit.16.12.3.48.lambda4,file="./numerical.simul_picore/pi.core.fit.16.12.3.48.lambda4.rds")  

######### n=p/2 (96)
pi.core.fit.16.12.3.96.lambda4=optim.iter(96,para.16.12.3,lambda,100,100,1)
save(pi.core.fit.16.12.3.96.lambda4,file="./numerical.simul_picore/pi.core.fit.16.12.3.96.lambda4.rds")

######### n=p (192)
pi.core.fit.16.12.3.192.lambda4=optim.iter(192,para.16.12.3,lambda,100,100,1)
save(pi.core.fit.16.12.3.192.lambda4,file="./numerical.simul_picore/pi.core.fit.16.12.3.192.lambda4.rds")

######### n=2p (384)
pi.core.fit.16.12.3.384.lambda4=optim.iter(384,para.16.12.3,lambda,100,100,1)
save(pi.core.fit.16.12.3.384.lambda4,file="./numerical.simul_picore/pi.core.fit.16.12.3.384.lambda4.rds")

######### r=5
para.16.12.5=para.generate(5,16,12,100)
save(para.16.12.5,file="./numerical.simul_picore/para.16.12.5.rds")

######### lambda=0.2

lambda=0.2

######### n=p/8 (24)
pi.core.fit.16.12.5.24.lambda1=optim.iter(24,para.16.12.5,lambda,100,100,1)
save(pi.core.fit.16.12.5.24.lambda1,file="./numerical.simul_picore/pi.core.fit.16.12.5.24.lambda1.rds")

######### n=p/4 (48)
pi.core.fit.16.12.5.48.lambda1=optim.iter(48,para.16.12.5,lambda,100,100,1)
save(pi.core.fit.16.12.5.48.lambda1,file="./numerical.simul_picore/pi.core.fit.16.12.5.48.lambda1.rds")  

######### n=p/2 (96)
pi.core.fit.16.12.5.96.lambda1=optim.iter(96,para.16.12.5,lambda,100,100,1)
save(pi.core.fit.16.12.5.96.lambda1,file="./numerical.simul_picore/pi.core.fit.16.12.5.96.lambda1.rds")

######### n=p (192)
pi.core.fit.16.12.5.192.lambda1=optim.iter(192,para.16.12.5,lambda,100,100,1)
save(pi.core.fit.16.12.5.192.lambda1,file="./numerical.simul_picore/pi.core.fit.16.12.5.192.lambda1.rds")

######### n=2p (384)
pi.core.fit.16.12.5.384.lambda1=optim.iter(384,para.16.12.5,lambda,100,100,1)
save(pi.core.fit.16.12.5.384.lambda1,file="./numerical.simul_picore/pi.core.fit.16.12.5.384.lambda1.rds")

######### lambda=0.4

lambda=0.4

######### n=p/8 (24)
pi.core.fit.16.12.5.24.lambda2=optim.iter(24,para.16.12.5,lambda,100,100,1)
save(pi.core.fit.16.12.5.24.lambda2,file="./numerical.simul_picore/pi.core.fit.16.12.5.24.lambda2.rds")

######### n=p/4 (48)
pi.core.fit.16.12.5.48.lambda2=optim.iter(48,para.16.12.5,lambda,100,100,1)
save(pi.core.fit.16.12.5.48.lambda2,file="./numerical.simul_picore/pi.core.fit.16.12.5.48.lambda2.rds")  

######### n=p/2 (96)
pi.core.fit.16.12.5.96.lambda2=optim.iter(96,para.16.12.5,lambda,100,100,1)
save(pi.core.fit.16.12.5.96.lambda2,file="./numerical.simul_picore/pi.core.fit.16.12.5.96.lambda2.rds")

######### n=p (192)
pi.core.fit.16.12.5.192.lambda2=optim.iter(192,para.16.12.5,lambda,100,100,1)
save(pi.core.fit.16.12.5.192.lambda2,file="./numerical.simul_picore/pi.core.fit.16.12.5.192.lambda2.rds")

######### n=2p (384)
pi.core.fit.16.12.5.384.lambda2=optim.iter(384,para.16.12.5,lambda,100,100,1)
save(pi.core.fit.16.12.5.384.lambda2,file="./numerical.simul_picore/pi.core.fit.16.12.5.384.lambda2.rds")

######### lambda=0.6

lambda=0.6

######### n=p/8 (24)
pi.core.fit.16.12.5.24.lambda3=optim.iter(24,para.16.12.5,lambda,100,100,1)
save(pi.core.fit.16.12.5.24.lambda3,file="./numerical.simul_picore/pi.core.fit.16.12.5.24.lambda3.rds")

######### n=p/4 (48)
pi.core.fit.16.12.5.48.lambda3=optim.iter(48,para.16.12.5,lambda,100,100,1)
save(pi.core.fit.16.12.5.48.lambda3,file="./numerical.simul_picore/pi.core.fit.16.12.5.48.lambda3.rds")  

######### n=p/2 (96)
pi.core.fit.16.12.5.96.lambda3=optim.iter(96,para.16.12.5,lambda,100,100,1)
save(pi.core.fit.16.12.5.96.lambda3,file="./numerical.simul_picore/pi.core.fit.16.12.5.96.lambda3.rds")

######### n=p (192)
pi.core.fit.16.12.5.192.lambda3=optim.iter(192,para.16.12.5,lambda,100,100,1)
save(pi.core.fit.16.12.5.192.lambda3,file="./numerical.simul_picore/pi.core.fit.16.12.5.192.lambda3.rds")

######### n=2p (384)
pi.core.fit.16.12.5.384.lambda3=optim.iter(384,para.16.12.5,lambda,100,100,1)
save(pi.core.fit.16.12.5.384.lambda3,file="./numerical.simul_picore/pi.core.fit.16.12.5.384.lambda3.rds")

######### lambda=0.8

lambda=0.8

######### n=p/8 (24)
pi.core.fit.16.12.5.24.lambda4=optim.iter(24,para.16.12.5,lambda,100,100,1)
save(pi.core.fit.16.12.5.24.lambda4,file="./numerical.simul_picore/pi.core.fit.16.12.5.24.lambda4.rds")

######### n=p/4 (48)
pi.core.fit.16.12.5.48.lambda4=optim.iter(48,para.16.12.5,lambda,100,100,1)
save(pi.core.fit.16.12.5.48.lambda4,file="./numerical.simul_picore/pi.core.fit.16.12.5.48.lambda4.rds")  

######### n=p/2 (96)
pi.core.fit.16.12.5.96.lambda4=optim.iter(96,para.16.12.5,lambda,100,100,1)
save(pi.core.fit.16.12.5.96.lambda4,file="./numerical.simul_picore/pi.core.fit.16.12.5.96.lambda4.rds")

######### n=p (192)
pi.core.fit.16.12.5.192.lambda4=optim.iter(192,para.16.12.5,lambda,100,100,1)
save(pi.core.fit.16.12.5.192.lambda4,file="./numerical.simul_picore/pi.core.fit.16.12.5.192.lambda4.rds")

######### n=2p (384)
pi.core.fit.16.12.5.384.lambda4=optim.iter(384,para.16.12.5,lambda,100,100,1)
save(pi.core.fit.16.12.5.384.lambda4,file="./numerical.simul_picore/pi.core.fit.16.12.5.384.lambda4.rds")


######### (p1,p2)=(18,8)
######### r=3
para.18.8.3=para.generate(3,18,8,100)
save(para.18.8.3,file="./numerical.simul_picore/para.18.8.3.rds")

######### lambda=0.2

lambda=0.2

######### n=p/8 (18)
pi.core.fit.18.8.3.18.lambda1=optim.iter(18,para.18.8.3,lambda,100,100,1)
save(pi.core.fit.18.8.3.18.lambda1,file="./numerical.simul_picore/pi.core.fit.18.8.3.18.lambda1.rds")

######### n=p/4 (36)
pi.core.fit.18.8.3.36.lambda1=optim.iter(36,para.18.8.3,lambda,100,100,1)
save(pi.core.fit.18.8.3.36.lambda1,file="./numerical.simul_picore/pi.core.fit.18.8.3.36.lambda1.rds")

######### n=p/2 (72)
pi.core.fit.18.8.3.72.lambda1=optim.iter(72,para.18.8.3,lambda,100,100,1)
save(pi.core.fit.18.8.3.72.lambda1,file="./numerical.simul_picore/pi.core.fit.18.8.3.72.lambda1.rds")

######### n=p (144)
pi.core.fit.18.8.3.144.lambda1=optim.iter(144,para.18.8.3,lambda,100,100,1)
save(pi.core.fit.18.8.3.144.lambda1,file="./numerical.simul_picore/pi.core.fit.18.8.3.144.lambda1.rds")

######### n=2p (288)
pi.core.fit.18.8.3.288.lambda1=optim.iter(288,para.18.8.3,lambda,100,100,1)
save(pi.core.fit.18.8.3.288.lambda1,file="./numerical.simul_picore/pi.core.fit.18.8.3.288.lambda1.rds")

######### lambda=0.4

lambda=0.4

######### n=p/8 (18)
pi.core.fit.18.8.3.18.lambda2=optim.iter(18,para.18.8.3,lambda,100,100,1)
save(pi.core.fit.18.8.3.18.lambda2,file="./numerical.simul_picore/pi.core.fit.18.8.3.18.lambda2.rds")

######### n=p/4 (36)
pi.core.fit.18.8.3.36.lambda2=optim.iter(36,para.18.8.3,lambda,100,100,1)
save(pi.core.fit.18.8.3.36.lambda2,file="./numerical.simul_picore/pi.core.fit.18.8.3.36.lambda2.rds")

######### n=p/2 (72)
pi.core.fit.18.8.3.72.lambda2=optim.iter(72,para.18.8.3,lambda,100,100,1)
save(pi.core.fit.18.8.3.72.lambda2,file="./numerical.simul_picore/pi.core.fit.18.8.3.72.lambda2.rds")

######### n=p (144)
pi.core.fit.18.8.3.144.lambda2=optim.iter(144,para.18.8.3,lambda,100,100,1)
save(pi.core.fit.18.8.3.144.lambda2,file="./numerical.simul_picore/pi.core.fit.18.8.3.144.lambda2.rds")

######### n=2p (288)
pi.core.fit.18.8.3.288.lambda2=optim.iter(288,para.18.8.3,lambda,100,100,1)
save(pi.core.fit.18.8.3.288.lambda2,file="./numerical.simul_picore/pi.core.fit.18.8.3.288.lambda2.rds")

######### lambda=0.6

lambda=0.6

######### n=p/8 (18)
pi.core.fit.18.8.3.18.lambda3=optim.iter(18,para.18.8.3,lambda,100,100,1)
save(pi.core.fit.18.8.3.18.lambda3,file="./numerical.simul_picore/pi.core.fit.18.8.3.18.lambda3.rds")

######### n=p/4 (36)
pi.core.fit.18.8.3.36.lambda3=optim.iter(36,para.18.8.3,lambda,100,100,1)
save(pi.core.fit.18.8.3.36.lambda3,file="./numerical.simul_picore/pi.core.fit.18.8.3.36.lambda3.rds")

######### n=p/2 (72)
pi.core.fit.18.8.3.72.lambda3=optim.iter(72,para.18.8.3,lambda,100,100,1)
save(pi.core.fit.18.8.3.72.lambda3,file="./numerical.simul_picore/pi.core.fit.18.8.3.72.lambda3.rds")

######### n=p (144)
pi.core.fit.18.8.3.144.lambda3=optim.iter(144,para.18.8.3,lambda,100,100,1)
save(pi.core.fit.18.8.3.144.lambda3,file="./numerical.simul_picore/pi.core.fit.18.8.3.144.lambda3.rds")

######### n=2p (288)
pi.core.fit.18.8.3.288.lambda3=optim.iter(288,para.18.8.3,lambda,100,100,1)
save(pi.core.fit.18.8.3.288.lambda3,file="./numerical.simul_picore/pi.core.fit.18.8.3.288.lambda3.rds")

######### lambda=0.8

lambda=0.8

######### n=p/8 (18)
pi.core.fit.18.8.3.18.lambda4=optim.iter(18,para.18.8.3,lambda,100,100,1)
save(pi.core.fit.18.8.3.18.lambda4,file="./numerical.simul_picore/pi.core.fit.18.8.3.18.lambda4.rds")

######### n=p/4 (36)
pi.core.fit.18.8.3.36.lambda4=optim.iter(36,para.18.8.3,lambda,100,100,1)
save(pi.core.fit.18.8.3.36.lambda4,file="./numerical.simul_picore/pi.core.fit.18.8.3.36.lambda4.rds")

######### n=p/2 (72)
pi.core.fit.18.8.3.72.lambda4=optim.iter(72,para.18.8.3,lambda,100,100,1)
save(pi.core.fit.18.8.3.72.lambda4,file="./numerical.simul_picore/pi.core.fit.18.8.3.72.lambda4.rds")

######### n=p (144)
pi.core.fit.18.8.3.144.lambda4=optim.iter(144,para.18.8.3,lambda,100,100,1)
save(pi.core.fit.18.8.3.144.lambda4,file="./numerical.simul_picore/pi.core.fit.18.8.3.144.lambda4.rds")

######### n=2p (288)
pi.core.fit.18.8.3.288.lambda4=optim.iter(288,para.18.8.3,lambda,100,100,1)
save(pi.core.fit.18.8.3.288.lambda4,file="./numerical.simul_picore/pi.core.fit.18.8.3.288.lambda4.rds")

######### r=5
para.18.8.5=para.generate(5,18,8,100)
save(para.18.8.5,file="./numerical.simul_picore/para.18.8.5.rds")

######### lambda=0.2

lambda=0.2

######### n=p/8 (18)
pi.core.fit.18.8.5.18.lambda1=optim.iter(18,para.18.8.5,lambda,100,100,1)
save(pi.core.fit.18.8.5.18.lambda1,file="./numerical.simul_picore/pi.core.fit.18.8.5.18.lambda1.rds")

######### n=p/4 (36)
pi.core.fit.18.8.5.36.lambda1=optim.iter(36,para.18.8.5,lambda,100,100,1)
save(pi.core.fit.18.8.5.36.lambda1,file="./numerical.simul_picore/pi.core.fit.18.8.5.36.lambda1.rds")

######### n=p/2 (72)
pi.core.fit.18.8.5.72.lambda1=optim.iter(72,para.18.8.5,lambda,100,100,1)
save(pi.core.fit.18.8.5.72.lambda1,file="./numerical.simul_picore/pi.core.fit.18.8.5.72.lambda1.rds")

######### n=p (144)
pi.core.fit.18.8.5.144.lambda1=optim.iter(144,para.18.8.5,lambda,100,100,1)
save(pi.core.fit.18.8.5.144.lambda1,file="./numerical.simul_picore/pi.core.fit.18.8.5.144.lambda1.rds")

######### n=2p (288)
pi.core.fit.18.8.5.288.lambda1=optim.iter(288,para.18.8.5,lambda,100,100,1)
save(pi.core.fit.18.8.5.288.lambda1,file="./numerical.simul_picore/pi.core.fit.18.8.5.288.lambda1.rds")

######### lambda=0.4

lambda=0.4

######### n=p/8 (18)
pi.core.fit.18.8.5.18.lambda2=optim.iter(18,para.18.8.5,lambda,100,100,1)
save(pi.core.fit.18.8.5.18.lambda2,file="./numerical.simul_picore/pi.core.fit.18.8.5.18.lambda2.rds")

######### n=p/4 (36)
pi.core.fit.18.8.5.36.lambda2=optim.iter(36,para.18.8.5,lambda,100,100,1)
save(pi.core.fit.18.8.5.36.lambda2,file="./numerical.simul_picore/pi.core.fit.18.8.5.36.lambda2.rds")

######### n=p/2 (72)
pi.core.fit.18.8.5.72.lambda2=optim.iter(72,para.18.8.5,lambda,100,100,1)
save(pi.core.fit.18.8.5.72.lambda2,file="./numerical.simul_picore/pi.core.fit.18.8.5.72.lambda2.rds")

######### n=p (144)
pi.core.fit.18.8.5.144.lambda2=optim.iter(144,para.18.8.5,lambda,100,100,1)
save(pi.core.fit.18.8.5.144.lambda2,file="./numerical.simul_picore/pi.core.fit.18.8.5.144.lambda2.rds")

######### n=2p (288)
pi.core.fit.18.8.5.288.lambda2=optim.iter(288,para.18.8.5,lambda,100,100,1)
save(pi.core.fit.18.8.5.288.lambda2,file="./numerical.simul_picore/pi.core.fit.18.8.5.288.lambda2.rds")

######### lambda=0.6

lambda=0.6

######### n=p/8 (18)
pi.core.fit.18.8.5.18.lambda3=optim.iter(18,para.18.8.5,lambda,100,100,1)
save(pi.core.fit.18.8.5.18.lambda3,file="./numerical.simul_picore/pi.core.fit.18.8.5.18.lambda3.rds")

######### n=p/4 (36)
pi.core.fit.18.8.5.36.lambda3=optim.iter(36,para.18.8.5,lambda,100,100,1)
save(pi.core.fit.18.8.5.36.lambda3,file="./numerical.simul_picore/pi.core.fit.18.8.5.36.lambda3.rds")

######### n=p/2 (72)
pi.core.fit.18.8.5.72.lambda3=optim.iter(72,para.18.8.5,lambda,100,100,1)
save(pi.core.fit.18.8.5.72.lambda3,file="./numerical.simul_picore/pi.core.fit.18.8.5.72.lambda3.rds")

######### n=p (144)
pi.core.fit.18.8.5.144.lambda3=optim.iter(144,para.18.8.5,lambda,100,100,1)
save(pi.core.fit.18.8.5.144.lambda3,file="./numerical.simul_picore/pi.core.fit.18.8.5.144.lambda3.rds")

######### n=2p (288)
pi.core.fit.18.8.5.288.lambda3=optim.iter(288,para.18.8.5,lambda,100,100,1)
save(pi.core.fit.18.8.5.288.lambda3,file="./numerical.simul_picore/pi.core.fit.18.8.5.288.lambda3.rds")

######### lambda=0.8

lambda=0.8

######### n=p/8 (18)
pi.core.fit.18.8.5.18.lambda4=optim.iter(18,para.18.8.5,lambda,100,100,1)
save(pi.core.fit.18.8.5.18.lambda4,file="./numerical.simul_picore/pi.core.fit.18.8.5.18.lambda4.rds")

######### n=p/4 (36)
pi.core.fit.18.8.5.36.lambda4=optim.iter(36,para.18.8.5,lambda,100,100,1)
save(pi.core.fit.18.8.5.36.lambda4,file="./numerical.simul_picore/pi.core.fit.18.8.5.36.lambda4.rds")

######### n=p/2 (72)
pi.core.fit.18.8.5.72.lambda4=optim.iter(72,para.18.8.5,lambda,100,100,1)
save(pi.core.fit.18.8.5.72.lambda4,file="./numerical.simul_picore/pi.core.fit.18.8.5.72.lambda4.rds")

######### n=p (144)
pi.core.fit.18.8.5.144.lambda4=optim.iter(144,para.18.8.5,lambda,100,100,1)
save(pi.core.fit.18.8.5.144.lambda4,file="./numerical.simul_picore/pi.core.fit.18.8.5.144.lambda4.rds")

######### n=2p (288)
pi.core.fit.18.8.5.288.lambda4=optim.iter(288,para.18.8.5,lambda,100,100,1)
save(pi.core.fit.18.8.5.288.lambda4,file="./numerical.simul_picore/pi.core.fit.18.8.5.288.lambda4.rds")

