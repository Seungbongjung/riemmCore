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

##################### Simulation

######### r=4
para.16.10.4=para.generate(4,16,10,100)
save(para.16.10.4,file="./numerical.simul/para.16.10.4.rds")
lambda=0.4

######### n=p/8 (20)

######### n=p/4 (40)

######### n=p/2 (80)

######### n=p (160)

######### n=2p (320)

######### r=8
para.16.10.8=para.generate(8,16,10,100)
save(para.16.10.8,file="./numerical.simul/para.16.10.8.rds")
lambda=0.4

######### n=p/8 (20)

######### n=p/4 (40)

######### n=p/2 (80)

######### n=p (160)

######### n=2p (320)

