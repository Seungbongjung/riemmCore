######## Utensils

library(job)
library(RSpectra)
library(covKCD)

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
  return(quantile(thres,1-alpha))
}

dat.center=function(dat){
  return(apply(dat,MARGIN=c(2,3),FUN=scale,scale=FALSE))
}

num.spike1=function(dat,C){
  
  # dat : centered data array
  
  dat.dim=dim(dat)
  n=dat.dim[1]; p1=dat.dim[2]; p2=dat.dim[3]
  p=p1*p2
  
  S=mcov(dat)
  S=covKCD(S,p1,p2)$C
  
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


spike.test=function(dat,n.spike,level){
  
  dat.dim=dim(dat)
  n=dat.dim[1]; p1=dat.dim[2]; p2=dat.dim[3]
  p=p1*p2

  y=p/n
  c=1+sqrt(y)
  
  S=mcov(dat)
  S.core=covKCD(S,p1,p2)$C 
  
  S.core.eig=Re(eigen(S.core)$values)
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
  
  S=mcov(dat)
  S.core=covKCD(S,p1,p2)$C 
  
  S.core.spike=eigs_sym(S.core,n.spike)$values
  test.stat=S.core.spike[n.spike]/((p-sum(S.core.spike))/(p-n.spike))
  
  test.stat.boot=NULL
  
  for(i in 1:boot.size){
    
    samp=sample(1:n,n,replace=TRUE)
    temp=dat[samp,,]
    S.temp=mcov(temp)
    S.core.temp=covKCD(S.temp,p1,p2)$C 
    
    temp.spike=eigs_sym(S.core.temp,n.spike)$values
    test.stat.boot=c(test.stat.boot,temp.spike[n.spike]/((p-sum(temp.spike))/(p-n.spike)))
    
    if(i%%10==0){
      print(paste(i,"th iteration",sep=""))
    }
  }
  
  return((sum(test.stat.boot<=test.stat)+1)/(boot.size+1))
  
}

######## data loading 

load("./real_data/dat.train.RData")
load("./real_data/dat.test.RData")

dat.array=array(0,dim=c(9298,16,16))
dat.array[1:7291,,]=dat.train$data
dat.array[7292:9298,,]=dat.test$data

dat.ind=c(dat.train$target,dat.test$target)

ind=list()
for(i in 0:9){
  ind[[i+1]]=which(dat.ind==i)
}

set.seed(4200)

ind.samp=list()

for(i in 0:9){
  ind.samp[[i+1]]=sort(sample(ind[[i+1]],200))
}

core.array=array(0,dim=c(10,256,256))

for(i in 1:10){
  temp=mcov(dat.array[ind.samp[[i]],,])
  temp.core=covKCD(temp,16,16)$C
  core.array[i,,]=temp.core
}

thres=thres.C(256,200,2000,1,0.02)
thres
# thres 
# 0.3690731 

n.spike=NULL
for(i in 1:10){
  temp=dat.array[ind.samp[[i]],,]
  n.spike=c(n.spike,num.spike1(temp,thres))
}
n.spike
# 7; 13; 5; 6; 8; 5; 7; 15; 8; 12

# i=0
class0=list()
k=0
for(r in 5:9){
  k=k+1
  temp=dat.array[ind.samp[[1]],,]
  class0[[k]]=spike.test(temp,r,0.05)
}

# i=1
class1=list()
k=0
for(r in 11:15){
  k=k+1
  temp=dat.array[ind.samp[[2]],,]
  class1[[k]]=spike.test(temp,r,0.05)
}

# i=2
class2=list()
k=0
for(r in 3:7){
  k=k+1
  temp=dat.array[ind.samp[[3]],,]
  class2[[k]]=spike.test(temp,r,0.05)
}

# i=3
class3=list()
k=0
for(r in 4:8){
  k=k+1
  temp=dat.array[ind.samp[[4]],,]
  class3[[k]]=spike.test(temp,r,0.05)
}

# i=4
class4=list()
k=0
for(r in 6:10){
  k=k+1
  temp=dat.array[ind.samp[[5]],,]
  class4[[k]]=spike.test(temp,r,0.05)
}

# i=5
class5=list()
k=0
for(r in 3:7){
  k=k+1
  temp=dat.array[ind.samp[[6]],,]
  class5[[k]]=spike.test(temp,r,0.05)
}

# i=6
class6=list()
k=0
for(r in 5:9){
  k=k+1
  temp=dat.array[ind.samp[[7]],,]
  class6[[k]]=spike.test(temp,r,0.05)
}

# i=7
class7=list()
k=0
for(r in 13:17){
  k=k+1
  temp=dat.array[ind.samp[[8]],,]
  class7[[k]]=spike.test(temp,r,0.05)
}

# i=8
class8=list()
k=0
for(r in 6:10){
  k=k+1
  temp=dat.array[ind.samp[[9]],,]
  class8[[k]]=spike.test(temp,r,0.05)
}

# i=9
class9=list()
k=0
for(r in 8:14){
  k=k+1
  temp=dat.array[ind.samp[[10]],,]
  class9[[k]]=spike.test(temp,r,0.05)
}

save(class0,file="./real_data/class0.RData")
save(class1,file="./real_data/class1.RData")
save(class2,file="./real_data/class2.RData")
save(class3,file="./real_data/class3.RData")
save(class4,file="./real_data/class4.RData")
save(class5,file="./real_data/class5.RData")
save(class6,file="./real_data/class6.RData")
save(class7,file="./real_data/class7.RData")
save(class8,file="./real_data/class8.RData")
save(class9,file="./real_data/class9.RData")

#### Bootstrap p-value

pval=NULL
n.spike=c(7,13,5,6,8,5,7,15,8,12)
for(i in c(1,3:10)){
  temp=dat.array[ind.samp[[i]],,]
  temp=spike.pval(temp,n.spike[i],10000,1)
  pval=c(pval,temp)
  print(paste(i,"th class done"))
}
save(pval,file="./real_data/pval.RData")

#### Plots of eigenspectrum 

par(mfrow=c(2,2))

plot(Re(eigen(core.array[1,,])$values),xlab="Index",ylab="Eigenvalues",main="Eigenvalue of Sample Core for Digit 0")
abline(v=7,col="blue",lty=2)
plot(Re(eigen(core.array[2,,])$values),xlab="Index",ylab="Eigenvalues",main="Eigenvalue of Sample Core for Digit 1")
abline(v=13,col="blue",lty=2)
plot(Re(eigen(core.array[3,,])$values),xlab="Index",ylab="Eigenvalues",main="Eigenvalue of Sample Core for Digit 2")
abline(v=5,col="blue",lty=2)
plot(Re(eigen(core.array[4,,])$values),xlab="Index",ylab="Eigenvalues",main="Eigenvalue of Sample Core for Digit 3")
abline(v=6,col="blue",lty=2)

par(mfrow=c(2,2))
plot(Re(eigen(core.array[5,,])$values),xlab="Index",ylab="Eigenvalues",main="Eigenvalue of Sample Core for Digit 4")
abline(v=8,col="blue",lty=2)
plot(Re(eigen(core.array[6,,])$values),xlab="Index",ylab="Eigenvalues",main="Eigenvalue of Sample Core for Digit 5")
abline(v=5,col="blue",lty=2)
plot(Re(eigen(core.array[7,,])$values),xlab="Index",ylab="Eigenvalues",main="Eigenvalue of Sample Core for Digit 6")
abline(v=7,col="blue",lty=2)
plot(Re(eigen(core.array[8,,])$values),xlab="Index",ylab="Eigenvalues",main="Eigenvalue of Sample Core for Digit 7")
abline(v=15,col="blue",lty=2)

par(mfrow=c(1,2))
plot(Re(eigen(core.array[9,,])$values),xlab="Index",ylab="Eigenvalues",main="Eigenvalue of Sample Core for Digit 8")
abline(v=8,col="blue",lty=2)
plot(Re(eigen(core.array[10,,])$values),xlab="Index",ylab="Eigenvalues",main="Eigenvalue of Sample Core for Digit 9")
abline(v=12,col="blue",lty=2)


