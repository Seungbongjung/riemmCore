################ utensils

library(tidyverse)
library(ggpubr)

load("./numerical.simul/pi.core.fit.16.10.4.20.rds")
load("./numerical.simul/pi.core.fit.16.10.4.40.rds")
load("./numerical.simul/pi.core.fit.16.10.4.80.rds")
load("./numerical.simul/pi.core.fit.16.10.4.160.rds")
load("./numerical.simul/pi.core.fit.16.10.4.320.rds")
load("./numerical.simul/pi.core.fit.12.8.4.12.rds")
load("./numerical.simul/pi.core.fit.12.8.4.24.rds")
load("./numerical.simul/pi.core.fit.12.8.4.48.rds")
load("./numerical.simul/pi.core.fit.12.8.4.96.rds")
load("./numerical.simul/pi.core.fit.12.8.4.192.rds")

stat.analy=function(fit){
  
  kro.mle=fit$kro.mle
  cse=fit$cse
  init.ai=fit$pi.core.init.ai
  init.chol=fit$pi.core.init.chol
  fit.ai=fit$pi.core.ai
  fit.chol=fit$pi.core.chol
  
  kro.mle=cbind(init.ai[,1],kro.mle)
  cse=cbind(init.ai[,1],cse)
  init.chol=cbind(init.ai[,1],init.chol)
  
  colnames(kro.mle)=colnames(cse)=colnames(init.chol)=c("K","C","Sigma")
  
  fit.ai.analy=NULL; fit.chol.analy=NULL
  for(i in 1:100){
    fit.ai.analy=rbind(fit.ai.analy,fit.ai[[i]]$eval)
    fit.chol.analy=rbind(fit.chol.analy,fit.chol[[i]]$eval)
  }
  return(list(kro.mle=kro.mle,cse=cse,init.ai=init.ai,init.chol=init.chol,fit.ai.analy=fit.ai.analy,fit.chol.analy=fit.chol.analy))
}

conv.analy=function(fit1,fit2,fit3,fit4,fit5,name){
  
  fit.ai1=fit1$pi.core.ai
  fit.chol1=fit1$pi.core.chol
  fit.ai2=fit2$pi.core.ai
  fit.chol2=fit2$pi.core.chol
  fit.ai3=fit3$pi.core.ai
  fit.chol3=fit3$pi.core.chol
  fit.ai4=fit4$pi.core.ai
  fit.chol4=fit4$pi.core.chol
  fit.ai5=fit5$pi.core.ai
  fit.chol5=fit5$pi.core.chol
  
  rel.conv.ai=NULL
  conv.ai=NULL
  rel.conv.chol=NULL
  conv.chol=NULL
  
  for(i in 1:100){
    temp1=c(fit.ai1[[i]]$rel.conv,fit.ai2[[i]]$rel.conv,fit.ai3[[i]]$rel.conv,fit.ai4[[i]]$rel.conv,fit.ai5[[i]]$rel.conv)
    temp2=c(fit.ai1[[i]]$convergence,fit.ai2[[i]]$convergence,fit.ai3[[i]]$convergence,fit.ai4[[i]]$convergence,fit.ai5[[i]]$convergence)
    rel.conv.ai=cbind(rel.conv.ai,temp1)
    conv.ai=cbind(conv.ai,temp2)
    
    temp1=c(fit.chol1[[i]]$rel.conv,fit.chol2[[i]]$rel.conv,fit.chol3[[i]]$rel.conv,fit.chol4[[i]]$rel.conv,fit.chol5[[i]]$rel.conv)
    temp2=c(fit.chol1[[i]]$convergence,fit.chol2[[i]]$convergence,fit.chol3[[i]]$convergence,fit.chol4[[i]]$convergence,fit.chol5[[i]]$convergence)
    rel.conv.chol=cbind(rel.conv.chol,temp1)
    conv.chol=cbind(conv.chol,temp2)
  }
  rownames(rel.conv.ai)=rownames(conv.ai)=rownames(rel.conv.chol)=rownames(conv.chol)=name
  
  return(list(rel.conv.ai=rel.conv.ai,conv.ai=conv.ai,rel.conv.chol=rel.conv.chol,conv.chol=conv.chol))
}


################ analysis : mean & sd

fit.12.8.4.12=stat.analy(pi.core.fit.12.8.4.12)
colMeans(fit.12.8.4.12$kro.mle)
# 0.3786591 0.9642735 0.9100132
colMeans(fit.12.8.4.12$cse)
# 0.3786591 0.6069784 0.6678124 
colMeans(fit.12.8.4.12$init.ai)
# 0.3786591 0.5960378 0.7617314 
colMeans(fit.12.8.4.12$init.chol)
# 0.3786591 0.6265214 0.7617314 
colMeans(fit.12.8.4.12$fit.ai.analy)
# 0.3017722 0.5410270 0.5982598
colMeans(fit.12.8.4.12$fit.chol.analy)
# 0.3002198 0.5552815 0.586933

apply(fit.12.8.4.12$kro.mle,2,sd)
# 0.13255332 0.00000000 0.01292454 
apply(fit.12.8.4.12$cse,2,sd)
# 0.13255332 0.08911778 0.12833524 
apply(fit.12.8.4.12$init.ai,2,sd)
# 0.1325533 0.1001271 0.2519618  
apply(fit.12.8.4.12$init.chol,2,sd)
#0.1325533 0.1033793 0.2519618 
apply(fit.12.8.4.12$fit.ai.analy,2,sd)
#0.08116978 0.07364883 0.12791956
apply(fit.12.8.4.12$fit.chol.analy,2,sd)
#0.08177394 0.07829873 0.12402712 

fit.12.8.4.24=stat.analy(pi.core.fit.12.8.4.24)
colMeans(fit.12.8.4.24$kro.mle)
# 0.2415698 0.9642735 0.9097429
colMeans(fit.12.8.4.24$cse)
# 0.2415698 0.5269639 0.5339689
colMeans(fit.12.8.4.24$init.ai)
# 0.2415698 0.3924712 0.4903909
colMeans(fit.12.8.4.24$init.chol)
# 0.2415698 0.4122659 0.4903909
colMeans(fit.12.8.4.24$fit.ai.analy)
# 0.1976207 0.3124212 0.3524997
colMeans(fit.12.8.4.24$fit.chol.analy)
# 0.1974966 0.3239209 0.3519779 

apply(fit.12.8.4.24$kro.mle,2,sd)
#0.060482790 0.000000000 0.006904542 
apply(fit.12.8.4.24$cse,2,sd)
#0.06048279 0.05153313 0.06693021
apply(fit.12.8.4.24$init.ai,2,sd)
#0.06048279 0.04910626 0.10572909 
apply(fit.12.8.4.24$init.chol,2,sd)
#0.06048279 0.05089842 0.10572909
apply(fit.12.8.4.24$fit.ai.analy,2,sd)
#0.04532423 0.02912423 0.05511835 
apply(fit.12.8.4.24$fit.chol.analy,2,sd)
#0.04452993 0.03214142 0.05522742

fit.12.8.4.48=stat.analy(pi.core.fit.12.8.4.48)
colMeans(fit.12.8.4.48$kro.mle)
# 0.1619736 0.9642735 0.9104522
colMeans(fit.12.8.4.48$cse)
# 0.1619736 0.4601483 0.4766726 
colMeans(fit.12.8.4.48$init.ai)
# 0.1619736 0.2793752 0.3537015 
colMeans(fit.12.8.4.48$init.chol)
# 0.1619736 0.2929055 0.3537015 
colMeans(fit.12.8.4.48$fit.ai.analy)
# 0.1372027 0.2059765 0.2483067
colMeans(fit.12.8.4.48$fit.chol.analy)
# 0.1374848 0.2166128 0.2491791

apply(fit.12.8.4.48$kro.mle,2,sd)
#0.034689713 0.000000000 0.005470268
apply(fit.12.8.4.48$cse,2,sd)
#0.03468971 0.04651033 0.06355936 
apply(fit.12.8.4.48$init.ai,2,sd)
#0.03468971 0.03975086 0.06854848 
apply(fit.12.8.4.48$init.chol,2,sd)
#0.03468971 0.04120695 0.06854848 
apply(fit.12.8.4.48$fit.ai.analy,2,sd)
#0.02973382 0.01847062 0.03774694
apply(fit.12.8.4.48$fit.chol.analy,2,sd)
#0.02953145 0.02238409 0.04055895

fit.12.8.4.96=stat.analy(pi.core.fit.12.8.4.96)
colMeans(fit.12.8.4.96$kro.mle)
# 0.1370463 0.9449166 0.8338917
colMeans(fit.12.8.4.96$cse)
#0.1370463 0.3739768 0.3642135
colMeans(fit.12.8.4.96$init.ai)
#0.1370463 0.2307708 0.2778840 
colMeans(fit.12.8.4.96$init.chol)
# 0.1370463 0.2413310 0.2778840 
colMeans(fit.12.8.4.96$fit.ai.analy)
#0.1174490 0.1825283 0.2114893 
colMeans(fit.12.8.4.96$fit.chol.analy)
#0.1179751 0.1901567 0.2109889 

apply(fit.12.8.4.96$kro.mle,2,sd)
#0.034867514 0.000000000 0.008424431
apply(fit.12.8.4.96$cse,2,sd)
#0.03486751 0.04524801 0.06732658
apply(fit.12.8.4.96$init.ai,2,sd)
#0.03486751 0.02917023 0.05000803
apply(fit.12.8.4.96$init.chol,2,sd)
#0.03486751 0.02885691 0.05000803 
apply(fit.12.8.4.96$fit.ai.analy,2,sd)
#0.02608062 0.02051378 0.03166960
apply(fit.12.8.4.96$fit.chol.analy,2,sd)
#0.02663617 0.01896941 0.03172677

fit.12.8.4.192=stat.analy(pi.core.fit.12.8.4.192)
colMeans(fit.12.8.4.192$kro.mle)
#0.09856144 0.94491661 0.83332669
colMeans(fit.12.8.4.192$cse)
#0.09856144 0.27477500 0.26962505
colMeans(fit.12.8.4.192$init.ai)
#0.09856144 0.15931361 0.19421087 
colMeans(fit.12.8.4.192$init.chol)
#0.09856144 0.16618407 0.19421087
colMeans(fit.12.8.4.192$fit.ai.analy)
#0.08594988 0.12656208 0.15071566
colMeans(fit.12.8.4.192$fit.chol.analy)
#0.08635607 0.13287770 0.15100579

apply(fit.12.8.4.192$kro.mle,2,sd)
#0.026003586 0.000000000 0.007177489 
apply(fit.12.8.4.192$cse,2,sd)
#0.02600359 0.03406020 0.04661363 
apply(fit.12.8.4.192$init.ai,2,sd)
#0.02600359 0.02223817 0.04398408 
apply(fit.12.8.4.192$init.chol,2,sd)
#0.02600359 0.02210112 0.04398408
apply(fit.12.8.4.192$fit.ai.analy,2,sd)
#0.02155511 0.01346157 0.02734832 
apply(fit.12.8.4.192$fit.chol.analy,2,sd)
#0.02153682 0.01375523 0.02785751

################ analysis : convergence

conv.16.10=conv.analy(pi.core.fit.16.10.4.20,pi.core.fit.16.10.4.40,pi.core.fit.16.10.4.80,
                         pi.core.fit.16.10.4.160,pi.core.fit.16.10.4.320,c("20","40","80","160","320"))

apply(conv.16.10$rel.conv.ai,1,mean)
# 0.92 0.98 0.98 1.00 1.00
apply(conv.16.10$conv.ai,1,mean)
# 1   1   1   1   1 
apply(conv.16.10$rel.conv.ai,1,sd)
# 0.2726599 0.1407053 0.1407053 0.0000000 0.0000000 
apply(conv.16.10$conv.ai,1,sd)
# 0 0 0 0 0

apply(conv.16.10$rel.conv.chol,1,mean)
# 0.97 0.98 0.99 0.99 0.99 
apply(conv.16.10$conv.chol,1,mean)
# 1   1   1   1   1 
apply(conv.16.10$rel.conv.chol,1,sd)
# 0.1714466 0.1407053 0.1000000 0.1000000 0.1000000 
apply(conv.16.10$conv.chol,1,sd)
# 0 0 0 0 0 

conv.12.8=conv.analy(pi.core.fit.12.8.4.12,pi.core.fit.12.8.4.24,pi.core.fit.12.8.4.48,
                      pi.core.fit.12.8.4.96,pi.core.fit.12.8.4.192,c("12","24","48","96","192"))

apply(conv.12.8$rel.conv.ai,1,mean)
# 0.92 0.98 0.98 1.00 1.00
apply(conv.12.8$conv.ai,1,mean)
# 1   1   1   1   1 
apply(conv.12.8$rel.conv.ai,1,sd)
# 0.2726599 0.1407053 0.1407053 0.0000000 0.0000000 
apply(conv.12.8$conv.ai,1,sd)
# 0 0 0 0 0

apply(conv.12.8$rel.conv.chol,1,mean)
# 0.97 0.98 0.99 0.99 0.99 
apply(conv.12.8$conv.chol,1,mean)
# 1   1   1   1   1 
apply(conv.12.8$rel.conv.chol,1,sd)
# 0.1714466 0.1407053 0.1000000 0.1000000 0.1000000 
apply(conv.12.8$conv.chol,1,sd)
# 0 0 0 0 0 

################ analysis : convergence