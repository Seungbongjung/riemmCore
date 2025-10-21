################ utensils

library(tidyverse)
library(latex2exp)

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

plot.K.analy=function(fit1,fit2,fit3,fit4,fit5,samp.size){

  method=rep(c("Base-AI","PICSE-AI","PICSE-Chol"),each=100)
  
  dat1=tibble(stat=c(fit1$init.ai[,1],fit1$fit.ai.analy[,1],fit1$fit.chol.analy[,1]),Method=method,samp.size=rep(samp.size[1],300))
  dat2=tibble(stat=c(fit2$init.ai[,1],fit2$fit.ai.analy[,1],fit2$fit.chol.analy[,1]),Method=method,samp.size=rep(samp.size[2],300))
  dat3=tibble(stat=c(fit3$init.ai[,1],fit3$fit.ai.analy[,1],fit3$fit.chol.analy[,1]),Method=method,samp.size=rep(samp.size[3],300))
  dat4=tibble(stat=c(fit4$init.ai[,1],fit4$fit.ai.analy[,1],fit4$fit.chol.analy[,1]),Method=method,samp.size=rep(samp.size[4],300))
  dat5=tibble(stat=c(fit5$init.ai[,1],fit5$fit.ai.analy[,1],fit5$fit.chol.analy[,1]),Method=method,samp.size=rep(samp.size[5],300))
  
  dat=rbind(dat1,dat2,dat3,dat4,dat5)
  dat$samp.size=factor(dat$samp.size,levels=samp.size)
  dat$Method=factor(dat$Method,levels=c("Base-AI","PICSE-AI","PICSE-Chol"))
  return(dat)
}

plot.C.analy=function(fit1,fit2,fit3,fit4,fit5,samp.size){
  
  method=rep(c("KMLE","CSE","Base-AI","BASE-Chol","PICSE-AI","PICSE-Chol"),each=100)
  
  dat1=tibble(stat=c(fit1$kro.mle[,2],fit1$cse[,2],fit1$init.ai[,2],fit1$init.chol[,2],
                     fit1$fit.ai.analy[,2],fit1$fit.chol.analy[,2]),Method=method,samp.size=rep(samp.size[1],600))
  dat2=tibble(stat=c(fit2$kro.mle[,2],fit2$cse[,2],fit2$init.ai[,2],fit2$init.chol[,2],
                     fit2$fit.ai.analy[,2],fit2$fit.chol.analy[,2]),Method=method,samp.size=rep(samp.size[2],600))
  dat3=tibble(stat=c(fit3$kro.mle[,2],fit3$cse[,2],fit3$init.ai[,2],fit3$init.chol[,2],
                     fit3$fit.ai.analy[,2],fit3$fit.chol.analy[,2]),Method=method,samp.size=rep(samp.size[3],600))
  dat4=tibble(stat=c(fit4$kro.mle[,2],fit4$cse[,2],fit4$init.ai[,2],fit4$init.chol[,2],
                     fit4$fit.ai.analy[,2],fit4$fit.chol.analy[,2]),Method=method,samp.size=rep(samp.size[4],600))
  dat5=tibble(stat=c(fit5$kro.mle[,2],fit5$cse[,2],fit5$init.ai[,2],fit5$init.chol[,2],
                     fit5$fit.ai.analy[,2],fit5$fit.chol.analy[,2]),Method=method,samp.size=rep(samp.size[5],600))
  
  dat=rbind(dat1,dat2,dat3,dat4,dat5)
  dat$samp.size=factor(dat$samp.size,levels=samp.size)
  dat$Method=factor(dat$Method,levels=c("KMLE","CSE","Base-AI","BASE-Chol","PICSE-AI","PICSE-Chol"))
  return(dat)
  
}

plot.Sigma.analy=function(fit1,fit2,fit3,fit4,fit5,samp.size){
  
  method=rep(c("KMLE","CSE","Base-AI","PICSE-AI","PICSE-Chol"),each=100)
  
  dat1=tibble(stat=c(fit1$kro.mle[,3],fit1$cse[,3],fit1$init.ai[,3],
                     fit1$fit.ai.analy[,3],fit1$fit.chol.analy[,3]),Method=method,samp.size=rep(samp.size[1],500))
  dat2=tibble(stat=c(fit2$kro.mle[,3],fit2$cse[,3],fit2$init.ai[,3],
                     fit2$fit.ai.analy[,3],fit2$fit.chol.analy[,3]),Method=method,samp.size=rep(samp.size[2],500))
  dat3=tibble(stat=c(fit3$kro.mle[,3],fit3$cse[,3],fit3$init.ai[,3],
                     fit3$fit.ai.analy[,3],fit3$fit.chol.analy[,3]),Method=method,samp.size=rep(samp.size[3],500))
  dat4=tibble(stat=c(fit4$kro.mle[,3],fit4$cse[,3],fit4$init.ai[,3],
                     fit4$fit.ai.analy[,3],fit4$fit.chol.analy[,3]),Method=method,samp.size=rep(samp.size[4],500))
  dat5=tibble(stat=c(fit5$kro.mle[,3],fit5$cse[,3],fit5$init.ai[,3],
                     fit5$fit.ai.analy[,3],fit5$fit.chol.analy[,3]),Method=method,samp.size=rep(samp.size[5],500))
  
  dat=rbind(dat1,dat2,dat3,dat4,dat5)
  dat$samp.size=factor(dat$samp.size,levels=samp.size)
  dat$Method=factor(dat$Method,levels=c("KMLE","CSE","Base-AI","PICSE-AI","PICSE-Chol"))
  return(dat)
  
  
}

################ analysis : mean & sd

fit.16.10.4.20=stat.analy(pi.core.fit.16.10.4.20)
colMeans(fit.16.10.4.20$kro.mle)
# 0.3786591 0.9642735 0.9100132
colMeans(fit.16.10.4.20$cse)
# 0.3786591 0.6069784 0.6678124 
colMeans(fit.16.10.4.20$init.ai)
# 0.3786591 0.5960378 0.7617314 
colMeans(fit.16.10.4.20$init.chol)
# 0.3786591 0.6265214 0.7617314 
colMeans(fit.16.10.4.20$fit.ai.analy)
# 0.3017722 0.5410270 0.5982598
colMeans(fit.16.10.4.20$fit.chol.analy)
# 0.3002198 0.5552815 0.586933

apply(fit.16.10.4.20$kro.mle,2,sd)
# 0.13255332 0.00000000 0.01292454 
apply(fit.16.10.4.20$cse,2,sd)
# 0.13255332 0.08911778 0.12833524 
apply(fit.16.10.4.20$init.ai,2,sd)
# 0.1325533 0.1001271 0.2519618  
apply(fit.16.10.4.20$init.chol,2,sd)
#0.1325533 0.1033793 0.2519618 
apply(fit.16.10.4.20$fit.ai.analy,2,sd)
#0.08116978 0.07364883 0.12791956
apply(fit.16.10.4.20$fit.chol.analy,2,sd)
#0.08177394 0.07829873 0.12402712 

fit.16.10.4.40=stat.analy(pi.core.fit.16.10.4.40)
colMeans(fit.16.10.4.40$kro.mle)
# 0.2415698 0.9642735 0.9097429
colMeans(fit.16.10.4.40$cse)
# 0.2415698 0.5269639 0.5339689
colMeans(fit.16.10.4.40$init.ai)
# 0.2415698 0.3924712 0.4903909
colMeans(fit.16.10.4.40$init.chol)
# 0.2415698 0.4122659 0.4903909
colMeans(fit.16.10.4.40$fit.ai.analy)
# 0.1976207 0.3124212 0.3524997
colMeans(fit.16.10.4.40$fit.chol.analy)
# 0.1974966 0.3239209 0.3519779 

apply(fit.16.10.4.40$kro.mle,2,sd)
#0.060482790 0.000000000 0.006904542 
apply(fit.16.10.4.40$cse,2,sd)
#0.06048279 0.05153313 0.06693021
apply(fit.16.10.4.40$init.ai,2,sd)
#0.06048279 0.04910626 0.10572909 
apply(fit.16.10.4.40$init.chol,2,sd)
#0.06048279 0.05089842 0.10572909
apply(fit.16.10.4.40$fit.ai.analy,2,sd)
#0.04532423 0.02912423 0.05511835 
apply(fit.16.10.4.40$fit.chol.analy,2,sd)
#0.04452993 0.03214142 0.05522742

fit.16.10.4.80=stat.analy(pi.core.fit.16.10.4.80)
colMeans(fit.16.10.4.80$kro.mle)
# 0.1619736 0.9642735 0.9104522
colMeans(fit.16.10.4.80$cse)
# 0.1619736 0.4601483 0.4766726 
colMeans(fit.16.10.4.80$init.ai)
# 0.1619736 0.2793752 0.3537015 
colMeans(fit.16.10.4.80$init.chol)
# 0.1619736 0.2929055 0.3537015 
colMeans(fit.16.10.4.80$fit.ai.analy)
# 0.1372027 0.2059765 0.2483067
colMeans(fit.16.10.4.80$fit.chol.analy)
# 0.1378048 0.2166128 0.2491791

apply(fit.16.10.4.80$kro.mle,2,sd)
#0.034689713 0.000000000 0.005470268
apply(fit.16.10.4.80$cse,2,sd)
#0.03468971 0.04651033 0.06355936 
apply(fit.16.10.4.80$init.ai,2,sd)
#0.03468971 0.03975086 0.06858048 
apply(fit.16.10.4.80$init.chol,2,sd)
#0.03468971 0.04120695 0.06858048 
apply(fit.16.10.4.80$fit.ai.analy,2,sd)
#0.02973382 0.01847062 0.03774694
apply(fit.16.10.4.80$fit.chol.analy,2,sd)
#0.02953145 0.02238409 0.04055895

fit.16.10.4.160=stat.analy(pi.core.fit.16.10.4.160)
colMeans(fit.16.10.4.160$kro.mle)
#0.1173761 0.9642735 0.9104901
colMeans(fit.16.10.4.160$cse)
#0.1173761 0.3646184 0.3773209 
colMeans(fit.16.10.4.160$init.ai)
#0.1173761 0.1917740 0.2445652
colMeans(fit.16.10.4.160$init.chol)
#0.1173761 0.2010542 0.2445652 
colMeans(fit.16.10.4.160$fit.ai.analy)
#0.09779566 0.13898238 0.17324516
colMeans(fit.16.10.4.160$fit.chol.analy)
#0.09853161 0.14649609 0.17411761

apply(fit.16.10.4.160$kro.mle,2,sd)
#0.021367066 0.000000000 0.003757931 
apply(fit.16.10.4.160$cse,2,sd)
#0.02136707 0.03636015 0.05286314
apply(fit.16.10.4.160$init.ai,2,sd)
#0.02136707 0.02587712 0.04493711
apply(fit.16.10.4.160$init.chol,2,sd)
#0.02136707 0.02621295 0.04493711
apply(fit.16.10.4.160$fit.ai.analy,2,sd)
#0.02180036 0.01304000 0.02704266
apply(fit.16.10.4.160$fit.chol.analy,2,sd)
#0.02186110 0.01354360 0.02762492

fit.16.10.4.320=stat.analy(pi.core.fit.16.10.4.320)
colMeans(fit.16.10.4.320$kro.mle)
#0.08055973 0.96427346 0.91003987
colMeans(fit.16.10.4.320$cse)
#0.08055973 0.27130863 0.27922741 
colMeans(fit.16.10.4.320$init.ai)
#0.08055973 0.13247543 0.16893904 
colMeans(fit.16.10.4.320$init.chol)
#0.08055973 0.13812545 0.16893904
colMeans(fit.16.10.4.320$fit.ai.analy)
#0.06890777 0.09458343 0.11991833 
colMeans(fit.16.10.4.320$fit.chol.analy)
#0.06980748 0.10102651 0.12270935 

apply(fit.16.10.4.320$kro.mle,2,sd)
#0.016040112 0.000000000 0.002707311
apply(fit.16.10.4.320$cse,2,sd)
#0.01604011 0.02868382 0.03939750
apply(fit.16.10.4.320$init.ai,2,sd)
#0.01604011 0.01616103 0.03099864
apply(fit.16.10.4.320$init.chol,2,sd)
#0.01604011 0.01609051 0.03099864
apply(fit.16.10.4.320$fit.ai.analy,2,sd)
#0.014152782 0.008587168 0.019824027 
apply(fit.16.10.4.320$fit.chol.analy,2,sd)
#0.014201314 0.009187079 0.021035087

fit.12.8.4.12=stat.analy(pi.core.fit.12.8.4.12)
colMeans(fit.12.8.4.12$kro.mle)
# 0.4522738 0.9449166 0.8367261 
colMeans(fit.12.8.4.12$cse)
# 0.4522738 0.6395022 0.7708857  
colMeans(fit.12.8.4.12$init.ai)
#0.4522738 0.7740625 0.8778848
colMeans(fit.12.8.4.12$init.chol)
#0.4522738 0.8150695 0.8778848 
colMeans(fit.12.8.4.12$fit.ai.analy)
#0.3848148 0.8554445 0.8472965 
colMeans(fit.12.8.4.12$fit.chol.analy)
#0.3796479 0.8701837 0.8326581 

apply(fit.12.8.4.12$kro.mle,2,sd)
# 0.14106896 0.00000000 0.02669203
apply(fit.12.8.4.12$cse,2,sd)
#0.14106896 0.07754842 0.20549847 
apply(fit.12.8.4.12$init.ai,2,sd)
#0.14106896 0.09467952 0.26276207 
apply(fit.12.8.4.12$init.chol,2,sd)
#0.14106896 0.09812882 0.26276207 
apply(fit.12.8.4.12$fit.ai.analy,2,sd)
#0.08760867 0.15095545 0.23450719
apply(fit.12.8.4.12$fit.chol.analy,2,sd)
#0.08521818 0.17099171 0.25583249

fit.12.8.4.24=stat.analy(pi.core.fit.12.8.4.24)
colMeans(fit.12.8.4.24$kro.mle)
# 0.3039174 0.9449166 0.8338129
colMeans(fit.12.8.4.24$cse)
#0.3039174 0.5559627 0.5419085 
colMeans(fit.12.8.4.24$init.ai)
#0.3039174 0.5022790 0.6122633 
colMeans(fit.12.8.4.24$init.chol)
#0.3039174 0.5263172 0.6122633 
colMeans(fit.12.8.4.24$fit.ai.analy)
#0.2555481 0.4531483 0.4823386
colMeans(fit.12.8.4.24$fit.chol.analy)
#0.2530182 0.4637133 0.4767197 

apply(fit.12.8.4.24$kro.mle,2,sd)
#0.09125345 0.00000000 0.02031843 
apply(fit.12.8.4.24$cse,2,sd)
#0.09125345 0.07566858 0.09946387
apply(fit.12.8.4.24$init.ai,2,sd)
#0.09125345 0.07266116 0.17771548 
apply(fit.12.8.4.24$init.chol,2,sd)
#0.09125345 0.07379251 0.17771548
apply(fit.12.8.4.24$fit.ai.analy,2,sd)
#0.06400491 0.05069165 0.09845180 
apply(fit.12.8.4.24$fit.chol.analy,2,sd)
#0.06261652 0.04630585 0.09515542

fit.12.8.4.48=stat.analy(pi.core.fit.12.8.4.48)
colMeans(fit.12.8.4.48$kro.mle)
#0.2028973 0.9449166 0.8346989 
colMeans(fit.12.8.4.48$cse)
#0.2028973 0.4733454 0.4478934 
colMeans(fit.12.8.4.48$init.ai)
#0.2028973 0.3316726 0.3933783
colMeans(fit.12.8.4.48$init.chol)
#0.2028973 0.3469413 0.3933783 
colMeans(fit.12.8.4.48$fit.ai.analy)
#0.1748659 0.2771203 0.3120516
colMeans(fit.12.8.4.48$fit.chol.analy)
#0.1760349 0.2869399 0.3096419 

apply(fit.12.8.4.48$kro.mle,2,sd)
#0.05601882 0.00000000 0.01428471
apply(fit.12.8.4.48$cse,2,sd)
#0.05601882 0.05922803 0.08162166
apply(fit.12.8.4.48$init.ai,2,sd)
#0.05601882 0.04818168 0.07955302 
apply(fit.12.8.4.48$init.chol,2,sd)
#0.05601882 0.04805332 0.07955302 
apply(fit.12.8.4.48$fit.ai.analy,2,sd)
#0.04436547 0.03193375 0.05490327 
apply(fit.12.8.4.48$fit.chol.analy,2,sd)
#0.04631731 0.03150046 0.05446189 

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
#0.79 0.94 0.95 0.99 1.00 
apply(conv.12.8$conv.ai,1,mean)
# 1   1   1   1   1 
apply(conv.12.8$rel.conv.ai,1,sd)
#0.4093602 0.2386833 0.2190429 0.1000000 0.000000
apply(conv.12.8$conv.ai,1,sd)
# 0 0 0 0 0

apply(conv.12.8$rel.conv.chol,1,mean)
#0.83 0.98 0.99 1.00 1.00 
apply(conv.12.8$conv.chol,1,mean)
# 1   1   1   1   1 
apply(conv.12.8$rel.conv.chol,1,sd)
# 0.3775252 0.1407053 0.1000000 0.0000000 0.0000000 
apply(conv.12.8$conv.chol,1,sd)
# 0 0 0 0 0 

################ analysis : plot 

K.analy.16.10.4=plot.K.analy(fit.16.10.4.20,fit.16.10.4.40,fit.16.10.4.80,fit.16.10.4.160,fit.16.10.4.320,c("20","40","80","160","320"))
ggplot(K.analy.16.10.4)+geom_boxplot(mapping=aes(x=samp.size,y=stat,col=Method))+theme_bw()+xlab("Sample Size")+ylab("Relative Norm")+
  ggtitle(TeX("Separable Component Consistency : $(p_1,p_2,r)=(16,10,4)$"))+theme(plot.title=element_text(hjust=0.5))

C.analy.16.10.4=plot.C.analy(fit.16.10.4.20,fit.16.10.4.40,fit.16.10.4.80,fit.16.10.4.160,fit.16.10.4.320,c("20","40","80","160","320"))
ggplot(C.analy.16.10.4)+geom_boxplot(mapping=aes(x=samp.size,y=stat,col=Method))+theme_bw()+xlab("Sample Size")+ylab("Relative Norm")+
  ggtitle(TeX("Core Component Consistency : $(p_1,p_2,r)=(16,10,4)$"))+theme(plot.title=element_text(hjust=0.5))

Sigma.analy.16.10.4=plot.Sigma.analy(fit.16.10.4.20,fit.16.10.4.40,fit.16.10.4.80,fit.16.10.4.160,fit.16.10.4.320,c("20","40","80","160","320"))
ggplot(Sigma.analy.16.10.4)+geom_boxplot(mapping=aes(x=samp.size,y=stat,col=Method))+theme_bw()+xlab("Sample Size")+ylab("Relative Norm")+
  ggtitle(TeX("Sigma Consistency : $(p_1,p_2,r)=(16,10,4)$"))+theme(plot.title=element_text(hjust=0.5))

K.analy.12.8.4=plot.K.analy(fit.12.8.4.12,fit.12.8.4.24,fit.12.8.4.48,fit.12.8.4.96,fit.12.8.4.192,c("12","24","48","96","192"))
ggplot(K.analy.12.8.4)+geom_boxplot(mapping=aes(x=samp.size,y=stat,col=Method))+theme_bw()+xlab("Sample Size")+ylab("Relative Norm")+
  ggtitle(TeX("Separable Component Consistency : $(p_1,p_2,r)=(12,8,4)$"))+theme(plot.title=element_text(hjust=0.5))

C.analy.12.8.4=plot.C.analy(fit.12.8.4.12,fit.12.8.4.24,fit.12.8.4.48,fit.12.8.4.96,fit.12.8.4.192,c("12","24","48","96","192"))
ggplot(C.analy.12.8.4)+geom_boxplot(mapping=aes(x=samp.size,y=stat,col=Method))+theme_bw()+xlab("Sample Size")+ylab("Relative Norm")+
  ggtitle(TeX("Core Component Consistency : $(p_1,p_2,r)=(12,8,4)$"))+theme(plot.title=element_text(hjust=0.5))

Sigma.analy.12.8.4=plot.Sigma.analy(fit.12.8.4.12,fit.12.8.4.24,fit.12.8.4.48,fit.12.8.4.96,fit.12.8.4.192,c("12","24","48","96","192"))
ggplot(Sigma.analy.12.8.4)+geom_boxplot(mapping=aes(x=samp.size,y=stat,col=Method))+theme_bw()+xlab("Sample Size")+ylab("Relative Norm")+
  ggtitle(TeX("Sigma Consistency : $(p_1,p_2,r)=(12,8,4)$"))+theme(plot.title=element_text(hjust=0.5))
