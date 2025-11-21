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
# 0.3017697 0.5410282 0.5982650 
colMeans(fit.16.10.4.20$fit.chol.analy)
# 0.3002252 0.5553206 0.5869142 

apply(fit.16.10.4.20$kro.mle,2,sd)
# 0.13255332 0.00000000 0.01292454 
apply(fit.16.10.4.20$cse,2,sd)
# 0.13255332 0.08911778 0.12833524 
apply(fit.16.10.4.20$init.ai,2,sd)
# 0.1325533 0.1001271 0.2519618  
apply(fit.16.10.4.20$init.chol,2,sd)
#0.1325533 0.1033793 0.2519618 
apply(fit.16.10.4.20$fit.ai.analy,2,sd)
#0.08116240 0.07365024 0.12792416 
apply(fit.16.10.4.20$fit.chol.analy,2,sd)
#0.08177432 0.07832626 0.12402947 

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
# 0.1976192 0.3124459 0.3541230 
colMeans(fit.16.10.4.40$fit.chol.analy)
# 0.1974940 0.3239187 0.3519754 

apply(fit.16.10.4.40$kro.mle,2,sd)
#0.060482790 0.000000000 0.006904542 
apply(fit.16.10.4.40$cse,2,sd)
#0.06048279 0.05153313 0.06693021
apply(fit.16.10.4.40$init.ai,2,sd)
#0.06048279 0.04910626 0.10572909 
apply(fit.16.10.4.40$init.chol,2,sd)
#0.06048279 0.05089842 0.10572909
apply(fit.16.10.4.40$fit.ai.analy,2,sd)
#0.04534170 0.02912326 0.05512040 
apply(fit.16.10.4.40$fit.chol.analy,2,sd)
#0.04452966 0.03213925 0.05522503 

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
# 0.1370745 0.2062081 0.2482845 
colMeans(fit.16.10.4.80$fit.chol.analy)
# 0.1378079 0.2166137 0.2491809

apply(fit.16.10.4.80$kro.mle,2,sd)
#0.034689713 0.000000000 0.005470268
apply(fit.16.10.4.80$cse,2,sd)
#0.03468971 0.04651033 0.06355936 
apply(fit.16.10.4.80$init.ai,2,sd)
#0.03468971 0.03975086 0.06858048 
apply(fit.16.10.4.80$init.chol,2,sd)
#0.03468971 0.04120695 0.06858048 
apply(fit.16.10.4.80$fit.ai.analy,2,sd)
#0.02991936 0.01827500 0.03777477 
apply(fit.16.10.4.80$fit.chol.analy,2,sd)
#0.02953074 0.02238390 0.04055885

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
#0.09777752 0.13896509 0.17322089 
colMeans(fit.16.10.4.160$fit.chol.analy)
#0.09844996 0.14601696 0.17360286 

apply(fit.16.10.4.160$kro.mle,2,sd)
#0.021367066 0.000000000 0.003757931 
apply(fit.16.10.4.160$cse,2,sd)
#0.02136707 0.03636015 0.05286314
apply(fit.16.10.4.160$init.ai,2,sd)
#0.02136707 0.02587712 0.04493711
apply(fit.16.10.4.160$init.chol,2,sd)
#0.02136707 0.02621295 0.04493711
apply(fit.16.10.4.160$fit.ai.analy,2,sd)
#0.02181327 0.01303182 0.02705257 
apply(fit.16.10.4.160$fit.chol.analy,2,sd)
#0.02181127 0.01311251 0.02706567

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
#0.06890775 0.09458361 0.11991930 
colMeans(fit.16.10.4.320$fit.chol.analy)
#0.06980759 0.10103182 0.12270656  

apply(fit.16.10.4.320$kro.mle,2,sd)
#0.016040112 0.000000000 0.002707311
apply(fit.16.10.4.320$cse,2,sd)
#0.01604011 0.02868382 0.03939750
apply(fit.16.10.4.320$init.ai,2,sd)
#0.01604011 0.01616103 0.03099864
apply(fit.16.10.4.320$init.chol,2,sd)
#0.01604011 0.01609051 0.03099864
apply(fit.16.10.4.320$fit.ai.analy,2,sd)
#0.014152768 0.008587013 0.019824213 
apply(fit.16.10.4.320$fit.chol.analy,2,sd)
#0.014201094 0.009184269 0.021037265 

fit.12.8.4.12=stat.analy(pi.core.fit.12.8.4.12)
colMeans(fit.12.8.4.12$kro.mle)
#0.5155598 0.9392441 0.8724060 
colMeans(fit.12.8.4.12$cse)
#0.5155598 0.7032931 0.8031305
colMeans(fit.12.8.4.12$init.ai)
#0.5155598 0.8385021 0.9165219 
colMeans(fit.12.8.4.12$init.chol)
#0.5155598 0.8835863 0.9165219 
colMeans(fit.12.8.4.12$fit.ai.analy)
#0.4447631 0.9395184 0.8666650
colMeans(fit.12.8.4.12$fit.chol.analy)
#0.4386686 0.9599945 0.8524755 

apply(fit.12.8.4.12$kro.mle,2,sd)
#0.17906644 0.00000000 0.01925994
apply(fit.12.8.4.12$cse,2,sd)
#0.1790664 0.0857128 0.2049946 
apply(fit.12.8.4.12$init.ai,2,sd)
#0.1790664 0.1021148 0.2645013
apply(fit.12.8.4.12$init.chol,2,sd)
#0.1790664 0.1080899 0.2645013 
apply(fit.12.8.4.12$fit.ai.analy,2,sd)
#0.1209173 0.1497593 0.2148102 
apply(fit.12.8.4.12$fit.chol.analy,2,sd)
#0.1175681 0.1454069 0.2015733 

fit.12.8.4.24=stat.analy(pi.core.fit.12.8.4.24)
colMeans(fit.12.8.4.24$kro.mle)
#0.3394362 0.9392441 0.8705569 
colMeans(fit.12.8.4.24$cse)
#0.3394362 0.6082128 0.5797804 
colMeans(fit.12.8.4.24$init.ai)
#0.3394362 0.5419639 0.6013973
colMeans(fit.12.8.4.24$init.chol)
#0.3394362 0.5689567 0.6013973 
colMeans(fit.12.8.4.24$fit.ai.analy)
#0.2735300 0.4857656 0.4718893 
colMeans(fit.12.8.4.24$fit.chol.analy)
#0.2751602 0.5002325 0.4694646

apply(fit.12.8.4.24$kro.mle,2,sd)
#0.0995071 0.0000000 0.0153211 
apply(fit.12.8.4.24$cse,2,sd)
#0.09950710 0.06311715 0.09536213
apply(fit.12.8.4.24$init.ai,2,sd)
#0.09950710 0.07394252 0.16642627
apply(fit.12.8.4.24$init.chol,2,sd)
#0.09950710 0.07706086 0.16642627 
apply(fit.12.8.4.24$fit.ai.analy,2,sd)
#0.05662426 0.06186861 0.08631568 
apply(fit.12.8.4.24$fit.chol.analy,2,sd)
#0.05794195 0.06666551 0.09047187

fit.12.8.4.48=stat.analy(pi.core.fit.12.8.4.48)
colMeans(fit.12.8.4.48$kro.mle)
#0.2337188 0.9392441 0.8710452 
colMeans(fit.12.8.4.48$cse)
#0.2337188 0.5113725 0.4843611 
colMeans(fit.12.8.4.48$init.ai)
#0.2337188 0.3733815 0.4219599 
colMeans(fit.12.8.4.48$init.chol)
#0.2337188 0.3894310 0.4219599 
colMeans(fit.12.8.4.48$fit.ai.analy)
#0.1948336 0.3008409 0.3221433
colMeans(fit.12.8.4.48$fit.chol.analy)
#0.1956306 0.3122652 0.3227018 

apply(fit.12.8.4.48$kro.mle,2,sd)
#0.05826907 0.00000000 0.01144727
apply(fit.12.8.4.48$cse,2,sd)
#0.05826907 0.06185817 0.08891447 
apply(fit.12.8.4.48$init.ai,2,sd)
#0.05826907 0.05405061 0.09715245
apply(fit.12.8.4.48$init.chol,2,sd)
#0.05826907 0.05681609 0.09715245
apply(fit.12.8.4.48$fit.ai.analy,2,sd)
#0.04429187 0.03266255 0.06528784  
apply(fit.12.8.4.48$fit.chol.analy,2,sd)
#0.04345458 0.03302290 0.06444436

fit.12.8.4.96=stat.analy(pi.core.fit.12.8.4.96)
colMeans(fit.12.8.4.96$kro.mle)
#0.1606135 0.9392441 0.8693768 
colMeans(fit.12.8.4.96$cse)
#0.1606135 0.4023968 0.3848465
colMeans(fit.12.8.4.96$init.ai)
#0.1606135 0.2564021 0.3033561
colMeans(fit.12.8.4.96$init.chol)
#0.1606135 0.2674503 0.3033561  
colMeans(fit.12.8.4.96$fit.ai.analy)
#0.1371613 0.2025365 0.2257079 
colMeans(fit.12.8.4.96$fit.chol.analy)
#0.1375550 0.2105208 0.2267474 

apply(fit.12.8.4.96$kro.mle,2,sd)
#0.041772925 0.000000000 0.007709417
apply(fit.12.8.4.96$cse,2,sd)
#0.04177292 0.04788673 0.06494516
apply(fit.12.8.4.96$init.ai,2,sd)
#0.04177292 0.03889531 0.06918943
apply(fit.12.8.4.96$init.chol,2,sd)
#0.04177292 0.03990911 0.06918943 
apply(fit.12.8.4.96$fit.ai.analy,2,sd)
#0.03540282 0.02238946 0.03937123 
apply(fit.12.8.4.96$fit.chol.analy,2,sd)
#0.03574683 0.02029851 0.03924695

fit.12.8.4.192=stat.analy(pi.core.fit.12.8.4.192)
colMeans(fit.12.8.4.192$kro.mle)
#0.1128971 0.9392441 0.8690010 
colMeans(fit.12.8.4.192$cse)
#0.1128971 0.3014231 0.2870328
colMeans(fit.12.8.4.192$init.ai)
#0.1128971 0.1774998 0.2062094 
colMeans(fit.12.8.4.192$init.chol)
#0.1128971 0.1846181 0.2062094
colMeans(fit.12.8.4.192$fit.ai.analy)
#0.09797826 0.13637381 0.15633188 
colMeans(fit.12.8.4.192$fit.chol.analy)
#0.09866275 0.14367309 0.15782878 

apply(fit.12.8.4.192$kro.mle,2,sd)
#0.026569756 0.000000000 0.005467947 
apply(fit.12.8.4.192$cse,2,sd)
#0.02656976 0.03387858 0.04985670
apply(fit.12.8.4.192$init.ai,2,sd)
#0.02656976 0.02708878 0.04351865 
apply(fit.12.8.4.192$init.chol,2,sd)
#0.02656976 0.02816389 0.04351865 
apply(fit.12.8.4.192$fit.ai.analy,2,sd)
#0.02434107 0.01667196 0.02744835
apply(fit.12.8.4.192$fit.chol.analy,2,sd)
#0.02458359 0.01750181 0.02859850 

################ analysis : convergence

conv.16.10=conv.analy(pi.core.fit.16.10.4.20,pi.core.fit.16.10.4.40,pi.core.fit.16.10.4.80,
                         pi.core.fit.16.10.4.160,pi.core.fit.16.10.4.320,c("20","40","80","160","320"))

apply(conv.16.10$rel.conv.ai,1,mean)
# 0.92 0.98 0.99 1.00 1.00
apply(conv.16.10$conv.ai,1,mean)
# 1   1   1   1   1 

apply(conv.16.10$rel.conv.chol,1,mean)
# 0.97 0.98 0.99 1.00 0.99 
apply(conv.16.10$conv.chol,1,mean)
# 1   1   1   1   1 

conv.12.8=conv.analy(pi.core.fit.12.8.4.12,pi.core.fit.12.8.4.24,pi.core.fit.12.8.4.48,
                      pi.core.fit.12.8.4.96,pi.core.fit.12.8.4.192,c("12","24","48","96","192"))

apply(conv.12.8$rel.conv.ai,1,mean)
#0.81 0.98 0.99 0.99 1.00 
apply(conv.12.8$conv.ai,1,mean)
# 1   1   1   1   1 

apply(conv.12.8$rel.conv.chol,1,mean)
#0.86 0.98 0.98 1.00 1.00 
apply(conv.12.8$conv.chol,1,mean)
# 1   1   1   1   1 

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
