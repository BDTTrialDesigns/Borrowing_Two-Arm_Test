#This script generates the frequentist OCs for the different simulation setups.

rm(list=ls())

library(parallel)
library(ggplot2)
library(patchwork)
library(latex2exp)
library(pwr)

setwd("~/")
source("Functions_Subm_N.R")

# alpha.up=0.1
# alpha.low=0.001
alpha.up=0.075
alpha.low=0.01
#ct=4
ct=2
cp=4
outcome='binomial'
#outcome='normal'

###normal outcomes
if (outcome=='normal')
{
mupc=0                        #control arm informative prior mean
w.mix.v=c(0.25,0.3,0.5,0.75)  #mixture prior weights
conflict=seq(-2,2,0.02)
sig=1                         #data standard deviation of one observation

setups=lapply(w.mix.v, function(w.mix) 
{
setups.l=data.frame(rbind(
c(nc=20,nt=20,n0c=20,mupc=mupc,w.mix=w.mix,alpha.up=alpha.up,alpha.low=alpha.low,outcome=outcome),
c(nc=100,nt=100,n0c=100,mupc=mupc,w.mix=w.mix,alpha.up=alpha.up,alpha.low=alpha.low,outcome=outcome),
c(nc=20,nt=20,n0c=10,mupc=mupc,w.mix=w.mix,alpha.up=alpha.up,alpha.low=alpha.low,outcome=outcome),
c(nc=20,nt=20,n0c=40,mupc=mupc,w.mix=w.mix,alpha.up=alpha.up,alpha.low=alpha.low,outcome=outcome),
c(nc=100,nt=100,n0c=50,mupc=mupc,w.mix=w.mix,alpha.up=alpha.up,alpha.low=alpha.low,outcome=outcome),
c(nc=100,nt=100,n0c=200,mupc=mupc,w.mix=w.mix,alpha.up=alpha.up,alpha.low=alpha.low,outcome=outcome),
c(nc=50,nt=100,n0c=50,mupc=mupc,w.mix=w.mix,alpha.up=alpha.up,alpha.low=alpha.low,outcome=outcome),
c(nc=20,nt=100,n0c=80,mupc=mupc,w.mix=w.mix,alpha.up=alpha.up,alpha.low=alpha.low,outcome=outcome),
c(nc=10,nt=20,n0c=10,mupc=mupc,w.mix=w.mix,alpha.up=alpha.up,alpha.low=alpha.low,outcome=outcome),
c(nc=4,nt=20,n0c=16,mupc=mupc,w.mix=w.mix,alpha.up=alpha.up,alpha.low=alpha.low,outcome=outcome)
))
setups.l[,1:7] <- lapply(setups.l[,1:7], as.numeric)
r=setups.l$nc/setups.l$nt
setups.l$delta=sqrt(((qnorm(1-0.025)+qnorm(0.8))^2*((1/r)+1)*sig^2)/setups.l$nt) #effect for power evaluations                                                                      #number of Monte Carlo simulations

setups.l
})

setups=do.call(rbind,setups)
nmc=2*10^5 
}

###binary outcomes
if (outcome=='binomial')
{
mupc=0.35 #control arm informative prior mean
w.mix=0.3 #mixture prior weight
conflict=seq(-mupc,1-mupc,0.01) #conflict evaluation range

setups=data.frame(rbind(
c(nc=20,nt=20,n0c=20,mupc=mupc,w.mix=w.mix,alpha.up=alpha.up,alpha.low=alpha.low,outcome=outcome),
c(nc=100,nt=100,n0c=100,mupc=mupc,w.mix=w.mix,alpha.up=alpha.up,alpha.low=alpha.low,outcome=outcome),
c(nc=20,nt=20,n0c=10,mupc=mupc,w.mix=w.mix,alpha.up=alpha.up,alpha.low=alpha.low,outcome=outcome),
c(nc=20,nt=20,n0c=40,mupc=mupc,w.mix=w.mix,alpha.up=alpha.up,alpha.low=alpha.low,outcome=outcome),
c(nc=100,nt=100,n0c=50,mupc=mupc,w.mix=w.mix,alpha.up=alpha.up,alpha.low=alpha.low,outcome=outcome),
c(nc=100,nt=100,n0c=200,mupc=mupc,w.mix=w.mix,alpha.up=alpha.up,alpha.low=alpha.low,outcome=outcome),
c(nc=50,nt=100,n0c=50,mupc=mupc,w.mix=w.mix,alpha.up=alpha.up,alpha.low=alpha.low,outcome=outcome),
c(nc=20,nt=100,n0c=80,mupc=mupc,w.mix=w.mix,alpha.up=alpha.up,alpha.low=alpha.low,outcome=outcome),
c(nc=10,nt=20,n0c=10,mupc=mupc,w.mix=w.mix,alpha.up=alpha.up,alpha.low=alpha.low,outcome=outcome),
c(nc=4,nt=20,n0c=16,mupc=mupc,w.mix=w.mix,alpha.up=alpha.up,alpha.low=alpha.low,outcome=outcome)
))

setups[,1:7] <- lapply(setups[,1:7], as.numeric)
setups$delta=sapply(1:nrow(setups), function(s) #treatment effect for power evaluation
{
  h=pwr.2p2n.test(h=NULL, n1=setups$nc[s], n2=setups$nt[s], sig.level=0.025, power=0.8, alternative='greater')$h
  p2=(sin((2*asin(sqrt(setups$mupc[s]))+h)/2))^2
  return(p2-setups$mupc[s])
}
  )
}

out=mclapply(1:nrow(setups), function(z)
{
  if(setups$outcome[z]=='normal')
  {
    priorPars=rbind(c(setups$mupc[z],setups$mupc[z]),c(sig/sqrt(setups$n0c[z]),Inf)) #prior parameters (control and treated arm)
    cp=cp #steepness of CD-Discard
    cf=ct*sqrt(sig^2/setups$nc[z]+sig^2/setups$n0c[z]) #point of full discard of CD-Discard
    
  }
  
  if(setups$outcome[z]=='binomial')
  {
    priorPars=cbind(c(setups$mupc[z]*setups$n0c[z],(1-setups$mupc[z])*setups$n0c[z]),c(0,0)) #prior parameters (control and treated arm)
    cp=cp #steepness of CD-Discard
    cf=ct #SDs for full discard of CD-Discard
  }
  
  ocs=compute.oc(outcome=setups$outcome[z], 
                 nt=setups$nt[z], #treatment arm sample size
                 nc=setups$nc[z], #control arm sample size
                 n0c=setups$n0c[z], #historical info effective sample size
                 delta=setups$delta[z], #treatment effect for power evaluation
                 conflict=conflict,
                 priorPars=priorPars,
                 sigma=sig, #sd data (normal outcomes only)
                 alpha.low=setups$alpha.low[z], #min TIE for CDs
                 alpha.up=setups$alpha.up[z], #max TIE for CDs
                 alpha.b = 0.025, #significance level of FD
                 w.mix=setups$w.mix[z], #mixture prior weight
                 cp=cp,
                 cf=cf,
                 bound=FALSE, #bound EB informativeness? if TRUE, PP is discounted to current sample size as maximum
                 nmc=nmc, #number of MC samples (only for outcome='normal')
                 seed=1234,
                 ncores=28)
  
  
  out=data.frame('Freq.type.I'=as.vector(ocs$TIE),
                 'Freq.Power'=as.vector(ocs$Pow),
                 'Type.I.cal'=as.vector(ocs$TIEcal),
                 'Power.cal'=as.vector(ocs$Powcal),
                 'Conflict'=rep(as.numeric(rownames(ocs$TIE)),ncol(ocs$TIE)),
                 'delta'=rep(setups$delta[z],ncol(ocs$TIE)*nrow(ocs$TIE)),
                 'nc'=rep(setups$nc[z],ncol(ocs$TIE)*nrow(ocs$TIE)),
                 'nt'=rep(setups$nt[z],ncol(ocs$TIE)*nrow(ocs$TIE)),
                 'alpha.up'=rep(setups$alpha.up[z],ncol(ocs$TIE)*nrow(ocs$TIE)),
                 'alpha.low'=rep(setups$alpha.low[z],ncol(ocs$TIE)*nrow(ocs$TIE)),
                 'mupc'=rep(setups$mupc[z],ncol(ocs$TIE)*nrow(ocs$TIE)),
                 'n0c'=rep(setups$n0c[z],ncol(ocs$TIE)*nrow(ocs$TIE)),
                 'Decision'=rep(colnames(ocs$TIE),each=nrow(ocs$TIE)))
  
  colors = c(
    "FD" = "#D55E00",
    "BD" = "#0072B2",
    "RMD-Unit" = "#56B4E9",
    "EBPowD" = "#009E73",
    "CD-Constraint" = "#CC79A7",
    "CD-Discard" = "#7B3294",
    "CD-Discard (1)" = "#1B9E77",
    "CD-Discard (2)" = "#A6761D",
    "CD-Discard (3)" = "#E7298A",
    "CD-Discard (4)" = "#7B3294"
  )
  
  theme=ggplot()+theme_light() +
    scale_colour_manual(values=colors) +
    xlab(expression(theta[C]-mu[C]))

  out$Decision <- factor(out$Decision, levels = c("CD-Constraint", "CD-Discard",
                                                  "RMD-Unit","EBPowD","BD","FD"))
  
  if(outcome=='normal')
  {
    g1 <- theme +
      geom_segment(aes(x=-2, xend=2,y = out$alpha.up[1],yend = out$alpha.up[1]), col='red', linetype='dashed', size=0.3)+
      geom_segment(aes(x=-2, xend=2,y = out$alpha.low[1],yend = out$alpha.low[1]), col='red', linetype='dashed', size=0.3)+ 
      geom_line(aes(Conflict, Freq.type.I,  color=Decision),size=0.5, data=out) + ylab("TIE rate") + 
      ylim(c(0,0.2)) + guides(fill = "none")
    
    g2 <- theme+ geom_line(aes(Conflict, Freq.Power,  color=Decision),
                           size=0.5, data=out) + ylab(TeX('Power')) + ylim(c(0,1)) + theme(legend.position = "none")
    g4 <- theme+ geom_line(aes(Conflict, Power.cal,  color=Decision),
                           size=0.5, data=out) + ylab(TeX('Power')) + theme(legend.position = "none")
  }
  
  
  if(outcome=='binomial')
  {
    
    g1 <- theme +
      geom_segment(aes(x=-out$mupc[1], xend=(1-out$mupc[1]),y = out$alpha.up[1],yend = out$alpha.up[1]), col='red', linetype='dashed', size=0.3)+
      geom_segment(aes(x=-out$mupc[1], xend=(1-out$mupc[1]),y = out$alpha.low[1],yend = out$alpha.low[1]), col='red', linetype='dashed', size=0.3)+ 
      geom_line(aes(Conflict, Freq.type.I,  color=Decision),size=0.5, data=out) + ylab("TIE rate") + 
      ylim(c(0,0.2)) + guides(fill = "none")
    
    g2 <- theme+ geom_line(aes(Conflict, Freq.Power,  color=Decision),
                           size=0.5, data=out) + ylab(TeX('Power')) + ylim(c(0,1))  + theme(legend.position = "none")
    
    g4 <- theme+ geom_line(aes(Conflict, Power.cal,  color=Decision),
                           size=0.5, data=out) + ylab(TeX('Power')) + theme(legend.position = "none")
  }
  
  g3 <- theme+ guides(fill = "none") + geom_line(aes(Conflict, Type.I.cal,  color=Decision),
                         size=0.5, data=out) + ylab("TIE rate") + ylim(c(0,0.2)) 
  
  
  
  a=(g1 + guides(col = guide_legend(nrow = 1, byrow = T))+theme(legend.title=element_blank()) | 
       g2+ guides(color = "none")) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top")
  
  b=(g3 + guides(col = guide_legend(nrow = 1, byrow = T))+theme(legend.title=element_blank()) | 
       g4+ guides(color = "none")) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top")
  
  
  pdf(file = paste0("~/results/",setups$outcome[z],"_n0c:",setups$n0c[z],
                    "_nc:",setups$nc[z],
                    "_nt:",setups$nt[z],
                    "_alpha.up:",setups$alpha.up[z],
                    "_alpha.low:",setups$alpha.low[z],
                    "_cp:",cp,
                    "_ct:",ct,
                    "_w.mix:",setups$w.mix[z],"_OC.pdf"), width=8, height=4, pointsize=12, onefile=TRUE)
  print(a)
  dev.off()
  
  pdf(file = paste0("~/results/",setups$outcome[z],"_n0c:",setups$n0c[z],
                    "_nc:",setups$nc[z],
                    "_nt:",setups$nt[z],
                    "_alpha.up:",setups$alpha.up[z],
                    "_alpha.low:",setups$alpha.low[z],
                    "_cp:",cp,
                    "_ct:",ct,
                    "_w.mix:",setups$w.mix[z],
                    "_OC_calibrated.pdf"), width=8, height=4, pointsize=12, onefile=TRUE)
  print(b)
  dev.off()
  
  save(out, file=paste0("~/results/",setups$outcome[z],"_n0c:",setups$n0c[z],
                        "_nc:",setups$nc[z],
                        "_nt:",setups$nt[z],
                        "_alpha.up:",setups$alpha.up[z],
                        "_alpha.low:",setups$alpha.low[z],
                        "_cp:",cp,
                        "_ct:",ct,
                        "_w.mix:",setups$w.mix[z],"_OCdat.RData",sep=""))
  out
},mc.cores=1)

# ################################################################################################################
# ####Tables
# 
# library(xtable)
# 
# load("normal_n0c:50_n:100_alpha.up:0.075_cf:0.692820323027551_cp:4_w.mix:0.3_OCdat.RData")
# 
# out_NoConflict=subset(out,out$Conflict==0)
# tabOCs=out_NoConflict[,c('Freq.type.I','Freq.Power')]*100
# rownames(tabOCs)=out_NoConflict$Decision
# maxTIEcal=tapply(out[,c('Type.I.cal')], out$Decision, function(x) max(x)*100)[rownames(tabOCs)]
# maxPowcal=tapply(out[,c('Power.cal')], out$Decision, function(x) max(x)*100)[rownames(tabOCs)]
# tabOCs100=cbind(tabOCs,maxTIEcal-maxTIEcal['FD'],maxPowcal-maxPowcal['FD'])
# 
# load("normal_n0c:10_n:20_alpha.up:0.075_cf:1.54919333848297_cp:4_w.mix:0.3_OCdat.RData")
# 
# out_NoConflict20=subset(out,out$Conflict==0)
# tabOCs=out_NoConflict20[,c('Freq.type.I','Freq.Power')]*100
# rownames(tabOCs)=out_NoConflict$Decision
# maxTIEcal=tapply(out[,c('Type.I.cal')], out$Decision, function(x) max(x)*100)[rownames(tabOCs)]
# maxPowcal=tapply(out[,c('Power.cal')], out$Decision, function(x) max(x)*100)[rownames(tabOCs)]
# tabOCs20=cbind(tabOCs,maxTIEcal-maxTIEcal['FD'],maxPowcal-maxPowcal['FD'])
# 
# xtable(cbind(tabOCs20, tabOCs100))
# 
# 
