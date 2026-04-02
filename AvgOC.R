#This script computes average OCs and power function

library(RBesT)
library(ggplot2)
library(latex2exp)
library(patchwork)

rm(list=ls())

setwd("~/")
source("Functions_Subm_N.R")

#####parameters

outcome='normal'
sigma=1     #variance of one observation
mupc=0      #prior mean control

alpha.up=0.075         #max TIE 
alpha.low=0.01         #min TIE
alpha.b=0.025          #frequentist TIE
th0=0                  #null hypothesis border
nmc=2*10^5             #number of Monte Carlo simulations
seed=1234              #random seed

delta=seq(-0.3,1,0.025) #effect values grid

setups=data.frame(rbind(
  c(nc=20,nt=20,n0c=20,mupc=mupc,alpha.up=alpha.up,alpha.low=alpha.low,outcome=outcome),
  c(nc=100,nt=100,n0c=100,mupc=mupc,alpha.up=alpha.up,alpha.low=alpha.low,outcome=outcome),
  c(nc=20,nt=20,n0c=10,mupc=mupc,alpha.up=alpha.up,alpha.low=alpha.low,outcome=outcome),
  c(nc=20,nt=20,n0c=40,mupc=mupc,alpha.up=alpha.up,alpha.low=alpha.low,outcome=outcome),
  c(nc=100,nt=100,n0c=50,mupc=mupc,alpha.up=alpha.up,alpha.low=alpha.low,outcome=outcome),
  c(nc=100,nt=100,n0c=200,mupc=mupc,alpha.up=alpha.up,alpha.low=alpha.low,outcome=outcome),
  c(nc=50,nt=100,n0c=50,mupc=mupc,alpha.up=alpha.up,alpha.low=alpha.low,outcome=outcome),
  c(nc=20,nt=100,n0c=80,mupc=mupc,alpha.up=alpha.up,alpha.low=alpha.low,outcome=outcome),
  c(nc=10,nt=20,n0c=10,mupc=mupc,alpha.up=alpha.up,alpha.low=alpha.low,outcome=outcome),
  c(nc=4,nt=20,n0c=16,mupc=mupc,alpha.up=alpha.up,alpha.low=alpha.low,outcome=outcome)
))
setups[,1:6] <- lapply(setups[,1:6], as.numeric)

dat.AvgOC=mclapply(1:nrow(setups), function(z)
{
  #borrowing parameters
  cp=4  #steepness
  ct=4
  cf=ct*sqrt(1/setups$nc[z]+1/setups$n0c[z]) #point of full discard
  
#obtain alpha.up and alpha.low for the calibrated comparisons, i.e. Robust mixture and EB with same 
#max and min TIE rate as the compromise

load(paste0("~/results/",setups$outcome[z],"_n0c:",setups$n0c[z],
                        "_nc:",setups$nc[z],
                        "_nt:",setups$nt[z],
                        "_alpha.up:",setups$alpha.up[z],
                        "_alpha.low:",setups$alpha.low[z],
                        "_cp:",cp,
                        "_ct:",ct,
                        "_w.mix:",0.3,"_OCdat.RData",sep=""))
  
alpha.low.mix.07=min(out$Freq.type.I[out$Decision=='RMD-Unit'])     #min TIE
alpha.up.mix.07=max(out$Freq.type.I[out$Decision=='RMD-Unit'])         #max TIE

load(paste0("~/results/",setups$outcome[z],"_n0c:",setups$n0c[z],
            "_nc:",setups$nc[z],
            "_nt:",setups$nt[z],
            "_alpha.up:",setups$alpha.up[z],
            "_alpha.low:",setups$alpha.low[z],
            "_cp:",cp,
            "_ct:",ct,
            "_w.mix:",0.25,"_OCdat.RData",sep=""))
alpha.low.mix.075=min(out$Freq.type.I[out$Decision=='RMD-Unit'])     #min TIE
alpha.up.mix.075=max(out$Freq.type.I[out$Decision=='RMD-Unit'])         #max TIE

load(paste0("~/results/",setups$outcome[z],"_n0c:",setups$n0c[z],
            "_nc:",setups$nc[z],
            "_nt:",setups$nt[z],
            "_alpha.up:",setups$alpha.up[z],
            "_alpha.low:",setups$alpha.low[z],
            "_cp:",cp,
            "_ct:",ct,
            "_w.mix:",0.5,"_OCdat.RData",sep=""))
alpha.low.mix.05=min(out$Freq.type.I[out$Decision=='RMD-Unit'])     #min TIE
alpha.up.mix.05=max(out$Freq.type.I[out$Decision=='RMD-Unit'])         #max TIE

load(paste0("~/results/",setups$outcome[z],"_n0c:",setups$n0c[z],
            "_nc:",setups$nc[z],
            "_nt:",setups$nt[z],
            "_alpha.up:",setups$alpha.up[z],
            "_alpha.low:",setups$alpha.low[z],
            "_cp:",cp,
            "_ct:",ct,
            "_w.mix:",0.75,"_OCdat.RData",sep=""))
alpha.low.mix.025=min(out$Freq.type.I[out$Decision=='RMD-Unit'])     #min TIE
alpha.up.mix.025=max(out$Freq.type.I[out$Decision=='RMD-Unit'])         #max TIE
alpha.low.eb=min(out$Freq.type.I[out$Decision=='EBPowD'])
alpha.up.eb=max(out$Freq.type.I[out$Decision=='EBPowD'])


#prior parameters
priorPars=rbind(c(setups$mupc[z],setups$mupc[z]),c(1/sqrt(setups$n0c[z]),Inf))


Avg.OC=sapply(delta, function(x)
{
  set.seed(seed)
  
  #sample control param.
  thc=rnorm(nmc,priorPars[1,1],priorPars[2,1]) 
  
  #sample data
  yc=rnorm(nmc,thc,sigma/sqrt(setups$nc[z]))
  yt=rnorm(nmc,thc+x,sigma/sqrt(setups$nt[z]))
  
  #Weight vector, 1 for monte carlo
  wmatT=rep(1/nmc,nmc)
  
  priorPars.v=priorPars
  priorPars.v[2,]=c(Inf,Inf)
  p.h1.unif=pH1(yc=yc,yt=yt,nc=setups$nc[z],nt=setups$nt[z],
                priorPars=priorPars.v,
                outcome=setups$outcome[z],th0=th0)
  
  p.h1.info=pH1(yc=yc,yt=yt,nc=setups$nc[z],nt=setups$nt[z],
                priorPars=priorPars,
                outcome=setups$outcome[z],th0=th0)
  
  p.h1.mix=pH1(prior='mix',yc=yc,yt=yt,nc=setups$nc[z],nt=setups$nt[z],
               priorPars=priorPars,w.mix=0.3,
               outcome=setups$outcome[z],th0=th0)
  
  p.h1.mix.075=pH1(prior='mix',yc=yc,yt=yt,nc=setups$nc[z],nt=setups$nt[z],
                   priorPars=priorPars,w.mix=0.25,
                   outcome=setups$outcome[z],th0=th0)
  
  p.h1.mix.05=pH1(prior='mix',yc=yc,yt=yt,nc=setups$nc[z],nt=setups$nt[z],
                  priorPars=priorPars,w.mix=0.5,
                  outcome=setups$outcome[z],th0=th0)
  
  p.h1.mix.025=pH1(prior='mix',yc=yc,yt=yt,nc=setups$nc[z],nt=setups$nt[z],
                   priorPars=priorPars,w.mix=0.75,
                   outcome=setups$outcome[z],th0=th0)
  
  p.h1.EB=pH1(prior='EB',yc=yc,yt=yt,nc=setups$nc[z],nt=setups$nt[z],
              priorPars=priorPars,w.mix=w.mix,
              outcome=setups$outcome[z],th0=th0)
  
  
  kappa.f=cThreshold(yc=yc,yt=yt,nc=setups$nc[z],nt=setups$nt[z],
                     priorPars=priorPars,
                     outcome=setups$outcome[z], mc.results=mc.results, 
                     alpha.b=alpha.b, alpha.up=alpha.up, 
                     alpha.low=alpha.low, th0=th0,
                     cf=cf,cp=cp,CD='Discard')
  
  
  kappa.f.c=cThreshold(yc=yc,yt=yt,nc=setups$nc[z],nt=setups$nt[z],
                       priorPars=priorPars,
                       outcome=setups$outcome[z], mc.results=mc.results, 
                       alpha.b=alpha.b, alpha.up=alpha.up, 
                       alpha.low=alpha.low, th0=th0,
                       cf=cf,cp=cp,CD='Constrain')
  
  kappa.f.mix.075=cThreshold(yc=yc,yt=yt,nc=setups$nc[z],nt=setups$nt[z],
                             priorPars=priorPars,
                             outcome=setups$outcome[z], mc.results=mc.results,
                             alpha.b=alpha.b, alpha.up=alpha.up.mix.075, alpha.low=alpha.low.mix.075, th0=th0,
                             cf=cf,cp=cp,CD='Discard')
  
  kappa.f.c.mix.075=cThreshold(yc=yc,yt=yt,nc=setups$nc[z],nt=setups$nt[z],
                               priorPars=priorPars,
                               outcome=setups$outcome[z], mc.results=mc.results,
                               alpha.b=alpha.b, alpha.up=alpha.up.mix.075, alpha.low=alpha.low.mix.075, th0=th0,
                               cf=cf,cp=cp,CD='Constrain')
  
  kappa.f.mix.05=cThreshold(yc=yc,yt=yt,nc=setups$nc[z],nt=setups$nt[z],
                            priorPars=priorPars,
                            outcome=setups$outcome[z], mc.results=mc.results,
                            alpha.b=alpha.b, alpha.up=alpha.up.mix.05, alpha.low=alpha.low.mix.05, th0=th0,
                            cf=cf,cp=cp,CD='Discard')
  
  kappa.f.c.mix.05=cThreshold(yc=yc,yt=yt,nc=setups$nc[z],nt=setups$nt[z],
                              priorPars=priorPars,
                              outcome=setups$outcome[z], mc.results=mc.results,
                              alpha.b=alpha.b, alpha.up=alpha.up.mix.05, alpha.low=alpha.low.mix.05, th0=th0,
                              cf=cf,cp=cp,CD='Constrain')
  
  kappa.f.mix.025=cThreshold(yc=yc,yt=yt,nc=setups$nc[z],nt=setups$nt[z],
                             priorPars=priorPars,
                             outcome=setups$outcome[z], mc.results=mc.results,
                             alpha.b=alpha.b, alpha.up=alpha.up.mix.025, alpha.low=alpha.low.mix.025, th0=th0,
                             cf=cf,cp=cp,CD='Discard')
  
  kappa.f.c.mix.025=cThreshold(yc=yc,yt=yt,nc=setups$nc[z],nt=setups$nt[z],
                               priorPars=priorPars,
                               outcome=setups$outcome[z], mc.results=mc.results,
                               alpha.b=alpha.b, alpha.up=alpha.up.mix.025, alpha.low=alpha.low.mix.025, th0=th0,
                               cf=cf,cp=cp,CD='Constrain')
  
  kappa.f.eb=cThreshold(yc=yc,yt=yt,nc=setups$nc[z],nt=setups$nt[z],
                        priorPars=priorPars,
                        outcome=setups$outcome[z], mc.results=mc.results,
                        alpha.b=alpha.b, alpha.up=alpha.up.eb, alpha.low=alpha.low.eb, th0=th0,
                        cf=cf,cp=cp,CD='Discard')
  
  kappa.f.c.eb=cThreshold(yc=yc,yt=yt,nc=setups$nc[z],nt=setups$nt[z],
                          priorPars=priorPars,
                          outcome=setups$outcome[z], mc.results=mc.results,
                          alpha.b=alpha.b, alpha.up=alpha.up.eb, alpha.low=alpha.low.eb, th0=th0,
                          cf=cf,cp=cp,CD='Constrain')
  
  #test decision
  p.h1.list=cbind(p.h1.info,p.h1.unif,p.h1.mix,
                  p.h1.EB)
  
  dec.list = cbind(p.h1.info>(1-alpha.b),
                   p.h1.unif>(1-alpha.b),
                   p.h1.mix>(1-alpha.b),
                   p.h1.EB>(1-alpha.b),
                   p.h1.unif>(1-kappa.f),
                   p.h1.unif>(1-kappa.f.c),
                   p.h1.mix.075>(1-alpha.b),
                   p.h1.mix.05>(1-alpha.b),
                   p.h1.mix.025>(1-alpha.b),
                   p.h1.unif>(1-kappa.f.mix.075),
                   p.h1.unif>(1-kappa.f.c.mix.075),
                   p.h1.unif>(1-kappa.f.mix.05),
                   p.h1.unif>(1-kappa.f.c.mix.05),
                   p.h1.unif>(1-kappa.f.mix.025),
                   p.h1.unif>(1-kappa.f.c.mix.025),
                   p.h1.unif>(1-kappa.f.eb),
                   p.h1.unif>(1-kappa.f.c.eb))
  
  
  AvgOC= wmatT %*% dec.list
  AvgOC
})


Decision=c("BD","FD","RMD-Unit", "EBPowD","CD-Discard","CD-Constraint",
           "RMD-Unit.075","RMD-Unit.05","RMD-Unit.025",
           "CD-Discard.075","CD-Constraint.075",
           "CD-Discard.05","CD-Constraint.05","CD-Discard.025","CD-Constraint.025",
           "CD-Discard.EB","CD-Constraint.EB")

data.frame(PrRej=as.vector(t(Avg.OC)),delta=rep(delta,nrow(Avg.OC)),
                     Decision=rep(Decision,each=length(delta)))
},mc.cores=nrow(setups))

save.image("~/results/AvgOCs.RData")
###

load("~/results/AvgOCs.RData")

j=1

AvgOC=dat.AvgOC[[j]]

dat.AvgOC.plot=subset(AvgOC, AvgOC$Decision=="BD" |
                        AvgOC$Decision=="FD"|
                        AvgOC$Decision=="RMD-Unit"|
                        AvgOC$Decision=="EBPowD"|
                        AvgOC$Decision=="CD-Discard"|
                        AvgOC$Decision=="CD-Constraint")


dat.AvgOC.plot$Decision <- factor(dat.AvgOC.plot$Decision, levels = c( "CD-Constraint", "CD-Discard",
                                                                      "RMD-Unit","EBPowD","BD","FD"))

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


avgTIE <-  
  ggplot(data=subset(dat.AvgOC.plot,dat.AvgOC.plot$delta>-0.25 & dat.AvgOC.plot$delta<=0)) +
  theme_light() +
  aes(x = delta, y = PrRej, colour = Decision) +
  scale_colour_manual(values=colors) +
  xlab(TeX("$\\delta$"))+ ylab("Avg. TIE rate")+
  geom_line(linewidth = 0.5) +
  geom_hline(yintercept = 0.025, linetype='dashed')+
  theme(legend.title=element_blank(), legend.position = "top") 

avgPow1 <-  
  ggplot(data=subset(dat.AvgOC.plot,dat.AvgOC.plot$delta<0.15 & dat.AvgOC.plot$delta>=0)) +
  theme_light() +
  aes(x = delta, y = PrRej, colour = Decision) +
  scale_colour_manual(values=colors) +
  xlab(TeX("$\\delta$"))+ ylab("Avg. Power")+
  geom_line(linewidth = 0.5) +
  geom_hline(yintercept = 0.025, linetype='dashed')+
  theme(legend.title=element_blank(), legend.position = "top") 

avgPow2 <-  
  ggplot(data=subset(dat.AvgOC.plot,dat.AvgOC.plot$delta<=0.35 & dat.AvgOC.plot$delta>=0.15)) +
  theme_light() +
  aes(x = delta, y = PrRej, colour = Decision) +
  scale_colour_manual(values=colors) +
  xlab(TeX("$\\delta$"))+ ylab("Avg. Power")+
  geom_line(linewidth = 0.5) +
  theme(legend.title=element_blank(), legend.position = "top") 


avgPow3 <-  
  ggplot(data=subset(dat.AvgOC.plot,dat.AvgOC.plot$delta<1 & dat.AvgOC.plot$delta>=0.35)) +
  theme_light() +
  aes(x = delta, y = PrRej, colour = Decision) +
  scale_colour_manual(values=colors) +
  xlab(TeX("$\\delta$"))+ ylab("Avg. Power")+
  geom_line(linewidth = 0.5) +
  theme(legend.title=element_blank(), legend.position = "top") 


avg.power=((avgTIE+ guides(col = guide_legend(nrow = 1, byrow = T))+theme(legend.title=element_blank()) | 
     avgPow1 + guides(color = "none")) /
     (avgPow2 + guides(color = "none") | avgPow3+ guides(color = "none"))) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")


pdf(file = paste0("~/results/AvgOC_",j,".pdf"), width=8, height=6, pointsize=10, onefile=TRUE)
print(avg.power)
dev.off()

####Table

#half external control sample size compared current

library(xtable)

load("~/results/AvgOCs.RData")

j=3
AvgOC=dat.AvgOC[[j]]

AvgOCtab=subset(AvgOC, (round(AvgOC$delta,2)==0 | round(AvgOC$delta,2)==0.25 | 
                                  round(AvgOC$delta,2)==0.50 |
                                  round(AvgOC$delta,2)==0.75 | round(AvgOC$delta,2)==1))

AvgOC.tab=data.frame(matrix(AvgOCtab$PrRej,17,5, byrow = TRUE))
rownames(AvgOC.tab)=AvgOC$Decision[round(AvgOC$delta,2)==0 ]
colnames(AvgOC.tab)=c('0','0.25','0.5','0.75','1')
AvgOC.tab.20=AvgOC.tab

j=5
AvgOC=dat.AvgOC[[j]]

AvgOCtab=subset(AvgOC, (round(AvgOC$delta,2)==0 | round(AvgOC$delta,2)==0.10 | 
                                  round(AvgOC$delta,2)==0.2 |
                                  round(AvgOC$delta,2)==0.3 | round(AvgOC$delta,2)==0.4))

AvgOC.tab=data.frame(matrix(AvgOCtab$PrRej,17,5, byrow = TRUE))
rownames(AvgOC.tab)=AvgOC$Decision[round(AvgOC$delta,2)==0 ]
colnames(AvgOC.tab)=c('0','0.1','0.2','0.3','0.4')
AvgOC.tab=cbind(AvgOC.tab.20,AvgOC.tab)


xtable(AvgOC.tab*100)

#double external control sample size compared to current

j=4
AvgOC=dat.AvgOC[[j]]

AvgOCtab=subset(AvgOC, (round(AvgOC$delta,2)==0 | round(AvgOC$delta,2)==0.25 | 
                          round(AvgOC$delta,2)==0.50 |
                          round(AvgOC$delta,2)==0.75 | round(AvgOC$delta,2)==1))

AvgOC.tab=data.frame(matrix(AvgOCtab$PrRej,17,5, byrow = TRUE))
rownames(AvgOC.tab)=AvgOC$Decision[round(AvgOC$delta,2)==0 ]
colnames(AvgOC.tab)=c('0','0.25','0.5','0.75','1')
AvgOC.tab.20=AvgOC.tab

j=6
AvgOC=dat.AvgOC[[j]]

AvgOCtab=subset(AvgOC, (round(AvgOC$delta,2)==0 | round(AvgOC$delta,2)==0.10 | 
                          round(AvgOC$delta,2)==0.2 |
                          round(AvgOC$delta,2)==0.3 | round(AvgOC$delta,2)==0.4))

AvgOC.tab=data.frame(matrix(AvgOCtab$PrRej,17,5, byrow = TRUE))
rownames(AvgOC.tab)=AvgOC$Decision[round(AvgOC$delta,2)==0 ]
colnames(AvgOC.tab)=c('0','0.1','0.2','0.3','0.4')
AvgOC.tab=cbind(AvgOC.tab.20,AvgOC.tab)


xtable(AvgOC.tab*100)


#external control sample size same as current (sum equal to treatment arm sample size)

j=9
AvgOC=dat.AvgOC[[j]]

AvgOCtab=subset(AvgOC, (round(AvgOC$delta,2)==0 | round(AvgOC$delta,2)==0.25 | 
                          round(AvgOC$delta,2)==0.50 |
                          round(AvgOC$delta,2)==0.75 | round(AvgOC$delta,2)==1))

AvgOC.tab=data.frame(matrix(AvgOCtab$PrRej,17,5, byrow = TRUE))
rownames(AvgOC.tab)=AvgOC$Decision[round(AvgOC$delta,2)==0 ]
colnames(AvgOC.tab)=c('0','0.25','0.5','0.75','1')
AvgOC.tab.20=AvgOC.tab

j=7
AvgOC=dat.AvgOC[[j]]

AvgOCtab=subset(AvgOC, (round(AvgOC$delta,2)==0 | round(AvgOC$delta,2)==0.10 | 
                          round(AvgOC$delta,2)==0.2 |
                          round(AvgOC$delta,2)==0.3 | round(AvgOC$delta,2)==0.4))

AvgOC.tab=data.frame(matrix(AvgOCtab$PrRej,17,5, byrow = TRUE))
rownames(AvgOC.tab)=AvgOC$Decision[round(AvgOC$delta,2)==0 ]
colnames(AvgOC.tab)=c('0','0.1','0.2','0.3','0.4')
AvgOC.tab=cbind(AvgOC.tab.20,AvgOC.tab)


xtable(AvgOC.tab*100)