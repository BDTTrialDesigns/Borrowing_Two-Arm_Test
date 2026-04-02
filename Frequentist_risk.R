#This script computes the power function and frequentist risk shown in the supplementary material

library(RBesT)
library(ggplot2)
library(latex2exp)

rm(list=ls())

setwd("~/")
source("Functions_Subm_N.R")


# Parameters --------------------------------------------------------------

sigma=1 #variance of one observation
nt=20      #sample size in treatment arm
nc=20    #sample size in control arm
n0c=10    #external data ESS
sigmapc=sigma/sqrt(n0c) #SD control arm prior
muc=0                 #prior mean control

cf=4*sqrt(sigma^2/nc+sigma^2/n0c)  #point of full discard
cp=4                  #steepness of change

alpha.up=0.075        #max TIE for CDs
alpha.low=0.01        #min TIE for CDs
alpha.b=0.025         #FD significance level
th0=0

priorPars=rbind(c(muc,muc),c(sigmapc,Inf))

w.mix=0.3
bound=FALSE

nmc=2*10^5
seed=1234

#grid for control arm parameter
thc=round(seq(muc-0.4,muc+0.6,0.2),1)
#grid of effect sizes
delta=c(seq(-0.6,-0.05,0.05),-0.035,seq(-0.01,0.01,0.001),0.035,seq(0.05,1.2,0.05))
eval.grid=expand.grid(thc,delta)


# Probability to reject ---------------------------------------------------

pRej=sapply(1:nrow(eval.grid), function(x)
{
  set.seed(seed)
  yc=rnorm(nmc,eval.grid[x,1],sigma/sqrt(nc))
  yt=rnorm(nmc,eval.grid[x,2]+eval.grid[x,1],sigma/sqrt(nt))
  
  #Weight vector, 1 for monte carlo
  wmatT=rep(1/nmc,nmc)
  
  priorPars.v=priorPars
  priorPars.v[2,]=c(Inf,Inf)
  p.h1.unif=pH1(yc=yc,yt=yt,nc=nc,nt=nt,
                priorPars=priorPars.v,w.mix=w.mix,
                outcome='normal',th0=th0)
  
  p.h1.info=pH1(yc=yc,yt=yt,nc=nc,nt=nt,
                priorPars=priorPars,w.mix=w.mix,
                outcome='normal',th0=th0)
  
  p.h1.mix=pH1(prior='mix',yc=yc,yt=yt,nc=nc,nt=nt,
               priorPars=priorPars,w.mix=w.mix,
               outcome='normal',th0=th0)
  
  p.h1.EB=pH1(prior='EB',yc=yc,yt=yt,nc=nc,nt=nt,
              priorPars=priorPars,w.mix=w.mix,
              outcome='normal',bound=bound,th0=th0)
  
  kappa.f=cThreshold(yc=yc,yt=yt,nc=nc,nt=nt,
                     priorPars=priorPars,
                     outcome='normal', mc.results=mc.results, 
                     alpha.b=alpha.b, alpha.up=alpha.up, 
                     alpha.low=alpha.low, th0=th0,
                     cf=cf,cp=cp,CD='Discard')
  
  
  kappa.f.c=cThreshold(yc=yc,yt=yt,nc=nc,nt=nt,
                       priorPars=priorPars,
                       outcome='normal', mc.results=mc.results, 
                       alpha.b=alpha.b, alpha.up=alpha.up, 
                       alpha.low=alpha.low, th0=th0,
                       cf=cf,cp=cp,CD='Constrain')
  
  #test decision
  dec.list = cbind(p.h1.info>(1-alpha.b),
                   p.h1.unif>(1-alpha.b),
                   p.h1.mix>(1-alpha.b),
                   p.h1.EB>(1-alpha.b),
                   p.h1.unif>(1-kappa.f),
                   p.h1.unif>(1-kappa.f.c))
  
  
  powf= wmatT %*% dec.list
  powf
  
})

decisions=c("BD","FD","RMD-Unit", "EBPowD","CD-Discard","CD-Constraint")

pRejD=data.frame('pRej'=as.vector(pRej), 'thc'=rep(eval.grid[,1],each=nrow(pRej)),
                 'delta'= rep(eval.grid[,2],each=nrow(pRej)),'Decision'= rep(decisions,nrow(eval.grid)))
pRejD$FRisk=ifelse(pRejD$delta<=0,0.975*pRejD$pRej,0.025*(1-pRejD$pRej))

pRejD$Decision <- factor(pRejD$Decision, levels = c("CD-Constraint", "CD-Discard", 
                                                  "RMD-Unit","EBPowD","BD","FD"))

save.image("~/results/freq_risk.RData")

# Plots -------------------------------------------------------------------

load("~/results/freq_risk.RData")

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


theme=ggplot()+theme_light()

powerfn <- theme+ geom_line(aes(delta, pRej,  color=Decision),
                       size=0.5, data=subset(pRejD, pRejD$thc<=0.8)) + ylab("Probability to reject") +
  theme(legend.title=element_blank(), legend.position = "top") +
  scale_colour_manual(values=colors) + guides(colour = guide_legend(nrow = 1)) +
  xlab(TeX("$\\delta$"))+ facet_wrap('thc', scales='free_y', nrow=2)

frisk <- theme+ geom_line(aes(delta, FRisk, group=Decision, color=Decision),
                       size=0.5, data=subset(pRejD, pRejD$thc<=0.8 & pRejD$delta<=0)) + 
  geom_line(aes(delta, FRisk, group=Decision, color=Decision),
            size=0.5, data=subset(pRejD, pRejD$thc<=0.8 & pRejD$delta>0)) + 
  ylab("Frequentist risk") +
  theme(legend.title=element_blank(), legend.position = "top") +
  scale_colour_manual(values=colors) + guides(colour = guide_legend(nrow = 1)) +
  geom_point(aes(x=delta, y=FRisk, color=Decision), 
             data=subset(pRejD, pRejD$thc<=0.8 & pRejD$delta==0)) +
  xlab(TeX("$\\delta$")) + facet_wrap('thc', scales='free_y', nrow=2)


powerfnzoom <- theme+ geom_line(aes(delta, pRej,  color=Decision),
                       size=0.5, data=subset(pRejD, pRejD$thc==0 & pRejD$delta>=0.2 & pRejD$delta<=0.3)) + ylab("Probability to reject") +
   geom_point(aes(delta, pRej,  color=Decision),
              size=0.5, data=subset(pRejD, pRejD$thc==0 & pRejD$delta>=0.2 & pRejD$delta<=0.3)) + ylab("Probability to reject") +
  theme(legend.title=element_blank(), legend.position = "top") +
  scale_colour_manual(values=colors) + guides(colour = guide_legend(nrow = 1)) +
  xlab(TeX("$\\delta$"))


cairo_pdf(file ="~/results/power_function_freq.pdf", width=7, height=5, pointsize=12, onefile=TRUE)
print(powerfn)
dev.off()

cairo_pdf(file ="~/results/freq_risk.pdf", width=7, height=5, pointsize=12, onefile=TRUE)
print(frisk)
dev.off()

cairo_pdf(file ="~/results/freq_risk_EBcase.pdf", width=7, height=5, pointsize=12, onefile=TRUE)
print(powerfnzoom)
dev.off()
