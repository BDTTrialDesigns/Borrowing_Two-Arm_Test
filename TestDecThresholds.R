#This script produces evaluations of the test decision thresholds (conditional on data outcomes)
#for the simulated data

library(RBesT)
library(ggplot2)
library(latex2exp)
require(RColorBrewer)
#library(GISTools)
library(patchwork)

rm(list=ls())

setwd("~/")
source("Functions_Subm_N.R")


# Parameters --------------------------------------------------------------

sigma=1         #variance of one observation
nt=20           #sample size in treatment arm
nc=20           #sample size in control arm
sigmapc=sigma/sqrt(10) #SD control arm prior
muc=0                   #prior mean control

cf=4*sqrt(sigma^2/nc+sigmapc^2) #point of full discard
cp=4                  #steepness of change

alpha.up=0.075        #max TIE for CDs
alpha.low=0.01        #min TIE for CDs
alpha.b=0.025         #frequentist TIE
th0=0
priorPars=rbind(c(muc,muc),c(sigmapc,Inf))
Ac=sigma^2/(nc*sigmapc^2)

w.mix=0.3
bound=FALSE


# Compute rejection region ------------------------------------------------

#grid
thc=seq(muc-2,muc+2,0.01)
delta=seq(-2,2,0.001)
eval.grid=expand.grid(thc,delta)

yc=eval.grid[,1]
dobs=eval.grid[,2]
yt=eval.grid[,1]+eval.grid[,2]

#z statistics
z=dobs/sqrt((sigma^2/nc)+ ((sigma^2/nt)))

#post probabilties of H1
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


p.h1.list=cbind(p.h1.info,p.h1.unif,p.h1.mix,p.h1.EB,p.h1.unif,p.h1.unif)

#Test decisions
dec.list = cbind(p.h1.info>(1-alpha.b),
                 p.h1.unif>(1-alpha.b),
                 p.h1.mix>(1-alpha.b),
                 p.h1.EB>(1-alpha.b),
                 p.h1.unif>(1-kappa.f),
                 p.h1.unif>(1-kappa.f.c))


dat.dec=data.frame(zval=z,yc=yc,p.h1=as.vector(p.h1.list),
                    dec=as.vector(dec.list),
                    Decision=rep(c('BD','FD','RMD-Unit','EBPowD',
                                   'CD-Discard','CD-Constraint'),each=length(yc)))

#Critical values
zstat=tapply(dat.dec$zval*dat.dec$dec,list(dat.dec$Decision,as.factor(dat.dec$yc)),
             function(x) x[abs(x)>0][which.min(x[abs(x)>0])])

# Test decision thresholds for the CDs - sensitivity --------------------------------

grid=seq(-2,2,0.001)

#threshold to obtain the BD under the frequentist test
kappa.pi=kappa.cal.f(yc=grid, alpha.b=alpha.b, At=0, Ac=Ac, mupc=muc, 
                     mupt=muc, sigma=sigma, nt=nt, nc=nc, th0=th0) 
#threshold to obtain the BD under the Bayes test (check)
kappa.bayes.pi=kappa.bayes.general(alpha.k=kappa.pi, eval=grid, Ac=Ac, 
                                   mupc=muc, sigma=sigma, nt=nt, nc=nc, th0=th0) 

#threshold for the CD-Constraint under the frequentist test
kappa.CDC=pmax(pmin(kappa.pi,alpha.up),alpha.low)
#threshold for the CD-Constraint under the Bayes test
kappa.bayes.CDC=kappa.bayes.general(alpha.k=kappa.CDC, eval=grid, Ac=Ac, 
                                    mupc=muc, sigma=sigma, nt=nt, nc=nc, th0=th0) 

#threshold to obtain the FD under the Bayes test
kappa.bayes.FD=kappa.bayes.general(alpha.k=alpha.b, eval=grid, Ac=Ac, 
                                   mupc=muc, sigma=sigma, nt=nt, nc=nc, th0=th0) 

##thresholds for the CD-Discard under different p

cp=0.5
#frequentist test
kappa.CDD.05=pmax(pmin(kappa.discount(alpha.k=alpha.b, alpha.k2=kappa.pi, 
                                      eval=grid, muc=muc, cf=cf, cp=cp),alpha.up),alpha.low)
#Bayes test
kappa.bayes.CDD.05=kappa.bayes.general(alpha.k=kappa.CDD.05, eval=grid, Ac=Ac, 
                                       mupc=muc, sigma=sigma, nt=nt, nc=nc,  th0=th0) 

cp=2
#frequentist test
kappa.CDD.2=pmax(pmin(kappa.discount(alpha.k=alpha.b, alpha.k2=kappa.pi, eval=grid, 
                                     muc=muc, cf=cf, cp=cp),alpha.up),alpha.low)
#Bayes test
kappa.bayes.CDD.2=kappa.bayes.general(alpha.k=kappa.CDD.2, eval=grid, Ac=Ac, 
                                      mupc=muc, sigma=sigma, nt=nt, nc=nc,  th0=th0) 

cp=4
#frequentist test
kappa.CDD.4=pmax(pmin(kappa.discount(alpha.k=alpha.b, alpha.k2=kappa.pi, 
                                     eval=grid, muc=muc, cf=cf, cp=cp),alpha.up),alpha.low)
#Bayes test
kappa.bayes.CDD.4=kappa.bayes.general(alpha.k=kappa.CDD.4, eval=grid, Ac=Ac, 
                                      mupc=muc, sigma=sigma, nt=nt, nc=nc,  th0=th0) 


# Datasets for plotting ---------------------------------------------------

#critical values
dat.plot.dec=data.frame(yc=as.numeric(rep(colnames(zstat),nrow(zstat))),
                        crit=as.vector(t(zstat)),
                        kappa=as.vector(1-pnorm(t(zstat))),
                        kappa.bayes=
                          as.vector(kappa.bayes.general(1-pnorm(t(zstat)),
                                                        eval=thc, Ac=Ac, mupc=muc, 
                                                        sigma=sigma, nt=nt, nc=nc,  
                                                        th0=th0)),
                        dstat=as.vector(t(zstat)*sqrt((sigma^2/nc)+((sigma^2/nt)))),
                        Decision=rep(rownames(zstat),each=ncol(zstat)))

#posterior probabilities
dat.plot.probs=data.frame(ph1=p.h1.list[,'p.h1.info'],
                          ph1f=p.h1.list[,'p.h1.unif'],
                          yc=yc,delta=eval.grid[,2])


#thresholds for frquentist test (kappa)
boundPlot=data.frame(grid=rep(grid,6),kappa=c(kappa.pi,kappa.CDD.05,kappa.CDD.2,
                                              kappa.CDD.4,kappa.CDC, rep(alpha.b,length(grid))),
                     Decision=rep(c("BD","CD-Discard (0.5)","CD-Discard (2)","CD-Discard (4)",
                                    "CD-Constraint","FD"),each=length(grid)))

#thresholds for Bayes test (gamma)
boundPlot.bayes=data.frame(grid=rep(grid,6),kappa=c(kappa.bayes.pi,kappa.bayes.CDD.05,kappa.bayes.CDD.2,kappa.bayes.CDD.4,
                                                    kappa.bayes.CDC, kappa.bayes.FD),
                           Decision=rep(c("BD","CD-Discard (0.5)","CD-Discard (2)","CD-Discard (4)",
                                          "CD-Constraint","FD"),each=length(grid)))


#reorder factors
dat.plot.dec$Decision <- factor(dat.plot.dec$Decision, levels = c("CD-Constraint","CD-Discard","CD-Discard (4)",
                                                                  "CD-Discard (2)","CD-Discard (0.5)",
                                                                  "RMD-Unit","EBPowD","BD","FD"))

boundPlot$Decision <- factor(boundPlot$Decision, levels = c("CD-Constraint","CD-Discard", "CD-Discard (4)",
                                                            "CD-Discard (2)","CD-Discard (0.5)",
                                                            "RMD-Unit","EBPowD","BD","FD"))

boundPlot.bayes$Decision <- factor(boundPlot.bayes$Decision, levels = c("CD-Constraint","CD-Discard","CD-Discard (4)",
                                                                        "CD-Discard (2)","CD-Discard (0.5)",
                                                                        "RMD-Unit","EBPowD","BD","FD"))


save.image("~/results/TestDecThresholds.RData")

# Plots -------------------------------------------------------------------

load("~/results/TestDecThresholds.RData")

colors = c(
  "FD" = "#D55E00",
  "BD" = "#0072B2",
  "RMD-Unit" = "#56B4E9",
  "EBPowD" = "#009E73",
  "CD-Constraint" = "#CC79A7",
  "CD-Discard" = "#7B3294",
  "CD-Discard (1)" = "#1B9E77",
  "CD-Discard (2)" = "#A6761D",
  "CD-Discard (0.5)" = "#E7298A",
  "CD-Discard (4)" = "#7B3294"
)


#Figure 1

kappa.bayes.plot.0= ggplot() + theme_light() +
  geom_line(data = subset(boundPlot.bayes, boundPlot.bayes$Decision=='FD' | boundPlot.bayes$Decision=='BD'),
            aes(x = grid, y = kappa, color=Decision), size=1) +
  xlab(TeX("$\\bar{y}_C-\\mu_{C}$"))+
  ylab(TeX("$\\gamma(\\bar{y}_C)$"))+
  ggtitle(TeX("Bayes test: reject $H_0$ if $ $ $\ \\frac{\\mu_{\\delta|y}-\\delta_0}{\\sigma_{\\delta|y}} \\geq z_{1-\\gamma(\\bar{y}_C)}$")) +
  ylim(c(0,0.4))+ 
  scale_colour_manual(values=colors)+
  theme(legend.title=element_blank(),legend.position="top")+ guides(colour = guide_legend(nrow = 1))+
  theme(plot.title = element_text(size=12)) 


kappa.freq.plot.0= ggplot() + theme_light() +
  geom_line(data = subset(boundPlot, boundPlot$Decision=='FD' | boundPlot$Decision=='BD'),
            aes(x = grid, y = kappa, color=Decision), size=1) +
  xlab(TeX("$\\bar{y}_C-\\mu_{C}$"))+
  ylab(TeX("$\\kappa(\\bar{y}_C)$"))+
  ggtitle(TeX("Frequentist test: reject $H_0$ if $ $ $\ \\frac{\\bar{y}_T-\\bar{y}_C-\\delta_0}{\\sqrt{\\sigma^2/n_T + \\sigma^2/n_C}} \\geq z_{1-\\kappa(\\bar{y}_C)}$")) +
  ylim(c(0,0.4))+ 
  scale_colour_manual(values=colors)+
  theme(legend.title=element_blank(),legend.position="top")+ guides(colour = guide_legend(nrow = 1))+
  theme(plot.title = element_text(size=12)) 

figure1=
  (kappa.freq.plot.0 + guides(col = guide_legend(nrow = 1, byrow = T))) | 
  kappa.bayes.plot.0 + guides(color = "none") +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

pdf(file ="~/results/figure1.pdf", width=9, height=6.5, pointsize=12, onefile=TRUE)
print(figure1)
dev.off()

#Figure 2

kappa.plot=ggplot() + theme_light() +
  geom_hline(yintercept = alpha.up, linetype='dashed', color='grey30') +
  geom_hline(yintercept = alpha.low, linetype='dashed', color='grey30') +
  geom_line(data = dat.plot.dec,
            aes(x=yc, y=kappa, color=Decision), size=1) +
  xlab(TeX("$\\bar{y}_C-\\mu_{C}$"))+
  ylab(TeX("$\\kappa(\\bar{y}_C)$"))+
  ylim(c(0,0.2))+ 
  scale_colour_manual(values=colors)+
  geom_hline(yintercept = alpha.b, linetype='dashed', colour="#D55E00")+ 
  annotate("text", x = min(dat.plot.dec$yc)-0.6, y = 0.005, label = TeX("$\\alpha^{LOW}$"), 
           hjust = 0, vjust = -0.5, size = 3.5) +
  annotate("text", x = min(dat.plot.dec$yc)-0.6, y = 0.07, label = TeX("$\\alpha^{UP}$"), 
           hjust = 0, vjust = -0.5, size = 3.5) +
  theme(legend.title=element_blank(),legend.position="top")


crit.plot = ggplot() +theme_light() +
  geom_hline(yintercept = qnorm(1-alpha.up), linetype='dashed', color='grey30') +
  geom_hline(yintercept = qnorm(1-alpha.low), linetype='dashed', color='grey30') +
  geom_line(data=subset(dat.plot.dec, dat.plot.dec$Decision!='CD-Discard'), 
                        aes(x=yc, y=crit, color=Decision),linewidth=0.8) + 
  geom_line(data=subset(dat.plot.dec, dat.plot.dec$Decision=='CD-Discard'), 
            aes(x=yc, y=crit, color=Decision),linewidth=0.8) + 
  scale_colour_manual(values=colors) +
  ylab(TeX("$z_{1-\\kappa(\\bar{y}_{C})}$")) + ylim(c(0.5,3))+
  xlab(TeX("$\\bar{y}_{C}-\\mu_C$")) +
  theme(legend.title=element_blank(),legend.position="top")
  
gamma.plot= ggplot() + theme_light() +
  geom_line(data = dat.plot.dec,
            aes(x = yc, y = kappa.bayes, color=Decision), size=1) +
  xlab(TeX("$\\bar{y}_C-\\mu_{C}$"))+
  ylab(TeX("$\\gamma(\\bar{y}_C)$"))+
  ylim(c(0,0.4))+ 
  scale_colour_manual(values=colors)+
  theme(legend.title=element_blank(),legend.position="top")


figure2 = (
  ((crit.plot + guides(col = guide_legend(nrow = 1, byrow = TRUE))) /
      (kappa.plot + guides(color = "none"))) |
    (gamma.plot + guides(color = "none"))) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

pdf(file ="~/results/figure2.pdf", width=9, height=6.5, pointsize=12, onefile=TRUE)
print(figure2)
dev.off()

#Sensitivity (Supplementary Figure 5 figure)

gamma.plot.sensitivity= ggplot() + theme_light() +
  geom_line(data = boundPlot.bayes,
            aes(x = grid, y = kappa, color=Decision), size=1) +
  xlab(TeX("$\\bar{y}_C-\\mu_{C}$"))+
  ylab(TeX("$\\gamma(\\bar{y}_C)$"))+
  ylim(c(0,0.4))+ 
  scale_colour_manual(values=colors)+
theme(legend.title=element_blank(),legend.position="top")+ guides(colour = guide_legend(nrow = 1))+
  theme(plot.title = element_text(size=18)) 


kappa.plot.sensitivity= ggplot() + theme_light() +
  geom_line(data = boundPlot,
            aes(x = grid, y = kappa, color=Decision), size=1) +
  xlab(TeX("$\\bar{y}_C-\\mu_{C}$"))+
  ylab(TeX("$\\kappa(\\bar{y}_C)$"))+
  ylim(c(0,0.4))+ 
  scale_colour_manual(values=colors)+
  theme(legend.title=element_blank(),legend.position="top")+ guides(colour = guide_legend(nrow = 1))+
  theme(plot.title = element_text(size=18)) 


sens=
  ((kappa.plot.sensitivity + guides(col = guide_legend(nrow = 1, byrow = T))) | 
  gamma.plot.sensitivity + guides(color = "none")) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

pdf(file ="~/results/sens.pdf", width=9, height=6.5, pointsize=12, onefile=TRUE)
print(sens)
dev.off()


#Posterior probability and critical values (supplementary figure 3)
ph1.crit <- ggplot() +
  geom_contour_filled(aes(x = yc, y = delta, z =(1-ph1)),
                      data=subset(dat.plot.probs, dat.plot.probs$d>0 & 
                                    dat.plot.probs$d<=1 & dat.plot.probs$yc<=2),
                      breaks = (1-c(0,0.8,0.9,0.95,0.96,0.97,0.98,0.99,0.995,1)), alpha=0.5) + 
  scale_fill_viridis_d(name=TeX("$P^{\\pi}(H_0|y)$"))+
  theme_minimal() + 
  geom_line(aes(x = yc, y = dstat, color = Decision), size=1.2, 
            data=subset(dat.plot.dec,dat.plot.dec$dstat>0 & 
                          dat.plot.dec$dstat<=1 & dat.plot.dec$yc<=2))+ 
  scale_colour_manual(values=colors)+
  xlab(TeX("$\\bar{y}_{C}-\\mu_C$"))+ 
  ylab(TeX("$\\bar{y}_{T}-\\bar{y}_{C}")) 

cairo_pdf(file ="~/results/ph1_thr.pdf", width=8, height=5, pointsize=12, onefile=TRUE)
print(ph1.crit)
dev.off()

#Posterior probability difference and critical values (supplementary figure 4)
ph1Diff <-  ggplot() +
  geom_contour_filled(aes(x = yc, y = delta, z = ph1-ph1f),
                      subset(dat.plot.probs,dat.plot.probs$d>0 & dat.plot.probs$d<=1 & 
                               ((dat.plot.probs$ph1>0.975 & dat.plot.probs$ph1f<=0.975) | 
                               (dat.plot.probs$ph1<=0.975 & dat.plot.probs$ph1f>0.975))),
                      breaks = c(-1,-0.5,-0.2,-0.1,-0.05,0,0.05,0.1,0.2,0.5,1), alpha=0.5) + 
  scale_fill_viridis_d(name=TeX("$P^{\\pi_0}(H_0|y)-P^{\\pi}(H_0|y)$")) +
  theme_minimal()+ 
  geom_line(aes(x = yc, y = dstat, color = Decision), size=1.2, 
            data=subset(dat.plot.dec,dat.plot.dec$dstat>0 & dat.plot.dec$dstat<=1))+ 
  scale_colour_manual(values=colors)+
  xlab(TeX("$\\bar{y}_{C}-\\mu_C$"))+ 
  ylab(TeX("$\\bar{y}_{T}-\\bar{y}_{C}"))


cairo_pdf(file ="~/results/ph1_diff.pdf", width=8, height=5, pointsize=12, onefile=TRUE)
print(ph1Diff)
dev.off()
