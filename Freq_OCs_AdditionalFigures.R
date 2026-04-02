#This script produces plots of the frequentist OCs (TIE rate, power)
#for the simulation study

rm(list=ls())

library(ggplot2)
library(patchwork)
library(latex2exp)

dir="~/results/"
setwd("~/results/")

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

alpha.up=0.075
alpha.low=0.01
ct=4
cp=4

# Normal outcomes ---------------------------------------------------------

outcome='normal'
mupc=0                        #control arm informative prior mean
w.mix.v=c(0.3)  #mixture prior weights
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


# Figure 3: Normal outcomes simulation study - main text ------------------

z=3

load(file=paste0(dir,setups$outcome[z],"_n0c:",setups$n0c[z],
                      "_nc:",setups$nc[z],
                      "_nt:",setups$nt[z],
                      "_alpha.up:",setups$alpha.up[z],
                      "_alpha.low:",setups$alpha.low[z],
                      "_cp:",cp,
                      "_ct:",ct,
                      "_w.mix:",setups$w.mix[z],"_OCdat.RData",sep=""))



out$Decision <- factor(out$Decision, levels = c("CD-Constraint", "CD-Discard",
                                                "RMD-Unit","EBPowD","BD","FD"))

#TIE rate & Power (uncalibrated)

t1e <- theme +
  geom_segment(aes(x=-2, xend=2,y = out$alpha.up[1],yend = out$alpha.up[1]), col='red', linetype='dashed', size=0.3)+
  geom_segment(aes(x=-2, xend=2,y = out$alpha.low[1],yend = out$alpha.low[1]), col='red', linetype='dashed', size=0.3)+ 
    geom_line(aes(Conflict, Freq.type.I,  color=Decision),size=0.5, data=out) + 
  ylab("TIE rate") +  ylim(c(0,0.2)) + guides(fill = "none")+
  theme(legend.title=element_blank())

power <- theme+ 
  geom_line(aes(Conflict, Freq.Power,  color=Decision),size=0.5, data=out) + 
  ylab(TeX('Power')) + ylim(c(0,1)) + 
  theme(legend.position = "none")
  
#TIE rate & Power (calibrated)
  
t1ecal <- theme+ 
  guides(fill = "none") + 
  geom_line(aes(Conflict, Type.I.cal,  color=Decision),
            size=0.5, data=out) + 
  ylab("TIE rate (recalibrated)") + ylim(c(0,0.2))+
  theme(legend.title=element_blank())
  
powercal <- theme+ 
  geom_line(aes(Conflict, Power.cal,  color=Decision),size=0.5, data=out) + 
  ylab(TeX('Power (recalibrated)')) + 
  theme(legend.position = "none")


figure3 = (
  ((t1e + guides(col = guide_legend(nrow = 1, byrow = TRUE))) /
     (power + guides(color = "none"))) |
    ((t1ecal + guides(col = guide_legend(nrow = 1, byrow = TRUE))) /
       (powercal + guides(color = "none")))) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

pdf(file ="figure3.pdf", width=9, height=6.5, pointsize=12, onefile=TRUE)
print(figure3)
dev.off()


# Supplementary Figures ---------------------------------------------------

#Same setup as main text but large sample size

z=5

load(file=paste0(dir,setups$outcome[z],"_n0c:",setups$n0c[z],
                 "_nc:",setups$nc[z],
                 "_nt:",setups$nt[z],
                 "_alpha.up:",setups$alpha.up[z],
                 "_alpha.low:",setups$alpha.low[z],
                 "_cp:",cp,
                 "_ct:",ct,
                 "_w.mix:",setups$w.mix[z],"_OCdat.RData",sep=""))


out$Decision <- factor(out$Decision, levels = c("CD-Constraint", "CD-Discard",
                                                "RMD-Unit","EBPowD","BD","FD"))

#TIE rate & Power (uncalibrated)

t1e <- theme +
  geom_segment(aes(x=-2, xend=2,y = out$alpha.up[1],yend = out$alpha.up[1]), col='red', linetype='dashed', size=0.3)+
  geom_segment(aes(x=-2, xend=2,y = out$alpha.low[1],yend = out$alpha.low[1]), col='red', linetype='dashed', size=0.3)+ 
  geom_line(aes(Conflict, Freq.type.I,  color=Decision),size=0.5, data=out) + 
  ylab("TIE rate") +  ylim(c(0,0.2)) + guides(fill = "none")+
  theme(legend.title=element_blank())

power <- theme+ 
  geom_line(aes(Conflict, Freq.Power,  color=Decision),size=0.5, data=out) + 
  ylab(TeX('Power')) + ylim(c(0,1)) + 
  theme(legend.position = "none")

#TIE rate & Power (calibrated)

t1ecal <- theme+ 
  guides(fill = "none") + 
  geom_line(aes(Conflict, Type.I.cal,  color=Decision),
            size=0.5, data=out) + 
  ylab("TIE rate (recalibrated)") + ylim(c(0,0.2))+
  theme(legend.title=element_blank())

powercal <- theme+ 
  geom_line(aes(Conflict, Power.cal,  color=Decision),size=0.5, data=out) + 
  ylab(TeX('Power (recalibrated)')) + 
  theme(legend.position = "none")


large_ss_example = (
  ((t1e + guides(col = guide_legend(nrow = 1, byrow = TRUE))) /
     (power + guides(color = "none"))) |
    ((t1ecal + guides(col = guide_legend(nrow = 1, byrow = TRUE))) /
       (powercal + guides(color = "none")))) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

pdf(file ="large_ss_example.pdf", width=9, height=6.5, pointsize=12, onefile=TRUE)
print(large_ss_example)
dev.off()


#Additional small sample size scenarios (1,4,9,10)

out.plot=NULL
t=1

for (z in c(1,4,9,10))
{

load(file=paste0(dir,setups$outcome[z],"_n0c:",setups$n0c[z],
                 "_nc:",setups$nc[z],
                 "_nt:",setups$nt[z],
                 "_alpha.up:",setups$alpha.up[z],
                 "_alpha.low:",setups$alpha.low[z],
                 "_cp:",cp,
                 "_ct:",ct,
                 "_w.mix:",setups$w.mix[z],"_OCdat.RData",sep=""))

  

  out$Decision <- factor(out$Decision, levels = c("CD-Constraint", "CD-Discard",
                                                  "RMD-Unit","EBPowD","BD","FD"))
  
  #TIE rate & Power (uncalibrated)
  
  t1e <- theme +
    geom_segment(aes(x=-2, xend=2,y = out$alpha.up[1],yend = out$alpha.up[1]), col='red', linetype='dashed', size=0.3)+
    geom_segment(aes(x=-2, xend=2,y = out$alpha.low[1],yend = out$alpha.low[1]), col='red', linetype='dashed', size=0.3)+ 
    geom_line(aes(Conflict, Freq.type.I,  color=Decision),size=0.5, data=out) + 
    ylab("TIE rate") +  ylim(c(0,0.3)) + guides(fill = "none")+
    theme(legend.title=element_blank()) + 
    ggtitle(TeX(sprintf("$n_{0C} = %s,\\; n_{C} = %s,\\; n_{T} = %s$",
                setups$n0c[z],
                setups$nc[z],
                setups$nt[z])))+
  theme( plot.title = element_text(size = 8))
  
  power <- theme+ 
    geom_line(aes(Conflict, Freq.Power,  color=Decision),size=0.5, data=out) + 
    ylab(TeX('Power')) + ylim(c(0.4,1)) + 
    theme(legend.position = "none") + 
    ggtitle(TeX(sprintf("$n_{0C} = %s,\\; n_{C} = %s,\\; n_{T} = %s,\\; \\delta = %s$",
                        setups$n0c[z],
                        setups$nc[z],
                        setups$nt[z],
                        round(out$delta[1],2)))) +
    theme(plot.title = element_text(size = 8)) 
  
  out.plot[[t]]= t1e
  out.plot[[t+1]]= power
  t=t+2
}

small_ss_normal = (
  ((out.plot[[1]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
     (out.plot[[2]] + guides(color = "none"))) /
    ((out.plot[[3]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[4]] + guides(color = "none"))) / 
    ((out.plot[[5]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[6]] + guides(color = "none"))) /
    ((out.plot[[7]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[8]] + guides(color = "none"))))  +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

pdf(file ="small_ss_normal.pdf", width=9, height=12, pointsize=12, onefile=TRUE)
print(small_ss_normal)
dev.off()

#Additional large sample size scenarios (2,6,7,8)

out.plot=NULL
t=1

for (z in c(2,6,7,8))
{
  
  load(file=paste0(dir,setups$outcome[z],"_n0c:",setups$n0c[z],
                   "_nc:",setups$nc[z],
                   "_nt:",setups$nt[z],
                   "_alpha.up:",setups$alpha.up[z],
                   "_alpha.low:",setups$alpha.low[z],
                   "_cp:",cp,
                   "_ct:",ct,
                   "_w.mix:",setups$w.mix[z],"_OCdat.RData",sep=""))
  
  
  
  out$Decision <- factor(out$Decision, levels = c("CD-Constraint", "CD-Discard",
                                                  "RMD-Unit","EBPowD","BD","FD"))
  
  #TIE rate & Power (uncalibrated)
  
  t1e <- theme +
    geom_segment(aes(x=-2, xend=2,y = out$alpha.up[1],yend = out$alpha.up[1]), col='red', linetype='dashed', size=0.3)+
    geom_segment(aes(x=-2, xend=2,y = out$alpha.low[1],yend = out$alpha.low[1]), col='red', linetype='dashed', size=0.3)+ 
    geom_line(aes(Conflict, Freq.type.I,  color=Decision),size=0.5, data=out) + 
    ylab("TIE rate") +  ylim(c(0,0.3)) + guides(fill = "none")+
    theme(legend.title=element_blank()) + 
    ggtitle(TeX(sprintf("$n_{0C} = %s,\\; n_{C} = %s,\\; n_{T} = %s$",
                        setups$n0c[z],
                        setups$nc[z],
                        setups$nt[z])))+
    theme( plot.title = element_text(size = 8))
  
  power <- theme+ 
    geom_line(aes(Conflict, Freq.Power,  color=Decision),size=0.5, data=out) + 
    ylab(TeX('Power')) + ylim(c(0.4,1)) + 
    theme(legend.position = "none") + 
    ggtitle(TeX(sprintf("$n_{0C} = %s,\\; n_{C} = %s,\\; n_{T} = %s,\\; \\delta = %s$",
                        setups$n0c[z],
                        setups$nc[z],
                        setups$nt[z],
                        round(out$delta[1],2)))) +
    theme(plot.title = element_text(size = 8)) 
  
  out.plot[[t]]= t1e
  out.plot[[t+1]]= power
  t=t+2
}

large_ss_normal = (
  ((out.plot[[1]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
     (out.plot[[2]] + guides(color = "none"))) /
    ((out.plot[[3]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[4]] + guides(color = "none"))) / 
    ((out.plot[[5]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[6]] + guides(color = "none"))) /
    ((out.plot[[7]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[8]] + guides(color = "none"))))  +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

pdf(file ="large_ss_normal.pdf", width=9, height=12, pointsize=12, onefile=TRUE)
print(large_ss_normal)
dev.off()


# Binomial outcomes -------------------------------------------------------

outcome='binomial'
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

#Small sample size scenarios (1,3,4,9,10)

out.plot=NULL
t=1

for (z in c(1,3,4,9,10))
{
  
  load(file=paste0(dir,setups$outcome[z],"_n0c:",setups$n0c[z],
                   "_nc:",setups$nc[z],
                   "_nt:",setups$nt[z],
                   "_alpha.up:",setups$alpha.up[z],
                   "_alpha.low:",setups$alpha.low[z],
                   "_cp:",cp,
                   "_ct:",ct,
                   "_w.mix:",setups$w.mix[z],"_OCdat.RData",sep=""))
  
  
  
  out$Decision <- factor(out$Decision, levels = c("CD-Constraint", "CD-Discard",
                                                  "RMD-Unit","EBPowD","BD","FD"))
  
  #TIE rate & Power (uncalibrated)
  
  t1e <- theme +
    geom_segment(aes(x=-out$mupc[1], xend=(1-out$mupc[1]),y = out$alpha.up[1],yend = out$alpha.up[1]), col='red', linetype='dashed', size=0.3)+
    geom_segment(aes(x=-out$mupc[1], xend=(1-out$mupc[1]),y = out$alpha.low[1],yend = out$alpha.low[1]), col='red', linetype='dashed', size=0.3)+ 
    geom_line(aes(Conflict, Freq.type.I,  color=Decision),size=0.5, data=out) + 
    ylab("TIE rate") +  ylim(c(0,0.35)) + guides(fill = "none")+
    theme(legend.title=element_blank()) + 
    ggtitle(TeX(sprintf("$n_{0C} = %s,\\; n_{C} = %s,\\; n_{T} = %s$",
                        setups$n0c[z],
                        setups$nc[z],
                        setups$nt[z])))+
    theme( plot.title = element_text(size = 8))
  
  power <- theme+ 
    geom_line(aes(Conflict, Freq.Power,  color=Decision),size=0.5, data=out) + 
    ylab(TeX('Power')) + ylim(c(0.4,1)) + xlim(c(-out$mupc[1],1-out$mupc[1]-out$delta[1])) + 
    theme(legend.position = "none") + 
    ggtitle(TeX(sprintf("$n_{0C} = %s,\\; n_{C} = %s,\\; n_{T} = %s,\\; \\delta = %s$",
                        setups$n0c[z],
                        setups$nc[z],
                        setups$nt[z],
                        round(out$delta[1],2)))) +
    theme(plot.title = element_text(size = 8)) 
  
  out.plot[[t]]= t1e
  out.plot[[t+1]]= power
  t=t+2
}

small_ss_binomial = (
  ((out.plot[[1]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
     (out.plot[[2]] + guides(color = "none"))) /
    ((out.plot[[3]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[4]] + guides(color = "none"))) / 
    ((out.plot[[5]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[6]] + guides(color = "none"))) /
    ((out.plot[[7]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
    (out.plot[[8]] + guides(color = "none"))) /
    ((out.plot[[9]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[10]] + guides(color = "none"))))  +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

pdf(file ="small_ss_binomial.pdf", width=9, height=12, pointsize=12, onefile=TRUE)
print(small_ss_binomial)
dev.off()

#Additional large sample size scenarios (2,5,6,7,8)

out.plot=NULL
t=1

for (z in c(2,5,6,7,8))
{
  
  load(file=paste0(dir,setups$outcome[z],"_n0c:",setups$n0c[z],
                   "_nc:",setups$nc[z],
                   "_nt:",setups$nt[z],
                   "_alpha.up:",setups$alpha.up[z],
                   "_alpha.low:",setups$alpha.low[z],
                   "_cp:",cp,
                   "_ct:",ct,
                   "_w.mix:",setups$w.mix[z],"_OCdat.RData",sep=""))
  
  
  
  out$Decision <- factor(out$Decision, levels = c("CD-Constraint", "CD-Discard",
                                                  "RMD-Unit","EBPowD","BD","FD"))
  
  #TIE rate & Power (uncalibrated)
  
  t1e <- theme +
    geom_segment(aes(x=-out$mupc[1], xend=(1-out$mupc[1]),y = out$alpha.up[1],yend = out$alpha.up[1]), col='red', linetype='dashed', size=0.3)+
    geom_segment(aes(x=-out$mupc[1], xend=(1-out$mupc[1]),y = out$alpha.low[1],yend = out$alpha.low[1]), col='red', linetype='dashed', size=0.3)+ 
    geom_line(aes(Conflict, Freq.type.I,  color=Decision),size=0.5, data=out) + 
    ylab("TIE rate") +  ylim(c(0,0.35)) + guides(fill = "none")+
    theme(legend.title=element_blank()) + 
    ggtitle(TeX(sprintf("$n_{0C} = %s,\\; n_{C} = %s,\\; n_{T} = %s$",
                        setups$n0c[z],
                        setups$nc[z],
                        setups$nt[z])))+
    theme( plot.title = element_text(size = 8))
  
  power <- theme+ 
    geom_line(aes(Conflict, Freq.Power,  color=Decision),size=0.5, data=out) + 
    ylab(TeX('Power')) + ylim(c(0.4,1))  + xlim(c(-out$mupc[1],1-out$mupc[1]-out$delta[1])) + 
    theme(legend.position = "none") + 
    ggtitle(TeX(sprintf("$n_{0C} = %s,\\; n_{C} = %s,\\; n_{T} = %s,\\; \\delta = %s$",
                        setups$n0c[z],
                        setups$nc[z],
                        setups$nt[z],
                        round(out$delta[1],2)))) +
    theme(plot.title = element_text(size = 8)) 
  
  out.plot[[t]]= t1e
  out.plot[[t+1]]= power
  t=t+2
}

large_ss_binomial = (
    ((out.plot[[1]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
     (out.plot[[2]] + guides(color = "none"))) /
    ((out.plot[[3]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[4]] + guides(color = "none"))) / 
    ((out.plot[[5]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[6]] + guides(color = "none"))) /
    ((out.plot[[7]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
    (out.plot[[8]] + guides(color = "none"))) /
  ((out.plot[[9]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
     (out.plot[[10]] + guides(color = "none"))))  +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

pdf(file ="large_ss_binomial.pdf", width=9, height=12, pointsize=12, onefile=TRUE)
print(large_ss_binomial)
dev.off()


# Different alphaUP and alphaLOW ------------------------------------------

#Normal outcomes - small sample size scenarios (1,3,4,9,10)
outcome='normal'
alpha.up=0.1
alpha.low=0.001

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

out.plot=NULL
t=1

for (z in c(1,3,4,9,10))
{
  
  load(file=paste0(dir,setups$outcome[z],"_n0c:",setups$n0c[z],
                   "_nc:",setups$nc[z],
                   "_nt:",setups$nt[z],
                   "_alpha.up:",setups$alpha.up[z],
                   "_alpha.low:",setups$alpha.low[z],
                   "_cp:",cp,
                   "_ct:",ct,
                   "_w.mix:",setups$w.mix[z],"_OCdat.RData",sep=""))
  
  
  
  out$Decision <- factor(out$Decision, levels = c("CD-Constraint", "CD-Discard",
                                                  "RMD-Unit","EBPowD","BD","FD"))
  
  #TIE rate & Power (uncalibrated)
  
  t1e <- theme +
    geom_segment(aes(x=-2, xend=2,y = out$alpha.up[1],yend = out$alpha.up[1]), col='red', linetype='dashed', size=0.3)+
    geom_segment(aes(x=-2, xend=2,y = out$alpha.low[1],yend = out$alpha.low[1]), col='red', linetype='dashed', size=0.3)+ 
    geom_line(aes(Conflict, Freq.type.I,  color=Decision),size=0.5, data=out) + 
    ylab("TIE rate") +  ylim(c(0,0.3)) + guides(fill = "none")+
    theme(legend.title=element_blank()) + 
    ggtitle(TeX(sprintf("$n_{0C} = %s,\\; n_{C} = %s,\\; n_{T} = %s$",
                        setups$n0c[z],
                        setups$nc[z],
                        setups$nt[z])))+
    theme( plot.title = element_text(size = 8))
  
  power <- theme+ 
    geom_line(aes(Conflict, Freq.Power,  color=Decision),size=0.5, data=out) + 
    ylab(TeX('Power')) + ylim(c(0.2,1)) + 
    theme(legend.position = "none") + 
    ggtitle(TeX(sprintf("$n_{0C} = %s,\\; n_{C} = %s,\\; n_{T} = %s,\\; \\delta = %s$",
                        setups$n0c[z],
                        setups$nc[z],
                        setups$nt[z],
                        round(out$delta[1],2)))) +
    theme(plot.title = element_text(size = 8)) 
  
  out.plot[[t]]= t1e
  out.plot[[t+1]]= power
  t=t+2
}

small_ss_normal_wider = (
  ((out.plot[[1]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
     (out.plot[[2]] + guides(color = "none"))) /
    ((out.plot[[3]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[4]] + guides(color = "none"))) / 
    ((out.plot[[5]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[6]] + guides(color = "none"))) /
    ((out.plot[[7]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[8]] + guides(color = "none"))) /
    ((out.plot[[9]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[10]] + guides(color = "none"))))  +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

pdf(file ="small_ss_normal_wider.pdf", width=9, height=12, pointsize=12, onefile=TRUE)
print(small_ss_normal_wider)
dev.off()

#Additional large sample size scenarios (2,6,7,8)

out.plot=NULL
t=1

for (z in c(2,5,6,7,8))
{
  
  load(file=paste0(dir,setups$outcome[z],"_n0c:",setups$n0c[z],
                   "_nc:",setups$nc[z],
                   "_nt:",setups$nt[z],
                   "_alpha.up:",setups$alpha.up[z],
                   "_alpha.low:",setups$alpha.low[z],
                   "_cp:",cp,
                   "_ct:",ct,
                   "_w.mix:",setups$w.mix[z],"_OCdat.RData",sep=""))
  
  
  
  out$Decision <- factor(out$Decision, levels = c("CD-Constraint", "CD-Discard",
                                                  "RMD-Unit","EBPowD","BD","FD"))
  
  #TIE rate & Power (uncalibrated)
  
  t1e <- theme +
    geom_segment(aes(x=-2, xend=2,y = out$alpha.up[1],yend = out$alpha.up[1]), col='red', linetype='dashed', size=0.3)+
    geom_segment(aes(x=-2, xend=2,y = out$alpha.low[1],yend = out$alpha.low[1]), col='red', linetype='dashed', size=0.3)+ 
    geom_line(aes(Conflict, Freq.type.I,  color=Decision),size=0.5, data=out) + 
    ylab("TIE rate") +  ylim(c(0,0.3)) + guides(fill = "none")+
    theme(legend.title=element_blank()) + 
    ggtitle(TeX(sprintf("$n_{0C} = %s,\\; n_{C} = %s,\\; n_{T} = %s$",
                        setups$n0c[z],
                        setups$nc[z],
                        setups$nt[z])))+
    theme( plot.title = element_text(size = 8))
  
  power <- theme+ 
    geom_line(aes(Conflict, Freq.Power,  color=Decision),size=0.5, data=out) + 
    ylab(TeX('Power')) + ylim(c(0.2,1)) + 
    theme(legend.position = "none") + 
    ggtitle(TeX(sprintf("$n_{0C} = %s,\\; n_{C} = %s,\\; n_{T} = %s,\\; \\delta = %s$",
                        setups$n0c[z],
                        setups$nc[z],
                        setups$nt[z],
                        round(out$delta[1],2)))) +
    theme(plot.title = element_text(size = 8)) 
  
  out.plot[[t]]= t1e
  out.plot[[t+1]]= power
  t=t+2
}

large_ss_normal_wider = (
  ((out.plot[[1]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
     (out.plot[[2]] + guides(color = "none"))) /
    ((out.plot[[3]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[4]] + guides(color = "none"))) / 
    ((out.plot[[5]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[6]] + guides(color = "none"))) /
    ((out.plot[[7]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[8]] + guides(color = "none"))) /
    ((out.plot[[9]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[10]] + guides(color = "none"))))  +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

pdf(file ="large_ss_normal_wider.pdf", width=9, height=12, pointsize=12, onefile=TRUE)
print(large_ss_normal_wider)
dev.off()


# Binomial outcomes -------------------------------------------------------

outcome='binomial'
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

#Small sample size scenarios (1,3,4,9,10)

out.plot=NULL
t=1

for (z in c(1,3,4,9,10))
{
  
  load(file=paste0(dir,setups$outcome[z],"_n0c:",setups$n0c[z],
                   "_nc:",setups$nc[z],
                   "_nt:",setups$nt[z],
                   "_alpha.up:",setups$alpha.up[z],
                   "_alpha.low:",setups$alpha.low[z],
                   "_cp:",cp,
                   "_ct:",ct,
                   "_w.mix:",setups$w.mix[z],"_OCdat.RData",sep=""))
  
  
  
  out$Decision <- factor(out$Decision, levels = c("CD-Constraint", "CD-Discard",
                                                  "RMD-Unit","EBPowD","BD","FD"))
  
  #TIE rate & Power (uncalibrated)
  
  t1e <- theme +
    geom_segment(aes(x=-out$mupc[1], xend=(1-out$mupc[1]),y = out$alpha.up[1],yend = out$alpha.up[1]), col='red', linetype='dashed', size=0.3)+
    geom_segment(aes(x=-out$mupc[1], xend=(1-out$mupc[1]),y = out$alpha.low[1],yend = out$alpha.low[1]), col='red', linetype='dashed', size=0.3)+ 
    geom_line(aes(Conflict, Freq.type.I,  color=Decision),size=0.5, data=out) + 
    ylab("TIE rate") +  ylim(c(0,0.35)) + guides(fill = "none")+
    theme(legend.title=element_blank()) + 
    ggtitle(TeX(sprintf("$n_{0C} = %s,\\; n_{C} = %s,\\; n_{T} = %s$",
                        setups$n0c[z],
                        setups$nc[z],
                        setups$nt[z])))+
    theme( plot.title = element_text(size = 8))
  
  power <- theme+ 
    geom_line(aes(Conflict, Freq.Power,  color=Decision),size=0.5, data=out) + 
    ylab(TeX('Power')) + ylim(c(0.4,1)) + xlim(c(-out$mupc[1],1-out$mupc[1]-out$delta[1])) + 
    theme(legend.position = "none") + 
    ggtitle(TeX(sprintf("$n_{0C} = %s,\\; n_{C} = %s,\\; n_{T} = %s,\\; \\delta = %s$",
                        setups$n0c[z],
                        setups$nc[z],
                        setups$nt[z],
                        round(out$delta[1],2)))) +
    theme(plot.title = element_text(size = 8)) 
  
  out.plot[[t]]= t1e
  out.plot[[t+1]]= power
  t=t+2
}

small_ss_binomial_wider = (
  ((out.plot[[1]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
     (out.plot[[2]] + guides(color = "none"))) /
    ((out.plot[[3]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[4]] + guides(color = "none"))) / 
    ((out.plot[[5]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[6]] + guides(color = "none"))) /
    ((out.plot[[7]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[8]] + guides(color = "none"))) /
    ((out.plot[[9]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[10]] + guides(color = "none"))))  +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

pdf(file ="small_ss_binomial_wider.pdf", width=9, height=12, pointsize=12, onefile=TRUE)
print(small_ss_binomial_wider)
dev.off()

#Additional large sample size scenarios (2,5,6,7,8)

out.plot=NULL
t=1

for (z in c(2,5,6,7,8))
{
  
  load(file=paste0(dir,setups$outcome[z],"_n0c:",setups$n0c[z],
                   "_nc:",setups$nc[z],
                   "_nt:",setups$nt[z],
                   "_alpha.up:",setups$alpha.up[z],
                   "_alpha.low:",setups$alpha.low[z],
                   "_cp:",cp,
                   "_ct:",ct,
                   "_w.mix:",setups$w.mix[z],"_OCdat.RData",sep=""))
  
  
  
  out$Decision <- factor(out$Decision, levels = c("CD-Constraint", "CD-Discard",
                                                  "RMD-Unit","EBPowD","BD","FD"))
  
  #TIE rate & Power (uncalibrated)
  
  t1e <- theme +
    geom_segment(aes(x=-out$mupc[1], xend=(1-out$mupc[1]),y = out$alpha.up[1],yend = out$alpha.up[1]), col='red', linetype='dashed', size=0.3)+
    geom_segment(aes(x=-out$mupc[1], xend=(1-out$mupc[1]),y = out$alpha.low[1],yend = out$alpha.low[1]), col='red', linetype='dashed', size=0.3)+ 
    geom_line(aes(Conflict, Freq.type.I,  color=Decision),size=0.5, data=out) + 
    ylab("TIE rate") +  ylim(c(0,0.35)) + guides(fill = "none")+
    theme(legend.title=element_blank()) + 
    ggtitle(TeX(sprintf("$n_{0C} = %s,\\; n_{C} = %s,\\; n_{T} = %s$",
                        setups$n0c[z],
                        setups$nc[z],
                        setups$nt[z])))+
    theme( plot.title = element_text(size = 8))
  
  power <- theme+ 
    geom_line(aes(Conflict, Freq.Power,  color=Decision),size=0.5, data=out) + 
    ylab(TeX('Power')) + ylim(c(0.4,1)) + xlim(c(-out$mupc[1],1-out$mupc[1]-out$delta[1])) +  
    theme(legend.position = "none") + 
    ggtitle(TeX(sprintf("$n_{0C} = %s,\\; n_{C} = %s,\\; n_{T} = %s,\\; \\delta = %s$",
                        setups$n0c[z],
                        setups$nc[z],
                        setups$nt[z],
                        round(out$delta[1],2)))) +
    theme(plot.title = element_text(size = 8)) 
  
  out.plot[[t]]= t1e
  out.plot[[t+1]]= power
  t=t+2
}

large_ss_binomial_wider = (
  ((out.plot[[1]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
     (out.plot[[2]] + guides(color = "none"))) /
    ((out.plot[[3]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[4]] + guides(color = "none"))) / 
    ((out.plot[[5]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[6]] + guides(color = "none"))) /
    ((out.plot[[7]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[8]] + guides(color = "none"))) /
    ((out.plot[[9]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[10]] + guides(color = "none"))))  +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

pdf(file ="large_ss_binomial_wider.pdf", width=9, height=12, pointsize=12, onefile=TRUE)
print(large_ss_binomial_wider)
dev.off()



# Different t ------------------------------------------

#Normal outcomes - small sample size scenarios (1,3,4,9,10)
outcome='normal'
alpha.up=0.075
alpha.low=0.01
ct=2

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

out.plot=NULL
t=1

for (z in c(1,3,4,9,10))
{
  
  load(file=paste0(dir,setups$outcome[z],"_n0c:",setups$n0c[z],
                   "_nc:",setups$nc[z],
                   "_nt:",setups$nt[z],
                   "_alpha.up:",setups$alpha.up[z],
                   "_alpha.low:",setups$alpha.low[z],
                   "_cp:",cp,
                   "_ct:",ct,
                   "_w.mix:",setups$w.mix[z],"_OCdat.RData",sep=""))
  
  
  
  out$Decision <- factor(out$Decision, levels = c("CD-Constraint", "CD-Discard",
                                                  "RMD-Unit","EBPowD","BD","FD"))
  
  #TIE rate & Power (uncalibrated)
  
  t1e <- theme +
    geom_segment(aes(x=-2, xend=2,y = out$alpha.up[1],yend = out$alpha.up[1]), col='red', linetype='dashed', size=0.3)+
    geom_segment(aes(x=-2, xend=2,y = out$alpha.low[1],yend = out$alpha.low[1]), col='red', linetype='dashed', size=0.3)+ 
    geom_line(aes(Conflict, Freq.type.I,  color=Decision),size=0.5, data=out) + 
    ylab("TIE rate") +  ylim(c(0,0.3)) + guides(fill = "none")+
    theme(legend.title=element_blank()) + 
    ggtitle(TeX(sprintf("$n_{0C} = %s,\\; n_{C} = %s,\\; n_{T} = %s$",
                        setups$n0c[z],
                        setups$nc[z],
                        setups$nt[z])))+
    theme( plot.title = element_text(size = 8))
  
  power <- theme+ 
    geom_line(aes(Conflict, Freq.Power,  color=Decision),size=0.5, data=out) + 
    ylab(TeX('Power')) + ylim(c(0.2,1)) + 
    theme(legend.position = "none") + 
    ggtitle(TeX(sprintf("$n_{0C} = %s,\\; n_{C} = %s,\\; n_{T} = %s,\\; \\delta = %s$",
                        setups$n0c[z],
                        setups$nc[z],
                        setups$nt[z],
                        round(out$delta[1],2)))) +
    theme(plot.title = element_text(size = 8)) 
  
  out.plot[[t]]= t1e
  out.plot[[t+1]]= power
  t=t+2
}

small_ss_normal_t2 = (
  ((out.plot[[1]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
     (out.plot[[2]] + guides(color = "none"))) /
    ((out.plot[[3]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[4]] + guides(color = "none"))) / 
    ((out.plot[[5]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[6]] + guides(color = "none"))) /
    ((out.plot[[7]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[8]] + guides(color = "none"))) /
    ((out.plot[[9]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[10]] + guides(color = "none"))))  +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

pdf(file ="small_ss_normal_t2.pdf", width=9, height=12, pointsize=12, onefile=TRUE)
print(small_ss_normal_t2)
dev.off()

#Additional large sample size scenarios (2,6,7,8)

out.plot=NULL
t=1

for (z in c(2,5,6,7,8))
{
  
  load(file=paste0(dir,setups$outcome[z],"_n0c:",setups$n0c[z],
                   "_nc:",setups$nc[z],
                   "_nt:",setups$nt[z],
                   "_alpha.up:",setups$alpha.up[z],
                   "_alpha.low:",setups$alpha.low[z],
                   "_cp:",cp,
                   "_ct:",ct,
                   "_w.mix:",setups$w.mix[z],"_OCdat.RData",sep=""))
  
  
  
  out$Decision <- factor(out$Decision, levels = c("CD-Constraint", "CD-Discard",
                                                  "RMD-Unit","EBPowD","BD","FD"))
  
  #TIE rate & Power (uncalibrated)
  
  t1e <- theme +
    geom_segment(aes(x=-2, xend=2,y = out$alpha.up[1],yend = out$alpha.up[1]), col='red', linetype='dashed', size=0.3)+
    geom_segment(aes(x=-2, xend=2,y = out$alpha.low[1],yend = out$alpha.low[1]), col='red', linetype='dashed', size=0.3)+ 
    geom_line(aes(Conflict, Freq.type.I,  color=Decision),size=0.5, data=out) + 
    ylab("TIE rate") +  ylim(c(0,0.3)) + guides(fill = "none")+
    theme(legend.title=element_blank()) + 
    ggtitle(TeX(sprintf("$n_{0C} = %s,\\; n_{C} = %s,\\; n_{T} = %s$",
                        setups$n0c[z],
                        setups$nc[z],
                        setups$nt[z])))+
    theme( plot.title = element_text(size = 8))
  
  power <- theme+ 
    geom_line(aes(Conflict, Freq.Power,  color=Decision),size=0.5, data=out) + 
    ylab(TeX('Power')) + ylim(c(0.2,1)) + 
    theme(legend.position = "none") + 
    ggtitle(TeX(sprintf("$n_{0C} = %s,\\; n_{C} = %s,\\; n_{T} = %s,\\; \\delta = %s$",
                        setups$n0c[z],
                        setups$nc[z],
                        setups$nt[z],
                        round(out$delta[1],2)))) +
    theme(plot.title = element_text(size = 8)) 
  
  out.plot[[t]]= t1e
  out.plot[[t+1]]= power
  t=t+2
}

large_ss_normal_t2 = (
  ((out.plot[[1]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
     (out.plot[[2]] + guides(color = "none"))) /
    ((out.plot[[3]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[4]] + guides(color = "none"))) / 
    ((out.plot[[5]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[6]] + guides(color = "none"))) /
    ((out.plot[[7]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[8]] + guides(color = "none"))) /
    ((out.plot[[9]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[10]] + guides(color = "none"))))  +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

pdf(file ="large_ss_normal_t2.pdf", width=9, height=12, pointsize=12, onefile=TRUE)
print(large_ss_normal_t2)
dev.off()

# Binomial outcomes -------------------------------------------------------

outcome='binomial'
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

#Small sample size scenarios (1,3,4,9,10)

out.plot=NULL
t=1

for (z in c(1,3,4,9,10))
{
  
  load(file=paste0(dir,setups$outcome[z],"_n0c:",setups$n0c[z],
                   "_nc:",setups$nc[z],
                   "_nt:",setups$nt[z],
                   "_alpha.up:",setups$alpha.up[z],
                   "_alpha.low:",setups$alpha.low[z],
                   "_cp:",cp,
                   "_ct:",ct,
                   "_w.mix:",setups$w.mix[z],"_OCdat.RData",sep=""))
  
  
  
  out$Decision <- factor(out$Decision, levels = c("CD-Constraint", "CD-Discard",
                                                  "RMD-Unit","EBPowD","BD","FD"))
  
  #TIE rate & Power (uncalibrated)
  
  t1e <- theme +
    geom_segment(aes(x=-out$mupc[1], xend=(1-out$mupc[1]),y = out$alpha.up[1],yend = out$alpha.up[1]), col='red', linetype='dashed', size=0.3)+
    geom_segment(aes(x=-out$mupc[1], xend=(1-out$mupc[1]),y = out$alpha.low[1],yend = out$alpha.low[1]), col='red', linetype='dashed', size=0.3)+ 
    geom_line(aes(Conflict, Freq.type.I,  color=Decision),size=0.5, data=out) + 
    ylab("TIE rate") +  ylim(c(0,0.35)) + guides(fill = "none")+
    theme(legend.title=element_blank()) + 
    ggtitle(TeX(sprintf("$n_{0C} = %s,\\; n_{C} = %s,\\; n_{T} = %s$",
                        setups$n0c[z],
                        setups$nc[z],
                        setups$nt[z])))+
    theme( plot.title = element_text(size = 8))
  
  power <- theme+ 
    geom_line(aes(Conflict, Freq.Power,  color=Decision),size=0.5, data=out) + 
    ylab(TeX('Power')) + ylim(c(0.4,1)) + xlim(c(-out$mupc[1],1-out$mupc[1]-out$delta[1])) + 
    theme(legend.position = "none") + 
    ggtitle(TeX(sprintf("$n_{0C} = %s,\\; n_{C} = %s,\\; n_{T} = %s,\\; \\delta = %s$",
                        setups$n0c[z],
                        setups$nc[z],
                        setups$nt[z],
                        round(out$delta[1],2)))) +
    theme(plot.title = element_text(size = 8)) 
  
  out.plot[[t]]= t1e
  out.plot[[t+1]]= power
  t=t+2
}

small_ss_binomial_t2 = (
  ((out.plot[[1]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
     (out.plot[[2]] + guides(color = "none"))) /
    ((out.plot[[3]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[4]] + guides(color = "none"))) / 
    ((out.plot[[5]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[6]] + guides(color = "none"))) /
    ((out.plot[[7]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[8]] + guides(color = "none"))) /
    ((out.plot[[9]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[10]] + guides(color = "none"))))  +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

pdf(file ="small_ss_binomial_t2.pdf", width=9, height=12, pointsize=12, onefile=TRUE)
print(small_ss_binomial_t2)
dev.off()

#Additional large sample size scenarios (2,5,6,7,8)

out.plot=NULL
t=1

for (z in c(2,5,6,7,8))
{
  
  load(file=paste0(dir,setups$outcome[z],"_n0c:",setups$n0c[z],
                   "_nc:",setups$nc[z],
                   "_nt:",setups$nt[z],
                   "_alpha.up:",setups$alpha.up[z],
                   "_alpha.low:",setups$alpha.low[z],
                   "_cp:",cp,
                   "_ct:",ct,
                   "_w.mix:",setups$w.mix[z],"_OCdat.RData",sep=""))
  
  
  
  out$Decision <- factor(out$Decision, levels = c("CD-Constraint", "CD-Discard",
                                                  "RMD-Unit","EBPowD","BD","FD"))
  
  #TIE rate & Power (uncalibrated)
  
  t1e <- theme +
    geom_segment(aes(x=-out$mupc[1], xend=(1-out$mupc[1]),y = out$alpha.up[1],yend = out$alpha.up[1]), col='red', linetype='dashed', size=0.3)+
    geom_segment(aes(x=-out$mupc[1], xend=(1-out$mupc[1]),y = out$alpha.low[1],yend = out$alpha.low[1]), col='red', linetype='dashed', size=0.3)+ 
    geom_line(aes(Conflict, Freq.type.I,  color=Decision),size=0.5, data=out) + 
    ylab("TIE rate") +  ylim(c(0,0.35)) + guides(fill = "none")+
    theme(legend.title=element_blank()) + 
    ggtitle(TeX(sprintf("$n_{0C} = %s,\\; n_{C} = %s,\\; n_{T} = %s$",
                        setups$n0c[z],
                        setups$nc[z],
                        setups$nt[z])))+
    theme( plot.title = element_text(size = 8))
  
  power <- theme+ 
    geom_line(aes(Conflict, Freq.Power,  color=Decision),size=0.5, data=out) + 
    ylab(TeX('Power')) + ylim(c(0.4,1))  + xlim(c(-out$mupc[1],1-out$mupc[1]-out$delta[1])) + 
    theme(legend.position = "none") + 
    ggtitle(TeX(sprintf("$n_{0C} = %s,\\; n_{C} = %s,\\; n_{T} = %s,\\; \\delta = %s$",
                        setups$n0c[z],
                        setups$nc[z],
                        setups$nt[z],
                        round(out$delta[1],2)))) +
    theme(plot.title = element_text(size = 8)) 
  
  out.plot[[t]]= t1e
  out.plot[[t+1]]= power
  t=t+2
}

large_ss_binomial_t2 = (
  ((out.plot[[1]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
     (out.plot[[2]] + guides(color = "none"))) /
    ((out.plot[[3]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[4]] + guides(color = "none"))) / 
    ((out.plot[[5]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[6]] + guides(color = "none"))) /
    ((out.plot[[7]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[8]] + guides(color = "none"))) /
    ((out.plot[[9]] + guides(col = guide_legend(nrow = 1, byrow = TRUE))) |
       (out.plot[[10]] + guides(color = "none"))))  +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

pdf(file ="large_ss_binomial_t2.pdf", width=9, height=12, pointsize=12, onefile=TRUE)
print(large_ss_binomial_t2)
dev.off()
