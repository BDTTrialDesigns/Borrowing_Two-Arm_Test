#This script can be used to reproduce the calculations and plots of the
#Normal and Binomial outcomes case-studies.

library(RBesT)
library(parallel)
library(ggplot2)
library(patchwork)

rm(list = ls())

setwd("~/")

source("Functions_Subm_N.R")

#####Case study - Best et al., 2025

#parameters
outcome = 'normal'
sig = 88     #variance of one observation

alpha.up = 0.07         #max TIE 80% prior CI
alpha.low = 0.002         #min TIE 80% prior CI
alpha.b = 0.025          #frequentist TIE
th0 = 0                  #null hypothesis border
nc = 20
nt = 40
conflict = -seq(-200, 0, 5) - 50
delta = 70
ct = 3
cp = 4


colors = c(
  "FD" = "#D55E00",
  "BD" = "#0072B2",
  "MAPD" = "#56B4E9",
  "EBPowD" = "#009E73",
  "MAPD-Constraint" = "#CC79A7",
  "RMAPD" = "#004488",         
  "RMAPD-Constraint" = "#882255",
  "CD-Discard" = "#7B3294",
  "CD-Discard (1)" = "#1B9E77",
  "CD-Discard (2)" = "#A6761D",
  "CD-Discard (3)" = "#E7298A",
  "CD-Discard (4)" = "#7B3294"
)

#simulation parameters
nmc = 2 * 10 ^ 5             #number of Monte Carlo simulations
seed = 1234

#robust prior specification (MAP and robust MAP)

#MAP prior approximation from Best et al.
mix.control = mixnorm(
  info1 = c(0.51, 51, 19.9),
  info2 = c(0.44, 46.8, 7.6),
  info3 = c(0.05, 54.1, 51.7),
  sigma = sig
)
mix.robust = mixnorm(c(1, 50, 88))
mix.control.rob = mixcombine(mix.control, mix.robust, weight = c(0.8, 0.2))

mix.treat = mixnorm(info1 = c(1, 50, 8800 ^ 2), sigma = sig)

#obtain test decision threshold to be applied to frequentist test to obtain the same test decisions
#as for the Bayes test under the mixture prior: this is needed to apply the TIE rate cap

ycv = -seq(-650, 500, 1)
deltav = -seq(-200, 200, 1)
eval.grid = expand.grid(ycv, deltav)

yc = eval.grid[, 1]
dobs = eval.grid[, 2]
yt = eval.grid[, 1] + eval.grid[, 2]

z = dobs / sqrt((sig ^ 2 / nc) + ((sig ^ 2 / nt)))

p.h1.mix.map = pH1(
  prior = 'mix',
  yc = yc,
  yt = yt,
  nc = nc,
  nt = nt,
  sigma = sig,
  priorPars = NULL,
  w.mix = w.mix,
  mix.control = mix.control,
  mix.treat = mix.treat,
  outcome = 'normal',
  th0 = th0
)

dec.map = (p.h1.mix.map > (1 - alpha.b))

p.h1.mix.map.robust = pH1(
  prior = 'mix',
  yc = yc,
  yt = yt,
  nc = nc,
  nt = nt,
  sigma = sig,
  priorPars = NULL,
  w.mix = w.mix,
  mix.control = mix.control.rob,
  mix.treat = mix.treat,
  outcome = 'normal',
  th0 = th0
)

dec.map.robust = (p.h1.mix.map.robust > (1 - alpha.b))

dat.dec = data.frame(
  zval = z,
  yc = yc,
  p.h1 = p.h1.mix.map,
  dec = dec.map,
  dec.robust = dec.map.robust
)

#obtain critical values
zstat.map = tapply(dat.dec$zval * dat.dec$dec, list(as.factor(dat.dec$yc)), function(x)
  x[abs(x) > 0][which.min(x[abs(x) > 0])])

zstat.map.robust = tapply(dat.dec$zval * dat.dec$dec.robust, list(as.factor(dat.dec$yc)), function(x)
  x[abs(x) > 0][which.min(x[abs(x) > 0])])

#obtain test decision threshold and record grid (for interpolation)
kappa.cal.mix.map = as.vector(1 - pnorm(t(zstat.map)))
kappa.cal.grid.map = as.numeric(names(zstat.map))

kappa.cal.mix.map.robust = as.vector(1 - pnorm(t(zstat.map.robust)))
kappa.cal.grid.map.robust = as.numeric(names(zstat.map.robust))

#corresponding critical value
z.constraint.map = qnorm(1 - pmax(pmin(kappa.cal.mix.map, alpha.up), alpha.low))
z.constraint.map.robust = qnorm(1 - pmax(pmin(kappa.cal.mix.map.robust, alpha.up), alpha.low))

#plot critical values

boundPlot = data.frame(
  grid = rep(-kappa.cal.grid.map, 5),
  zval = c(zstat.map, z.constraint.map, zstat.map.robust, z.constraint.map.robust, rep(
    qnorm(1 - alpha.b), length(kappa.cal.grid.map)
  )),
  Decision = rep(c("MAPD", "MAPD-Constraint","RMAPD", "RMAPD-Constraint", "FD"), each =
                   length(kappa.cal.grid.map))
)

boundPlot$Decision <- factor(boundPlot$Decision, levels = c("MAPD", "MAPD-Constraint",
                                                            "RMAPD", "RMAPD-Constraint", 
                                                            "FD"))

z.MAP = ggplot() +   theme_minimal() +
  geom_hline(
    yintercept = qnorm(1 - alpha.up),
    linetype = 'dashed',
    color = 'red',
    size=0.25
  ) +
  geom_hline(
    yintercept = qnorm(1 - alpha.low),
    linetype = 'dashed',
    color = 'red',
    size=0.25
  ) +
  geom_contour_filled(aes(x = -yc, y = z, z =(1-p.h1)*100),
                      data=dat.dec,
                      breaks = (1-c(0,0.95,0.975,0.99,1))*100, alpha=0.2) + 
  scale_fill_viridis_d(name=TeX("$P^{MAP}(H_0|y) \\%")) +
  geom_line(data = boundPlot,
            aes(x = grid, y = zval, color = Decision),
            size = 0.5) +
  geom_vline(xintercept = -50,
             linetype = 'dashed',
             color = 'grey40', size=0.25)+
  xlab(TeX("$\\bar{y}_C$")) +
  ylab(TeX("$z_{1-\\kappa(\\bar{y}_{C})}$")) +
  ylim(c(0, 4)) +
  xlim(c(-150, 50)) +
  scale_colour_manual(values = colors) +
  theme(legend.title = element_text(),
        legend.position = "top",
        plot.title = element_text(size = 18)
  ) +
  guides(
    colour = guide_legend(title = NULL, nrow = 2),
    fill = guide_legend(position = "bottom", title.theme = element_text(size = 8),
                        label.theme = element_text(size = 6),nrow = 1),
  )


#informative (non-robust) prior
#use study Gastr06 in crohn dataset

priorPars = rbind(c(-crohn[1, ]$y, -crohn[1, ]$y), c(sig / sqrt(crohn[1, ]$n), Inf))
cp = cp
cf = ct * sqrt(sig ^ 2 / nc + sig ^ 2 / crohn[1, ]$n) #point of full discard

#obtain critical values for different p

Ac = sig ^ 2 / (nc * (sig / sqrt(crohn[1, ]$n)) ^ 2)
kappa.cal = kappa.cal.f(
  yc = kappa.cal.grid.map,
  alpha.b = alpha.b,
  At = 0,
  Ac = Ac,
  mupc = (-crohn[1, ]$y),
  mupt = (-crohn[1, ]$y),
  sigma = sig,
  nt = nt,
  nc = nc,
  th0 = th0
)
kappa1 = pmax(pmin(
  kappa.discount(
    alpha.k = alpha.b,
    alpha.k2 = kappa.cal,
    eval = (kappa.cal.grid.map),
    muc = (-crohn[1, ]$y),
    cf = cf,
    cp = 1
  ),
  alpha.up
), alpha.low)#,
kappa2 = pmax(pmin(
  kappa.discount(
    alpha.k = alpha.b,
    alpha.k2 = kappa.cal,
    eval = (kappa.cal.grid.map),
    muc = (-crohn[1, ]$y),
    cf = cf,
    cp = 2
  ),
  alpha.up
), alpha.low)#,
kappa3 = pmax(pmin(
  kappa.discount(
    alpha.k = alpha.b,
    alpha.k2 = kappa.cal,
    eval = (kappa.cal.grid.map),
    muc = (-crohn[1, ]$y),
    cf = cf,
    cp = 3
  ),
  alpha.up
), alpha.low)#,
kappa4 = pmax(pmin(
  kappa.discount(
    alpha.k = alpha.b,
    alpha.k2 = kappa.cal,
    eval = (kappa.cal.grid.map),
    muc = (-crohn[1, ]$y),
    cf = cf,
    cp = 4
  ),
  alpha.up
), alpha.low)#,

boundPlotNonRobust = data.frame(
  grid = rep(-kappa.cal.grid.map, 6),
  zval = c(
    qnorm(1 - kappa.cal),
    qnorm(1 - kappa1),
    qnorm(1 - kappa2),
    qnorm(1 - kappa3),
    qnorm(1 - kappa4),
    rep(qnorm(1 - alpha.b), length(kappa.cal.grid.map))
  ),
  Decision = rep(
    c(
      "BD",
      "CD-Discard (1)",
      "CD-Discard (2)",
      "CD-Discard (3)",
      "CD-Discard (4)",
      "FD"
    ),
    each = length(kappa.cal.grid.map)
  )
)

boundPlotNonRobust$Decision <- factor(
  boundPlotNonRobust$Decision,
  levels = c(
    "BD",
    "CD-Discard (1)",
    "CD-Discard (2)",
    "CD-Discard (3)",
    "CD-Discard (4)",
    "FD"
  )
)

p.h1.info = pH1(
  prior = 'info',
  yc = yc,
  yt = yt,
  nc = nc,
  nt = nt,
  sigma = sig,
  priorPars = priorPars,
  outcome = 'normal',
  th0 = th0
)

dec = (p.h1.info > (1 - alpha.b))

dat.dec.info = data.frame(
  zval = z,
  yc = yc,
  p.h1 = p.h1.info,
  dec = dec
)

z.NonRobust = ggplot() +   theme_minimal() + 
  geom_contour_filled(aes(x = -yc, y = z, z =(1-p.h1)*100),
                      data=dat.dec.info,
                      breaks = (1-c(0,0.95,0.975,0.99,1))*100, alpha=0.2) + 
  scale_fill_viridis_d(name=TeX("$P^{\\pi}(H_0|y) \\%"))+
  geom_line(data = boundPlotNonRobust,
            aes(x = grid, y = zval, color = Decision),
            size = 0.5) +
  geom_vline(xintercept = -50,
             linetype = 'dashed',
             color = 'grey40', size=0.25) +
  geom_hline(
    yintercept = qnorm(1 - alpha.up),
    linetype = 'dashed',
    color = 'red',
    size=0.25
  ) +
  geom_hline(
    yintercept = qnorm(1 - alpha.low),
    linetype = 'dashed',
    color = 'red',
    size=0.25
  ) +
  xlab(TeX("$\\bar{y}_C$")) +
  ylab(TeX("$z_{1-\\kappa(\\bar{y}_{C})}$")) +
  ylim(c(0, 4)) +
  xlim(c(-150, 50)) +
  scale_colour_manual(values = colors)  +
  scale_colour_manual(values = colors) +
  theme(legend.title = element_text(),
        legend.position = "top",
        plot.title = element_text(size = 18)
  ) +
  guides(
    colour = guide_legend(title = NULL, nrow = 2),
    fill = guide_legend(position = "bottom", title.theme = element_text(size = 8),
                        label.theme = element_text(size = 6),nrow = 1),
  )

dec.thr = (z.MAP + guides(col = guide_legend(nrow = 2, byrow = T)) + theme(legend.title =
                                                                             element_blank()) |
             z.NonRobust)

#compute OCs, under vague prior, mixture prior, mixture prior + cap (CD constrain),
#informative prior, informative prior + cap and discard (CD Discard)

ocs = compute.oc(
  outcome = outcome,
  nt = nt,
  nc = nc,
  delta = delta,
  conflict = conflict-50,
  priorPars = priorPars,
  sigma = sig,
  alpha.low = alpha.low,
  alpha.up = alpha.up,
  alpha.b = alpha.b,
  mix.control = mix.control,
  mix.treat = mix.treat,
  kappa.cal = kappa.cal.mix.map,
  kappa.cal.grid = kappa.cal.grid.map,
  cp = cp,
  cf = cf,
  nmc = nmc,
  seed = 1234,
  ncores = 28
)


Decision = colnames(ocs$TIE)
Decision[c(3, 5, 6)] = c("MAPD", "CD-Discard (4)", "MAPD-Constraint")

out.map = data.frame(
  'Freq.type.I' = as.vector(ocs$TIE),
  'Freq.Power' = as.vector(ocs$Pow),
  'Type.I.cal' = as.vector(ocs$TIEcal),
  'Power.cal' = as.vector(ocs$Powcal),
  'Conflict' = rep(-as.numeric(rownames(ocs$TIE))-50, ncol(ocs$TIE)),
  'Decision' = rep(Decision, each = nrow(ocs$TIE))
)

#add robust map

ocs.robust = compute.oc(
  outcome = outcome,
  nt = nt,
  nc = nc,
  delta = delta,
  conflict = conflict-50,
  priorPars = priorPars,
  sigma = sig,
  alpha.low = alpha.low,
  alpha.up = alpha.up,
  alpha.b = alpha.b,
  mix.control = mix.control.rob,
  mix.treat = mix.treat,
  kappa.cal = kappa.cal.mix.map.robust,
  kappa.cal.grid = kappa.cal.grid.map.robust,
  cp = cp,
  cf = cf,
  nmc = nmc,
  seed = 1234,
  ncores = 28
)

Decision = colnames(ocs.robust$TIE)
Decision[c(3, 6)] = c("RMAPD", "RMAPD-Constraint")

out.robust = data.frame(
  'Freq.type.I' = as.vector(ocs.robust$TIE),
  'Freq.Power' = as.vector(ocs.robust$Pow),
  'Type.I.cal' = as.vector(ocs.robust$TIEcal),
  'Power.cal' = as.vector(ocs.robust$Powcal),
  'Conflict' = rep(-as.numeric(rownames(ocs.robust$TIE))-50, ncol(ocs.robust$TIE)),
  'Decision' = rep(Decision, each = nrow(ocs.robust$TIE))
)

out=rbind(out.map,out.robust[out.robust$Decision=='RMAPD' | 
                           out.robust$Decision=='RMAPD-Constraint', ])

#find alphaUP and alphaLOW
# limq=qmix(mix.control, c(0.1, 0.9))
# f = approxfun(x=out$Conflict[out$Decision=='MAPD'],y=out$Freq.type.I[out$Decision=='MAPD'],rule=1,method='linear')
# f(-limq)
# 
# limq=qnorm(c(0.1, 0.9),crohn[1, ]$y,sig / sqrt(crohn[1, ]$n))
# f = approxfun(x=out$Conflict[out$Decision=='BD'],y=out$Freq.type.I[out$Decision=='BD'],rule=1,method='linear')
# f(limq)

save.image("~/results/case_study_OCdat.RData")
##

theme = ggplot() + theme_light() +
  scale_colour_manual(values = colors) +
  xlab(expression(theta[C])) #+

out$Decision <- factor(out$Decision,
                       levels = c("MAPD-Constraint","RMAPD-Constraint", "CD-Discard (4)", 
                                  "MAPD","RMAPD","EBPowD", "BD", "FD"))


g1 <- theme +
  geom_segment(
    aes(
      x = -150,
      xend = 50,
      y = alpha.up,
      yend = alpha.up
    ),
    col = 'red',
    linetype = 'dashed',
    size = 0.25
  ) +
  geom_segment(
    aes(
      x = -150,
      xend = 50,
      y = alpha.low,
      yend = alpha.low
    ),
    ,
    col = 'red',
    linetype = 'dashed',
    size = 0.25
  ) +
  geom_line(
    aes(Conflict, Freq.type.I, color = Decision),
    size = 0.5,
    data = subset(
      out,
      out$Decision == 'FD' | out$Decision == 'MAPD' |
        out$Decision =='MAPD-Constraint'|
        out$Decision =='RMAPD-Constraint'|
        out$Decision =='RMAPD'
    )
  ) +
  ylab("TIE rate") +
  geom_vline(xintercept = -50,
             linetype = 'dashed',
             color = 'grey40', size=0.25) +
  ylim(c(0, 0.25)) + guides(fill = "none",
   colour = guide_legend(title = NULL, nrow = 2))

g2 <- theme + geom_line(
  aes(Conflict, Freq.Power, color = Decision),
  size = 0.5,
  data = subset(
    out,
    out$Decision == 'FD' | out$Decision == 'MAPD' |
      out$Decision =='MAPD-Constraint'|
      out$Decision =='RMAPD-Constraint'|
      out$Decision =='RMAPD'
  )
) + ylab(TeX('Power')) +
  geom_vline(xintercept = -50,
             linetype = 'dashed',
             color = 'grey40', size=0.25) +
  ylim(c(0.3, 1)) + theme(legend.position = "none")


g3 <- theme +
  geom_segment(
    aes(
      x = -150,
      xend = 50,
      y = alpha.up,
      yend = alpha.up
    ),
    col = 'red',
    linetype = 'dashed',
    size = 0.25
  ) +
  geom_segment(
    aes(
      x = -150,
      xend = 50,
      y = alpha.low,
      yend = alpha.low
    ),
    ,
    col = 'red',
    linetype = 'dashed',
    size = 0.25
  ) +
  geom_line(
    aes(Conflict, Freq.type.I, color = Decision),
    size = 0.5,
    data = subset(
      out,
      out$Decision == 'FD' | out$Decision == 'BD' |
        out$Decision ==
        'EBPowD' |
        out$Decision ==
        'CD-Discard (4)'
    )
  ) +
  geom_vline(xintercept = -50,
             linetype = 'dashed',
             color = 'grey40', size=0.25) +
  ylab("TIE rate") +
  ylim(c(0, 0.25)) + guides(fill = "none")

g4 <- theme + geom_line(
  aes(Conflict, Freq.Power, color = Decision),
  size = 0.5,
  data = subset(
    out,
    out$Decision == 'FD' | out$Decision == 'BD' |
      out$Decision == 'EBPowD' |
      out$Decision == 'CD-Discard (4)'
  )
) +
  geom_vline(xintercept = -50,
             linetype = 'dashed',
             color = 'grey40', size=0.25) +
  ylab(TeX('Power')) + ylim(c(0.3, 1)) + theme(legend.position = "none")



MAP_oc = (
  (g1 + guides(col = guide_legend(nrow = 2, byrow = T)) + theme(legend.title =
                                                                  element_blank())) / 
    g2 + guides(color = "none")
) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

NonRobust_oc = (
  (g3 + guides(col = guide_legend(nrow = 1, byrow = T)) + theme(legend.title =
                                                                  element_blank())) /
    g4 + guides(color = "none")
) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

all = (dec.thr /
         (MAP_oc |
            NonRobust_oc) +
         plot_layout(heights = c(1, 2))
)

pdf(
  file = paste0("~/results/case_study_dec.thr.pdf"),
  width = 8,
  height = 4,
  pointsize = 12,
  onefile = TRUE
)
dec.thr
dev.off()

pdf(
  file = paste0("~/results/case_study_MAP_OC.pdf"),
  width = 8,
  height = 4,
  pointsize = 12,
  onefile = TRUE
)
print(MAP_oc)
dev.off()

pdf(
  file = paste0("~/results/case_study_NonRobust_OC.pdf"),
  width = 8,
  height = 4,
  pointsize = 12,
  onefile = TRUE
)
print(NonRobust_oc)
dev.off()

pdf(
  file = paste0("~/results/case_study_MAP.pdf"),
  width = 9,
  height = 11,
  pointsize = 14,
  onefile = TRUE
)
print(all)
dev.off()

#####################################################################################
#Binomial outcomes - Viele et al.

rm(list=ls())

setwd("~/")
source("Functions_Subm_N.R")

colors = c(
  "FD" = "#D55E00",
  "BD" = "#0072B2",
  "EBPowD" = "#009E73",
  "RMD-Unit" = "#56B4E9",
  "CD-Constraint" = "#CC79A7",
  "CD-Discard" = "#7B3294",
  "CD-Discard (1)" = "#1B9E77",
  "CD-Discard (2)" = "#A6761D",
  "CD-Discard (3)" = "#E7298A",
  "CD-Discard (4)" = "#7B3294"
)

#####parameters

outcome="binomial"
nt=200      #sample size in treatment arm
nc=200    #sample size in control arm
n0c=100  #control arm prior sample size
mupc=0.65                #prior mean control

cf=4  #point of full discard
cp=4        #steepness of change

alpha.up=0.056         #max TIE (of BD at approx. 97.5% quantile of prior)
alpha.low=0.005       #min TIE (of BD at approx. 2.5% quantile of prior)
alpha.b=0.025         #frequentist TIE
th0=0
priorPars=cbind(c(mupc*n0c,(1-mupc)*n0c),c(0,0))

w.mix=0.3
bound=FALSE

a0t=a0c=b0t=b0c=0.5

#pre-computes post pr of H1 for informative and uniform/jeffreys
pH1array=lapply(c(0,1), function(p) 
{p.h1.precomp(pow=p,nc=nc,nt=nt,n0c=n0c,a0t=a0t,a0c=a0c,b0t=b0t,b0c=b0c,th0=th0)})

#all combinations
yc=0:nc
yt=0:nt
xv=expand.grid(yt,yc)

yt1=xv[,1]
yc1=xv[,2]


p.h1.unif=pH1(prior='vague',yc=yc,yt=yt,nc=nc,nt=nt,
              priorPars=priorPars,w.mix=w.mix,
              outcome='binomial',pH1array=pH1array,
              a0t=a0t,a0c=a0c,b0t=b0t,b0c=b0c,th0=th0)

p.h1.info=pH1(prior='info',yc=yc,yt=yt,nc=nc,nt=nt,
              priorPars=priorPars,w.mix=w.mix,
              outcome='binomial',pH1array=pH1array,
              a0t=a0t,a0c=a0c,b0t=b0t,b0c=b0c,th0=th0)
p.h1.mix=pH1(prior='mix',yc=yc,yt=yt,nc=nc,nt=nt,
             priorPars=priorPars,w.mix=w.mix,
             outcome='binomial',pH1array=pH1array,
             a0t=a0t,a0c=a0c,b0t=b0t,b0c=b0c,th0=th0)
p.h1.EB=pH1(prior='EB',yc=yc,yt=yt,nc=nc,nt=nt,
            priorPars=priorPars,w.mix=w.mix,
            outcome='binomial',pH1array=pH1array,
            bound=bound,a0t=a0t,a0c=a0c,b0t=b0t,b0c=b0c,th0=th0)

kappa.f=cThreshold(yc=yc,yt=yt,nc=nc,nt=nt,
                   priorPars=priorPars,
                   outcome='binomial',
                   alpha.b=alpha.b, alpha.up=alpha.up, pH1array=pH1array,
                   alpha.low=alpha.low, a0t=a0t, a0c=a0c,b0t=b0t,b0c=b0c,
                   th0=th0,cf=cf,cp=4,CD='Discard',
                   kappa.calv=NULL,kappa.cal.grid=NULL)

kappa.f3=cThreshold(yc=yc,yt=yt,nc=nc,nt=nt,
                    priorPars=priorPars,
                    outcome='binomial',
                    alpha.b=alpha.b, alpha.up=alpha.up, pH1array=pH1array,
                    alpha.low=alpha.low, a0t=a0t, a0c=a0c,b0t=b0t,b0c=b0c,
                    th0=th0,cf=cf,cp=3,CD='Discard',
                    kappa.calv=NULL,kappa.cal.grid=NULL)

kappa.f2=cThreshold(yc=yc,yt=yt,nc=nc,nt=nt,
                    priorPars=priorPars,
                    outcome='binomial',
                    alpha.b=alpha.b, alpha.up=alpha.up, pH1array=pH1array,
                    alpha.low=alpha.low, a0t=a0t, a0c=a0c,b0t=b0t,b0c=b0c,
                    th0=th0,cf=cf,cp=2,CD='Discard',
                    kappa.calv=NULL,kappa.cal.grid=NULL)

kappa.f1=cThreshold(yc=yc,yt=yt,nc=nc,nt=nt,
                    priorPars=priorPars,
                    outcome='binomial',
                    alpha.b=alpha.b, alpha.up=alpha.up, pH1array=pH1array,
                    alpha.low=alpha.low, a0t=a0t, a0c=a0c,b0t=b0t,b0c=b0c,
                    th0=th0,cf=cf,cp=1,CD='Discard',
                    kappa.calv=NULL,kappa.cal.grid=NULL)

kappa.f.c=cThreshold(yc=yc,yt=yt,nc=nc,nt=nt,
                     priorPars=priorPars,
                     outcome='binomial',
                     alpha.b=alpha.b, alpha.up=alpha.up, pH1array=pH1array,
                     alpha.low=alpha.low, a0t=a0t, a0c=a0c,b0t=b0t,b0c=b0c,
                     th0=th0,
                     cf=cf,cp=cp,CD='Constrain',
                     kappa.calv=kappa.cal,kappa.cal.grid=kappa.cal.grid)


#test decision
p.h1.list=cbind(p.h1.info,p.h1.unif,p.h1.mix,p.h1.EB,p.h1.unif,p.h1.unif,
                p.h1.unif,p.h1.unif,p.h1.unif)


dec.list= cbind(p.h1.info>(1-alpha.b),
                p.h1.unif>(1-alpha.b),
                p.h1.mix>(1-alpha.b), 
                p.h1.EB>(1-alpha.b),
                p.h1.unif>(1-kappa.f),
                p.h1.unif>(1-kappa.f.c),
                p.h1.unif>(1-kappa.f1),
                p.h1.unif>(1-kappa.f2),
                p.h1.unif>(1-kappa.f3))


dat.dec=data.frame(yt=rep(xv[,1],9),yc=rep(xv[,2],9),p.h1=as.vector(p.h1.list),
                   dec=as.vector(dec.list),
                   Decision=rep(c('BD','FD','RMD-Unit','EBPowD','CD-Discard (4)','CD-Constraint',
                                  'CD-Discard (1)','CD-Discard (2)','CD-Discard (3)'),each=length(xv[,2])))


ytr=tapply(dat.dec$yt*dat.dec$dec,list(dat.dec$Decision,as.factor(dat.dec$yc)),
           function(x) 
           {x=unlist(x)
           if(any(x>0))
           {
             out=x[abs(x)>0][which.min(x[abs(x)>0])]
           }
           else
           {
             out=200 
           }
           #check monotnicity
           print(sum(diff(x)<0)==0)
           out
           })

#check monotonicity in yt-yc

monotonicity.check=tapply(dat.dec$yt*dat.dec$dec,list(dat.dec$Decision,as.factor(dat.dec$yc)),
           function(x) 
           {x=unlist(x)
           #check monotnicity
           sum(diff(x)<0)==0})

dat.plot.dec=data.frame(yc=as.numeric(rep(colnames(ytr),nrow(ytr))),
                        crit=as.vector(t(ytr)),
                        Decision=rep(rownames(ytr),each=ncol(ytr)))

dat.plot.probs=data.frame(ph1=p.h1.list[,'p.h1.info'],
                          ph1f=p.h1.list[,'p.h1.unif'],
                          yc=xv[,2],yt=xv[,1])

dat.plot.dec$Decision <- factor(dat.plot.dec$Decision, levels = c("CD-Constraint","CD-Discard (4)",
                                                                  "CD-Discard (3)","CD-Discard (2)",
                                                                  "CD-Discard (1)",
                                                                  "RMD-Unit","EBPowD","BD","FD"))


dec.sens = ggplot() + theme_light() + 
  geom_line(data=subset(dat.plot.dec, dat.plot.dec$Decision=="CD-Discard (1)" | 
                          dat.plot.dec$Decision=="CD-Discard (2)" |
                          dat.plot.dec$Decision=="CD-Discard (3)"|
                          dat.plot.dec$Decision=="CD-Discard (4)"|
                          dat.plot.dec$Decision=="FD"|
                          dat.plot.dec$Decision=="BD"), 
            aes(x=yc, y=(crit-yc), color=Decision),linewidth=0.8) + 
  scale_colour_manual(values=colors) +
  ylab(TeX("$y_{T}-y_{C}$")) + ylim(c(0,40))+
  xlab(TeX("$y_{C}$"))+
  theme(legend.title = element_blank(), legend.position = "top") + guides(colour = guide_legend(nrow = 2)) +
  theme(plot.title = element_text(size = 18))

combined <- ggplot() + theme_light() + 
  geom_contour_filled(aes(x = yc, y = (yt-yc), z =(1-ph1)*100),
                      data=dat.plot.probs,
                      breaks = (1-c(0,0.9,0.95,0.975,0.98,0.99,1))*100, alpha=0.2) + 
  scale_fill_viridis_d(name=TeX("$P^{\\pi}(H_0|y) \\%"))+
  geom_line(aes(x = yc, y = (crit-yc), color = Decision), linewidth=0.8, 
            data=subset(dat.plot.dec, dat.plot.dec$Decision!="CD-Discard (1)" & 
                          dat.plot.dec$Decision!="CD-Discard (2)" &
                          dat.plot.dec$Decision!="CD-Discard (3)"))+ ylim(c(-12,60))+
  scale_colour_manual(values=colors)+
  xlab(TeX("$y_{C}$"))+ 
  ylab(TeX("$y_{T}-y_{C}$")) +
  theme(legend.title = element_text(),
    legend.position = "top",
    plot.title = element_text(size = 18)
  ) +
    guides(
    colour = guide_legend(title = NULL, nrow = 2),
    fill = guide_legend(position = "bottom", title.theme = element_text(size = 8),
                        label.theme = element_text(size = 6),nrow = 1),
  )


#compute OCs, under vague prior, mixture prior, mixture prior + cap (CD constrain),
#informative prior, informative prior + cap and discard (CD Discard)

#effect size for power
delta=0.12
#conflict range
conflict=seq(-0.65,1-0.65,0.01)

ocs = compute.oc(
  outcome = outcome,
  nt = nt,
  nc = nc,
  n0c=n0c,
  delta = delta,
  conflict = conflict,
  priorPars = priorPars,
  alpha.low = alpha.low,
  alpha.up = alpha.up,
  alpha.b = alpha.b,
  w.mix=w.mix,
  cp = cp,
  cf = cf,
  ncores = 28
)


Decision = colnames(ocs$TIE)

out = data.frame(
  'Freq.type.I' = as.vector(ocs$TIE),
  'Freq.Power' = as.vector(ocs$Pow),
  'Type.I.cal' = as.vector(ocs$TIEcal),
  'Power.cal' = as.vector(ocs$Powcal),
  'Conflict' = rep(as.numeric(rownames(ocs$TIE)), ncol(ocs$TIE)),
  'Decision' = rep(Decision, each = nrow(ocs$TIE))
)

#find alphaUP and alphaLOW
# limq=qbeta(c(0.1, 0.9),mupc*n0c,(1-mupc)*n0c)
# f = approxfun(x=out$Conflict[out$Decision=='BD']+0.65,y=out$Freq.type.I[out$Decision=='BD'],rule=1,method='linear')
# f(limq)

theme = ggplot() + theme_light() +
  scale_colour_manual(values = colors) +
  xlab(expression(theta[C])) #+

out$Decision <- factor(out$Decision,
                       levels = c("CD-Constraint", "CD-Discard", "RMD-Unit", 
                                  "EBPowD", "BD", "FD"))


g1 <- theme +
  geom_segment(
    aes(
      x = 0,
      xend = 1,
      y = alpha.up,
      yend = alpha.up
    ),
    col = 'red',
    linetype = 'dashed',
    size = 0.25
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = 1,
      y = alpha.low,
      yend = alpha.low
    ),
    ,
    col = 'red',
    linetype = 'dashed',
    size = 0.3
  ) +
  geom_line(
    aes(Conflict+0.65, Freq.type.I, color = Decision),
    size = 0.5,
    data = out
  ) +
  ylab("TIE rate") +
  geom_vline(xintercept = 0.65,
             linetype = 'dashed',
             color = 'grey80') +
  ylim(c(0, 0.25)) + guides(fill = "none")

g2 <- theme + geom_line(
  aes(Conflict+0.65, Freq.Power, color = Decision),
  size = 0.5,
  data = out) + ylab(TeX('Power')) +
  geom_vline(xintercept = 0.65,
             linetype = 'dashed',
             color = 'grey80') +
  ylim(c(0, 1)) + theme(legend.position = "none")



dec.thr = (combined |
             dec.sens)+
  plot_layout(widths = c(1.2, 1))

Binomial_oc = (
  (g1 + guides(col = guide_legend(nrow = 1, byrow = T)) + theme(legend.title =
                                                                  element_blank())) | 
    g2 + guides(color = "none")
) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")


all = (dec.thr /
         Binomial_oc +
         plot_layout(heights = c(1, 1))
)

pdf(
  file = paste0("~/results/case_study_Binomial.pdf"),
  width = 9,
  height = 9,
  pointsize = 14,
  onefile = TRUE
)
print(all)
dev.off()

save(out, file = paste0("~/results/case_study_OCdat_Binomial.RData", sep = ""))
