require(RBesT)
require(parallel)
require(VGAM)

#true TIE rate at thc=eval
tIe=function(alpha.ev, At=0, Ac, mupc, mupt=0, sigma, nt, nc, eval, th0=0)
{1-pnorm((th0+ mupc*(Ac/(Ac+1))-mupt*(At/(At+1))-eval*((1/(At+1))-(1/(Ac+1))) +  
            qnorm(1-alpha.ev)*sqrt(((sigma^2/nc)/(1+Ac)) + ((sigma^2/nt)/(1+At))))/
           sqrt(((sigma^2/nc)/(1+Ac)^2) + ((sigma^2/nt)/(1+At)^2)))}

#find delta=(thc-mupc) corresponding to alpha.target (true) TIE rate at thc
# deltafun=function(alpha.b, alpha.target, Ac, mupc,  sigma, nt, nc, th0=0) 
# {-((Ac+1)/(Ac))*(qnorm(1-alpha.target)*sqrt((sigma^2/nc)/(1+Ac)^2 + sigma^2/nt) - qnorm(1-alpha.b)*
#                    sqrt((sigma^2/nc)/(1+Ac) + sigma^2/nt) - th0)} 

#threshold to obtain the BD under the frequentist test
kappa.cal.f=function(yc, alpha.b, At, Ac, mupc, mupt, sigma, nt, nc, th0=th0) 
{1-pnorm((th0*(At/(At+1)) + mupc*(Ac/(Ac+1))-mupt*(At/(At+1))-yc*((1/(At+1))-(1/(Ac+1))) + 
            qnorm(1-alpha.b)*sqrt(((sigma^2/nc)/(1+Ac)) + ((sigma^2/nt)/(1+At))))/
           sqrt((sigma^2/nc)/((At+1)^2) + (sigma^2/nt)/((At+1)^2)))} 

#data value triggering change of rule at delta (=eval)
# ycfun=function(alpha.b, Ac, mupc,  sigma, nt, nc, eval, th0=0) 
# {((Ac+1)/(-Ac))*((th0 + eval*(Ac/(Ac+1)) + qnorm(1-alpha.b)*
#                     sqrt((sigma^2/nc)/(1+Ac) + sigma^2/nt))*sqrt(1/((sigma^2/nc)/((Ac+1)^2) + (sigma^2/nt)))*sqrt(sigma^2/nc + sigma^2/nt)+
#                    - mupc*(Ac/(Ac+1))-qnorm(1-alpha.b)*sqrt((sigma^2/nc)/(Ac+1) + (sigma^2/nt)))} 

#data value triggering change of rule at alpha.c 
# ycfun.alpha=function(alpha.b, alpha.c, Ac, mupc,  sigma, nt, nc) 
# {mupc-((Ac+1)/(Ac))*(qnorm(1-alpha.c)*sqrt(sigma^2/nc + sigma^2/nt)-qnorm(1-alpha.b)*sqrt((sigma^2/nc)/(1+Ac) + sigma^2/nt))} 

#threshold to be used for the Bayes test giving rule equivalent to that of threshold alpha.k under the FD
kappa.bayes.general=function(alpha.k, eval, Ac, mupc, sigma, nt, nc, th0=0)
{1-pnorm(((eval-mupc)*(Ac/(Ac+1)) +
            qnorm(1-alpha.k)*sqrt(sigma^2/nc + sigma^2/nt))/
           sqrt(((sigma^2/nc)/(1+Ac)) + sigma^2/nt))}


kappa.discount=function(alpha.k, alpha.k2, eval, muc, cf, cp) 
{ 
  w=pmin((abs(eval-muc)/cf)^cp,1)
  (alpha.k)^w*(alpha.k2)^(1-w)
}

p.h1.precomp <-function(pow,nc,nt,n0c,a0t=0.5,a0c=0.5,b0t=0.5,b0c=0.5,th0=0)
{
  sapply(0:(pow*n0c+nc),function(xS)
  {
    sapply(0:nt, function(xT)
    {
      post.treat=mixbeta(c(1, a0t + xT, b0t + nt - xT))
      post.control=mixbeta(c(1, a0c + xS, b0c + pow*n0c + nc -xS))
      out=try(pmixdiff(post.treat,post.control,th0,lower.tail = FALSE), silent=TRUE)
      if(inherits(out,"try-error"))
      {
        print("use freq. unconditional test")
        #out=1-fisher.test(matrix(c(xS,pow*n0c+nc-xS,xT,nt-xT),2,2),alternative = "less")$p.value
        out=1-pnorm((xS*(nt-xT)-xT*(pow*n0c+nc-xS))*sqrt((pow*n0c+nc+nt)/((pow*n0c+nc)*nt*(xS+xT)*(pow*n0c+nc+nt-xS-xT))))
      }
      as.numeric(out)
    })
  })
}



cThreshold=function(yc,yt,nc,nt,priorPars,outcome='normal',
                    mc.results=NULL,g=0.05,sigma=1,alpha.b=0.025,
                    alpha.up=1,alpha.low=0,a0t=0.5,a0c=0.5,b0t=0.5,b0c=0.5, th0=0, 
                    pH1array=NULL, cf=NULL,cp=NULL, CD="Constrain",
                    kappa.calv=NULL, kappa.cal.grid=NULL)
{
  if(outcome=='normal')
  {
    mupc=priorPars[1,1]
    mupt=priorPars[1,2]
    sigmapc=priorPars[2,1]
    sigmapt=priorPars[2,2]
    
    Ac=sigma^2/(nc*sigmapc^2)
    At=sigma^2/(nt*sigmapt^2)
    eval=yc 
    
    if(is.null(kappa.calv))
    {
    kappa.cal=kappa.cal.f(yc=eval, alpha.b=alpha.b, At=At, Ac=Ac, mupc=mupc, mupt=mupt, sigma=sigma, nt=nt, nc=nc, th0=th0) 
    }
    if(!is.null(kappa.calv))
    {
    kappa.calf = approxfun(x=kappa.cal.grid,y=kappa.calv,rule=1,method='linear')
    kappa.cal = kappa.calf(eval)
    }
    
    kappa= pmax(pmin(kappa.cal,alpha.up),alpha.low)
    
    if (CD=='Discard')
    {
      kappa=pmax(pmin(kappa.discount(alpha.k=alpha.b, alpha.k2=kappa.cal, eval=eval, muc=mupc, cf=cf, cp=cp),alpha.up),alpha.low)#, 
    }
    
  }
  
  if(outcome=='binomial')
  {
    ae=priorPars[1,1]
    at=priorPars[1,2]
    be=priorPars[2,1]
    bt=priorPars[2,2]
    n0c=ae+be
    n0=a0t+a0c
    
    xv=expand.grid(yt,yc)
    
    yt1=xv[,1]
    yc1=xv[,2]
    
    
    thc.e=yc/nc #assumed control arm prop
    wmat=sapply(thc.e,function(z) exp(dbinom(yc1,nc,z,log=TRUE)+dbinom(yt1,nt,z,log=TRUE))) #weights un H0
    
    p.h1.unif=as.vector(pH1array[[1]][(1+yt),(1+yc)]) #pmixdiff(post.treat,post.control,0,lower.tail = FALSE) #approximated if at is non-integer
    p.h1.info=as.vector(pH1array[[2]][(1+yt+at),(1+yc+ae)]) #pmixdiff(post.treat,post.control.info,0,lower.tail = FALSE) 
    
    #informative prior on control and tr.
    tIe.info=colSums(wmat*(p.h1.info>(1-alpha.b)))
    tIe.vague=colSums(wmat*(p.h1.unif>(1-alpha.b))) 
    
    #difference between noninformative and informative prior pr of the null
    pa=(p.h1.unif-p.h1.info)
  
    kappa=as.vector(sapply(1:length(thc.e), function(z)
    {
      out=matrix(pmax(pmin(alpha.b-pa,alpha.up+1e-10),alpha.low-1e-10),nt+1,nc+1)[,z]
      out
    }))
    
    if (CD=='Discard')
    {
      alpha.k2=pmax(pmin(alpha.b-pa,1-1e-10),1e-10)
      rhohat=(ae+yc1)/(nc+n0c)
      kappa=pmax(pmin(kappa.discount(alpha.k=alpha.b, alpha.k2=alpha.k2, eval=(yc1/nc), muc=(ae/n0c), cf=cf*sqrt(rhohat*(1-rhohat)*(1/nc+1/n0c)), cp=cp),alpha.up+1e-10),alpha.low-1e-10)#, 
    }
  }
  
  kappa
}

#adapted from studyPrior package
power_par = function(dat,n,sigma=1,priorPars,outcome,bound=FALSE, b0=0.5,a0=0.5) {
  if(outcome=="normal")
  {
    sD = sigma^2/n
    d = priorPars[2]^2 / (pmax((dat - priorPars[1]) ^ 2, sD + priorPars[2]^2) - sD)
    if(bound) d = pmin(d,priorPars[2]^2/(sigma^2/n))
  }
  if(outcome=="binomial")
  {
    optfn=function(datC)
    {
      lik.d = function(d) VGAM::dbetabinom.ab(x=datC, size=n, d * (priorPars[1]) + a0,  d *(priorPars[2]) + b0) #uniform historical prior
      
      opd = BB::spg(par = .005,
                    fn = lik.d,
                    lower=0,
                    upper=1,
                    control=list(maximize=TRUE,
                                 trace=FALSE))
      
      if(opd$convergence!=0) print(opd)
      
      ds = opd$par
      
      return(ds)
    }
    d=sapply(dat,optfn)
    if(bound) d=pmin(d,n/(priorPars[1]+priorPars[2]))
  }
  return(d)
}


pH1=function(yc,yt,nc,nt,priorPars,w.mix=0.3,outcome='normal',sigma=1,prior='info',bound=FALSE,pH1array=NULL,
             a0t=0.5,a0c=0.5,b0t=0.5,b0c=0.5,n0=1,th0=0,d1=NULL,delta.v=NULL,w=0.5,mix.control=NULL,mix.treat=NULL)
{
  if(outcome=="normal")
  {
    if(!is.null(priorPars))
    {
    mupc=priorPars[1,1]
    mupt=priorPars[1,2]
    sigmapc=priorPars[2,1]
    sigmapt=priorPars[2,2]
    }
    
    if(prior=='mix')
    {
      if(is.null(mix.control))
      {
      if(sigmapc<10^3)
      {mix.control=mixnorm(info=c((1-w.mix),mupc,sigmapc), unit=c(w.mix,mupc,sigma), sigma=sigma)}
      else
      {mix.control=mixnorm(info=c(1,mupc,10^3), sigma=sigma)}
      }
      
      if(is.null(mix.treat))
      {
      mix.treat=mixnorm(info=c(1,mupt,10^3), sigma=sigma)
      }
      
      p.h1=sapply(1:length(yc), function(z)
      {
        post.control.mix=suppressMessages(postmix(mix.control,m=yc[z],n=nc))
        post.treat=suppressMessages(postmix(mix.treat,m=yt[z],n=nt))
        pmixdiff(post.treat,post.control.mix,th0,lower.tail = FALSE)
      })
    }
    else
    {
      if(prior=='EB')
      {
        if(sigmapc<10^3)
        {  
          dc =  power_par(dat=yc, n=nc, sigma=sigma, priorPars=c(mupc,sigmapc), outcome='normal', bound=bound)  
          Ac=dc*sigma^2/(nc*sigmapc^2)}
        else {Ac=0}
      }
      else
      {
        Ac=sigma^2/(nc*sigmapc^2)
      }
      
      At=0
      
      pthr=ifelse(rep(sigmapt,length(yc))>0,
                  ((mupt*At + yt)/(At+1) - (mupc*Ac + yc)/(Ac+1))/sqrt(((sigma^2/nc)/(1+Ac)) + ((sigma^2/nt)/(1+At))),
                  (mupt - (mupc*Ac + yc)/(Ac+1))/sqrt(((sigma^2/nc)/(1+Ac))))
      
      p.h1=pnorm(th0,pthr,1,lower.tail=FALSE)
    }
  }
  
  if(outcome=="binomial")
  {
    ae=priorPars[1,1]
    at=priorPars[1,2]
    be=priorPars[2,1]
    bt=priorPars[2,2]
    n0c=ae+be
    
    xv=expand.grid(yt,yc)
    
    yt1=xv[,1]
    yc1=xv[,2]
    
    if(prior=='vague')
    {
      p.h1=as.vector(pH1array[[1]][1+yt,1+yc]) #pmixdiff(post.treat,post.control,0,lower.tail = FALSE) #approximated if at is non-integer
    }
    
    if(prior=='info')
    {
      p.h1=as.vector(pH1array[[2]][(1+yt+at),(1+yc+ae)]) #pmixdiff(post.treat,post.control.info,0,lower.tail = FALSE) 
    }
    
    if(prior=='mix')
    {
      if(is.null(mix.control))
      {
      if(n0c>0)
      {  
        mix.contr=mixbeta(info=c((1-w.mix),ae+a0c,be+b0c), unit=c(w.mix,a0c,b0c))
      }
      else
      {
        mix.contr=mixbeta(unit=c(1,a0c,b0c))
      }
      }
      
      if(is.null(mix.treat))
      {
      mix.treat=mixbeta(unit=c(1,a0t,b0t))
      }
      
      p.h1=sapply(1:length(yc1), function(z)
      {
        post.control.mix=suppressMessages(postmix(mix.contr,r=yc1[z],n=nc))
        post.treat=suppressMessages(postmix(mix.treat,r=yt1[z],n=nt))
        pmixdiff(post.treat,post.control.mix,th0,lower.tail = FALSE)})
    }
    
    if(prior=='EB')
    {
      if(n0c>0)
      {  
        dc =  power_par(dat=yc1, n=nc, priorPars=c(ae,be), b0=0.5, a0=0.5,outcome='binomial', bound=bound)  
      }
      else {dc=0}
      
      post.a.c =  dc * ae + a0c + yc1
      post.b.c =  dc * be + b0c + (nc-yc1)
      post.a.t =  a0t + yt1
      post.b.t =  b0t + (nt-yt1)
      
      p.h1=sapply(1:length(yc1), function(z)
      {
        pmixdiff(mixbeta(c(1,post.a.t[z],post.b.t[z])),
                 mixbeta(c(1,post.a.c[z],post.b.c[z])),th0,lower.tail = FALSE)
      })
    }
    
  }
  
  p.h1
}


compute.oc=function( 
    nt, # treat arm ss
    nc, # control arm ss
    n0c=NULL, #control prior EHSS
    delta, #treatm effect for power evaluation (ocs='fixed')
    conflict, #vector of bias for control arm
    priorPars, #2x2 matrix of prior parameters (row1: mean/num of succ, row2: sd/num of failures)
    sigma=1, #sd data (normal outcomes only)
    a0t=0.5, #default noninformative prior num of successes for binomial - unconditional approx. test
    a0c=0.5,
    b0t=0.5,
    b0c=0.5,
    alpha.low=0, #min TIE
    alpha.up=1, #max TIE
    alpha.b = 0.025, #significance level of FD
    w.mix=0.3, #mixture prior weight
    mix.control=NULL, #pre-specified (robust) mixture for the control arm
    mix.treat=NULL, #pre-specified (robust) mixture for the treatment arm
    outcome='normal',
    kappa.cal=NULL, #pre-specified test decision threshold to be applied to the frequentist test to obtain test under chosen prior 
    kappa.cal.grid=NULL, #pre-specified grid test decision threshold to be applied to the frequentist test to obtain test under chosen prior (for interpolation)
    bound=FALSE, #bound EB informativeness
    nmc=10^4, #num of MC samples (only for outcome='normal')
    seed=1234,
    th0=0,
    cf=NULL,
    cp=NULL,
    ncores=1)
{
  if(outcome=="binomial") 
  {
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
                       th0=th0,cf=cf,cp=cp,CD='Discard',
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
    p.h1.list=cbind(p.h1.info,p.h1.unif,p.h1.mix,p.h1.EB)
    
    
    dec.list= cbind(p.h1.info>(1-alpha.b),
                    p.h1.unif>(1-alpha.b),
                    p.h1.mix>(1-alpha.b), 
                    p.h1.EB>(1-alpha.b),
                    p.h1.unif>(1-kappa.f),
                    p.h1.unif>(1-kappa.f.c))
  }
  
  
  ocps.c=mclapply(conflict, function(c1)
  {
    print(c1)
    
    if(outcome=='binomial')
    {
      thc=(priorPars[1,1]/n0c)+c1
    }
    if(outcome=='normal')
    {
      thc=priorPars[1,1]+c1
    }
    tht.tie=thc
    tht.pow=thc+delta
    #print(c(tht.tie,tht.pow))
    
    if(outcome=="binomial")
    {
      dec.list.f=lapply(c(tht.tie,tht.pow), function(tht)
      {
        wmatT=exp(dbinom(yc1,nc,thc,log=TRUE)+dbinom(yt1,nt,tht,log=TRUE)) #weights for integration
        
        list(dec.list,wmatT,p.h1.list)
      })
      
    }
    if(outcome=="normal")
    {
      set.seed(seed)
      
      yc=rnorm(nmc,thc,sigma/sqrt(nc))
      dec.list.f=lapply(c(tht.tie,tht.pow), function(tht)
      {
        yt=rnorm(nmc,tht,sigma/sqrt(nt))
        
        #Weight vector, 1 for monte carlo
        wmatT=rep(1/nmc,nmc)
        
        priorPars.v=priorPars
        priorPars.v[2,]=c(Inf,Inf)
        p.h1.unif=pH1(yc=yc,yt=yt,nc=nc,nt=nt,sigma=sigma,
                      priorPars=priorPars.v,w.mix=w.mix,
                      outcome='normal',th0=th0)
        
        p.h1.info=pH1(yc=yc,yt=yt,nc=nc,nt=nt,sigma=sigma,
                      priorPars=priorPars,w.mix=w.mix,
                      outcome='normal',th0=th0)
        
        p.h1.mix=pH1(prior='mix',yc=yc,yt=yt,nc=nc,nt=nt,sigma=sigma,
                     priorPars=priorPars,w.mix=w.mix,
                     mix.treat=mix.treat,mix.control=mix.control,
                     outcome='normal',th0=th0)
        
        p.h1.EB=pH1(prior='EB',yc=yc,yt=yt,nc=nc,nt=nt,sigma=sigma,
                    priorPars=priorPars,w.mix=w.mix,
                    outcome='normal',bound=bound,th0=th0)
        
        
        kappa.f=cThreshold(yc=yc,yt=yt,nc=nc,nt=nt,sigma=sigma,
                           priorPars=priorPars,
                           outcome='normal', mc.results=mc.results, 
                           alpha.b=alpha.b, alpha.up=alpha.up, 
                           alpha.low=alpha.low, th0=th0,
                           cf=cf,cp=cp,CD='Discard',
                           kappa.calv=NULL,kappa.cal.grid=NULL)
        
        
        kappa.f.c=cThreshold(yc=yc,yt=yt,nc=nc,nt=nt,sigma=sigma,
                             priorPars=priorPars,
                             outcome='normal', mc.results=mc.results, 
                             alpha.b=alpha.b, alpha.up=alpha.up, 
                             alpha.low=alpha.low, th0=th0,
                             cf=cf,cp=cp,CD='Constrain',
                             kappa.calv=kappa.cal,kappa.cal.grid=kappa.cal.grid)
        
        #test decision
        p.h1.list=cbind(p.h1.info,p.h1.unif,p.h1.mix,p.h1.EB)
        
        dec.list = cbind(p.h1.info>(1-alpha.b),
                         p.h1.unif>(1-alpha.b),
                         p.h1.mix>(1-alpha.b), 
                         p.h1.EB>(1-alpha.b),
                         p.h1.unif>(1-kappa.f),
                         p.h1.unif>(1-kappa.f.c))
        
        list(dec.list,wmatT,p.h1.list)
      })
      
    }
    
    pow=1-(dec.list.f[[2]][[2]] %*% (1-dec.list.f[[2]][[1]]))
    tie=dec.list.f[[1]][[2]] %*% (dec.list.f[[1]][[1]])
    
    out.c=list(cbind(t(tie),t(pow)),dec.list.f)
    out.c
  },mc.cores=ncores)
  
  ############TIE rate and Power
  pow=do.call(rbind,lapply(ocps.c,function(x) x[[1]][,2]))
  tie=do.call(rbind,lapply(ocps.c,function(x) x[[1]][,1]))
  
  colnames(tie)=colnames(pow)=c('BD','FD','RMD-Unit','EBPowD','CD-Discard','CD-Constraint')
  rownames(tie)=rownames(pow)=conflict
  
  ############recalibrated TIE rate and Power
  
  #extract post prob for each conflict value (for binomial outcomes equal)
  ph1.mat=do.call(cbind,lapply(1:length(ocps.c),function(x) ocps.c[[x]][[2]][[1]][[3]]))
  ph1.mat.pow=do.call(cbind,lapply(1:length(ocps.c),function(x) ocps.c[[x]][[2]][[2]][[3]]))

  
  #extract integration weights for each conflict value (for normal outcomes equal)
  w.mat= do.call(cbind,lapply(1:length(ocps.c),function(x) ocps.c[[x]][[2]][[1]][[2]]))
  w.mat.pow= do.call(cbind,lapply(1:length(ocps.c),function(x) ocps.c[[x]][[2]][[2]][[2]]))
  
  #select threshold giving max TIE rate as the max of CD
  alpha.recal=max(tie[,'CD-Discard'])
  
  #function to select the threshold for target TIE rate
  tie.sim=function(thr,ph1vals,wvals)
  {
    abs(wvals%*%(ph1vals>(1-thr))-alpha.recal)
  }
  
  #find threshold for each conflict scenario
  thrC=sapply(1:ncol(ph1.mat), function(t) optim(0.025, function(thr) tie.sim(thr, ph1vals=ph1.mat[,t],wvals=w.mat[,ceiling(t/4)]),
                                                 lower=0, upper=1, method = "L-BFGS-B")$par)
  #take the minimum for each approach across conflict
  thrMB=tapply(thrC,rep(c('info','unif','mix','EB'),length(ocps.c)),min)[c(2,4,3,1)]
  
  #compute calibrated TIE rate and power
  dec.mat=t(apply(ph1.mat,1,function(x) {x>1-rep(thrMB,length(ocps.c))}))
  t1e.cal=matrix(sapply(1:ncol(dec.mat), function(t) w.mat[,ceiling(t/4)]%*%dec.mat[,t]),length(ocps.c),4, byrow = TRUE)
  
  dec.mat.pow=t(apply(ph1.mat.pow,1,function(x) {x>1-rep(thrMB,length(ocps.c))}))
  pow.cal=matrix(sapply(1:ncol(dec.mat.pow), function(t) w.mat.pow[,ceiling(t/4)]%*%dec.mat.pow[,t]),length(ocps.c),4, byrow = TRUE)
  
  #recalibrate CD-Constraint (if needed)
  
  if (abs(max(tie[,'CD-Discard'])-max(tie[,'CD-Constraint']))>0.0001)
  {
  ocps.c.rec=mclapply(conflict, function(c1)
  {
    print(c1)
    
    if(outcome=='binomial')
    {
      thc=(priorPars[1,1]/n0c)+c1
    }
    if(outcome=='normal')
    {
      thc=priorPars[1,1]+c1
    }
    tht.tie=thc
    tht.pow=thc+delta
    #print(c(tht.tie,tht.pow))
    
    if(outcome=="binomial")
    {
      kappa.f.c=cThreshold(yc=yc,yt=yt,nc=nc,nt=nt,
                           priorPars=priorPars,
                           outcome='binomial',
                           alpha.b=alpha.b, alpha.up=alpha.recal, pH1array=pH1array,
                           alpha.low=alpha.low, a0t=a0t, a0c=a0c,b0t=b0t,b0c=b0c,
                           th0=th0,
                           cf=cf,cp=cp,CD='Constrain',
                           kappa.calv=kappa.cal,kappa.cal.grid=kappa.cal.grid)
      
      
      #test decision
      p.h1.list=cbind(p.h1.unif)
      
      dec.list= cbind(p.h1.unif>(1-kappa.f.c))
      
      dec.list.f=lapply(c(tht.tie,tht.pow), function(tht)
      {
        wmatT=exp(dbinom(yc1,nc,thc,log=TRUE)+dbinom(yt1,nt,tht,log=TRUE)) #weights for integration
        
        list(dec.list,wmatT,p.h1.list)
      })
      
    }
    if(outcome=="normal")
    {
      set.seed(seed)
      
      yc=rnorm(nmc,thc,sigma/sqrt(nc))
      dec.list.f=lapply(c(tht.tie,tht.pow), function(tht)
      {
        yt=rnorm(nmc,tht,sigma/sqrt(nt))
        
        #Weight vector, 1 for monte carlo
        wmatT=rep(1/nmc,nmc)
        
        priorPars.v=priorPars
        priorPars.v[2,]=c(Inf,Inf)
        p.h1.unif=pH1(yc=yc,yt=yt,nc=nc,nt=nt,sigma=sigma,
                      priorPars=priorPars.v,w.mix=w.mix,
                      outcome='normal',th0=th0)
        
        kappa.f.c=cThreshold(yc=yc,yt=yt,nc=nc,nt=nt,sigma=sigma,
                             priorPars=priorPars,
                             outcome='normal', mc.results=mc.results, 
                             alpha.b=alpha.b, alpha.up=alpha.recal, 
                             alpha.low=alpha.low, th0=th0,
                             cf=cf,cp=cp,CD='Constrain',
                             kappa.calv=kappa.cal,kappa.cal.grid=kappa.cal.grid)
        
        #test decision
        p.h1.list=cbind(p.h1.unif)
        dec.list = cbind(p.h1.unif>(1-kappa.f.c))
        
        list(dec.list,wmatT,p.h1.list)
      })
      
    }
    
    pow=1-(dec.list.f[[2]][[2]] %*% (1-dec.list.f[[2]][[1]]))
    tie=dec.list.f[[1]][[2]] %*% (dec.list.f[[1]][[1]])
    
    out.c=list(cbind(t(tie),t(pow)),dec.list.f)
    out.c
  },mc.cores=ncores)

  pow.CDC=do.call(rbind,lapply(ocps.c.rec,function(x) x[[1]][,2]))
  tie.CDC=do.call(rbind,lapply(ocps.c.rec,function(x) x[[1]][,1]))
  
  #
  t1e.cal=cbind(t1e.cal,tie[,'CD-Discard'],tie.CDC)
  pow.cal=cbind(pow.cal,pow[,'CD-Discard'],pow.CDC)
  }
  else
{
  t1e.cal=cbind(t1e.cal,tie[,'CD-Discard'],tie[,'CD-Constraint'])
  pow.cal=cbind(pow.cal,pow[,'CD-Discard'],pow[,'CD-Constraint'])
}
  colnames(t1e.cal)=colnames(pow.cal)=c('BD','FD','RMD-Unit','EBPowD','CD-Discard','CD-Constraint')
  rownames(t1e.cal)=rownames(pow.cal)=conflict
  
  ####################
  
  list(TIE=tie,Pow=pow,TIEcal=t1e.cal,Powcal=pow.cal)
}

