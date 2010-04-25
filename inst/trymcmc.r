source('Devel/Colloc.MCMC.R')

walk.var = diag(c(0.005,0.005,0.05))
cscale = 1

prior = function(pars){ return( sum(dnorm(pars/c(1,1,10),log=TRUE))) }
nstep = 50


hmm = Colloc.MCMC(FhNtimes,FhNdata,FhNpars,Ires2$coefs,srklik,srkproc,prior,walk.var,cscale,nstep,in.meth='SplineEst',control.in=control.in)
