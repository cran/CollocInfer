source('Devel/Colloc.MCMC.R')

walk.var = diag(c(0.005,0.005,0.05))
cscale = 1

prior = function(pars){ return( sum(dnorm(pars/c(1,1,10),log=TRUE))) }
nstep = 20


hmm = Colloc.MCMC(FhNtimes,FhNdata,FhNpars,coefs,lik,proc,prior,walk.var,cscale,nstep,in.meth='SplineEst',control.in=control.in)
