##########################################
##       FitzHugh-Nagumo Example       ###
##########################################

# Obtain some pre-generated data

FhNdata

# Now we want to set up a basis expansion; this makes use of
# routines in the fda library

knots  = seq(0,20,0.5)
norder = 3
nbasis = length(knots) + norder - 2
range  = c(0,20)

bbasis = create.bspline.basis(range=range,nbasis=nbasis,
	norder=norder,breaks=knots)
	
	
# We'll start off by creating a smooth of each state variable to get initial
# values for the parameters.

bfdPar = fdPar(bbasis,lambda=1,int2Lfd(1))
DEfd = smooth.basis(FhNtimes,FhNdata,bfdPar,fdnames=list(NULL,c('V','R')) )$fd

coefs = DEfd$coefs
colnames(coefs) = FhNvarnames

#### Option 1: pre-defined set-up for squared error

lambda = 1000

# First just optimize the coefficients

Ires1	= Smooth.LS(make.fhn(),FhNdata,FhNtimes,FhNpars,coefs,bbasis,lambda,
  in.meth='nlminb')

# Let's have a look at this

IDEfd1 = fd(Ires1$coefs,bbasis)
plot(IDEfd1)
matplot(FhNtimes,FhNdata,add=TRUE)


# Now we can do the profiling

Ores1 = Profile.LS(make.fhn(),FhNdata,FhNtimes,FhNpars,coefs=coefs,
  basisvals=bbasis,lambda=lambda)

# And look at the result

ODEfd1 = fd(Ores1$coefs,bbasis)
plot(ODEfd1)
matplot(FhNtimes,FhNdata,add=TRUE)
  
# If we only coded the FitzHugh Nagumo funciton with no derivatives we would just
# have make.fhn()$fn, in which case we default to estimating the derivatives by
# finite differencing

Ires1a	= Smooth.LS(make.fhn()$fn,FhNdata,FhNtimes,FhNpars,coefs,bbasis,lambda,
  in.meth='nlminb')

# Now we can do the profiling

Ores1a = Profile.LS(make.fhn()$fn,FhNdata,FhNtimes,FhNpars,coefs=coefs,
  basisvals=bbasis,lambda=lambda)
  
  
  
#### Option 2:  set-up functions for proc and lik objects and then call
# optimizers

profile.obj = LS.setup(pars=FhNpars,fn=make.fhn(),lambda=lambda,times=FhNtimes,
      coefs=coefs,basisvals=bbasis)
lik = profile.obj$lik
proc= profile.obj$proc

# Now we can get initial parameter estimates from "gradient matching"

pres = ParsMatchOpt(FhNpars,coefs,proc)
npars = pres$pars

# Smoothing can be done more specifically with 'inneropt'

Ires2 = inneropt(FhNdata,times=FhNtimes,npars,coefs,lik,proc)

# And we can also do profiling

Ores2 = outeropt(FhNdata,FhNtimes,npars,coefs,lik,proc)

#### Option 3:  set everything up manually

## lik object

lik2 = make.SSElik()                      # we are using squared error
lik2$bvals = eval.basis(FhNtimes,bbasis)  # values of the basis at observation times
lik2$more = make.id()                     # identity transformation
lik2$more$weights = c(1,0)                # only use data for V

## proc object

qpts = knots[1:(length(knots)-1)]+diff(knots)/2  # Collocation points at midpoints
                                                 # between knots

proc2 = make.SSEproc()                    # we are using squared error
proc2$bvals = list(bvals = eval.basis(qpts,bbasis),    # values and derivative of
                  dbvals = eval.basis(qpts,bbasis,1))  # basis at collocation points
proc2$more = make.fhn()                   # FitzHugh-Nagumo right hand side
proc2$more$names = FhNvarnames            # State variable names
proc2$more$parnames = FhNparnames         # Parameter names
proc2$more$qpts = qpts                    # Collocation points
proc2$more$weights  = c(1000,1000)        # Weight relative to observations of
                                          # each collocation point.
                                          
# We can also specify optimization methods

in.meth = 'nlminb'            # Inner Optimization
control.in = list(rel.tol=1e-12,iter.max=1000,eval.max=2000,trace=0)  # Optimization control

out.meth = 'optim'            # Outer Optimization
control.out = list(trace=6,maxit=100,reltol=1e-8,meth='BFGS') # Optimization control
                                          
                                          
# Now we will also remove the 'R' part of the data (ie, assume we did not measure
# it) and set the corresponding coefficients to zero

new.FhNdata = FhNdata
new.FhNdata[,2] = NA

new.coefs = coefs
new.coefs[,2] = 0

# We can now fix the smooth or 'V' and try and choose 'R' to make the differential
# equation match as well as possible.

fres = FitMatchOpt(coefs=new.coefs,which=2,pars=FhNpars,proc2)
new.coefs2 = fres$coefs

# And we can call the same inner optimization as above, this time with our
# own control parameters and method

Ires3 = inneropt(new.FhNdata,FhNtimes,FhNpars,new.coefs2,lik2,proc2,
  in.meth=in.meth,control.in=control.in)

# And we can also do profiling with specified controls

Ores3 = outeropt(new.FhNdata,FhNtimes,FhNpars,coefs=new.coefs2,lik=lik2,proc=proc2,
  in.meth=in.meth,out.meth=out.meth,control.in=control.in,control.out=control.out)

# We can also use finite differencing to calculate derivatives of the right hand side of the
# FitzHugh-Nagumo equations, this is defined by modifying the proc object.

proc2a = make.SSEproc()                    # we are using squared error
proc2a$bvals = list(bvals = eval.basis(qpts,bbasis),    # values and derivative of
                  dbvals = eval.basis(qpts,bbasis,1))  # basis at collocation points
proc2a$more = make.findif.ode()                   # FitzHugh-Nagumo right hand side
proc2a$more$names = FhNvarnames            # State variable names
proc2a$more$parnames = FhNparnames         # Parameter names
proc2a$more$qpts = qpts                    # Collocation points
proc2a$more$weights  = c(1000,1000)        # Weight relative to observations of
                                           # each collocation point.
proc2a$more$more = list(fn = make.fhn()$fn, eps=1e-6) # Tell findif that the function
                                                      # to difference is fhn and the
                                                      # difference stepsize


Ires3a = inneropt(new.FhNdata,FhNtimes,FhNpars,new.coefs2,lik2,proc2a,
  in.meth=in.meth,control.in=control.in)

# And we can also do profiling with specified controls

Ores3a = outeropt(new.FhNdata,FhNtimes,FhNpars,coefs=new.coefs2,lik=lik2,proc=proc2a,
  in.meth=in.meth,out.meth=out.meth,control.in=control.in,control.out=control.out)
