# An example file fitting the FitzHugh-Nagumo equations to data in the
# new R Profiling Code. This will eventually be interfaced with the new R 
# "Partially Observed Markov Process (pomp)" class of objects. 

#  Last modifed 2 Feb 2010 by Jim

library("CollocInfer")

#########################################################
####   Data generation from a solution of the ODE #######
#########################################################

#  I wanted fewer and noisier observations

#  41 time values from 0 to 20 in steps of 0.5

t = seq(0,20,0.5)

n = length(t)

#  parameter values

pars = c(0.2,0.2,3)
names(pars) = c('a','b','c')

#  initial values and variable labels

x0 = c(-1,1)
names(x0) = c('V','R')

#  define the functions

fhn = make.fhn()

#  approximate the solution for these values over about three cycles

solution = lsoda(x0,times=t,func=fhn$fn.ode,pars)

#  extract the function values at time points

y = solution[,2:3]

#  add noise to the function values to define data

sigerr = 0.2
data = y + sigerr*array(rnorm(2*n),dim(y))

matplot(t, y, "l")
matpoints(t, data, pch="*")

###############################
####  Basis object bbasis  ####
###############################

knots  = seq(0,20,0.2)
norder = 4
nbasis = length(knots) + norder - 2
range  = c(0,20)

bbasis = create.bspline.basis(range=range,nbasis=nbasis,
	                            norder=norder,breaks=knots)

########################################
### Initial values for coefficients  ###
########################################

fd.data = array(data,c(dim(data)[1],1,dim(data)[2]))

bfdPar = fdPar(bbasis, 2, 1e-2)

DEfd = smooth.basis(t,fd.data,bfdPar,fdnames=list(NULL,NULL,c('V','R')) )$fd

par(mfcol=c(2,1))
plot(DEfd)

#  DEfd = smooth.basis(t,fd.data[,1],bbasis)$fd

coefs = DEfd$coefs[,1,]

######################################################################
### Set up likelihood and process functions for squared error loss ###
######################################################################

#  set up functions to evaluate values of ODE right side and values
#  of its derivatives

profile.obj = LS.setup(pars=pars, fn=fhn, lambda=10000, times=t,
                        fd.obj=DEfd)

lik  = profile.obj$lik

proc = profile.obj$proc

#############################################
####  set optimization control parameters ###
#############################################

#  control

control=list()                 
control$trace  = 1

#  control.in

control.in = control

#  control.out

control.out = control
control.out$trace = 2

###############################
#### Parameter Optimization ###
###############################

# Perturbed parameter starting values

spars = c(0.2,0.2,2)          
names(spars) = names(pars)

#  set smoothing level (value to be used for both variables)

lambda = 10000

####################################################
### Estimate parameters using squared error loss ###
####################################################

#  This doesn't work.

Ires1	= Smooth.LS(fhn,data,t,pars=spars,coefs,bbasis,lambda=lambda,
                   in.meth='nlminb',control.in=control.in,names=names(x0))
  
Ores1 = Profile.LS(fhn,data=data,times=t,pars=spars,coefs=coefs,
                    basisvals=bbasis,lambda=lambda,names=names(x0),
                    in.meth='nlminb',out.meth='nls',
	                 control.in=control.in,control.out=control.out)
	
####################################################
###         SSE with ProfileErr                  ###
####################################################

#  This does.

Ires2 = inneropt(data,times=t,pars,coefs,lik,proc,in.meth='nlminb',control.in)

Ores2 = outeropt(data=data,times=t,pars=pars,coefs=coefs,lik=lik,proc=proc,
                 in.meth="nlminb",out.meth="nls",
                 control.in=control.in,control.out=control.out)


####################################################
###          Lets look at the result             ###
####################################################

Ores = Ores2

DEfd = fd(Ores$coefs,bbasis)   # Data and reconstructed trajectory

par(mfrow=c(2,1))
plotfit.fd(data,t,DEfd)
 
#  Ores2 does not return lik and proc, 
#  so Ores$lik and Ores$proc changed to lik and proc, resp.

traj = as.matrix(proc$bvals$bvals %*% Ores$coefs)    # Look at how well the
colnames(traj) = proc$more$names                     # derivative of the
dtraj = as.matrix(proc$bvals$dbvals %*% Ores$coefs)  # trajectory fits the 
ftraj = proc$more$fn(t,traj,Ores$pars)               # right hand side. 

matplot(dtraj,type='l',col=2)
matplot(ftraj,type='l',col=4,add=TRUE) 
 
#  This doesn't work

Profile.covariance(pars=Ores$pars,times=t,data=data,coefs=Ores$coefs,lik=lik,proc=proc)
  
  
	
