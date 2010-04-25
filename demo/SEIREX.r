##########
# This demonstration illustrates the use of log transformations with
# CollocInfer. We use the SEIR equations with a seasonally varying
# infection rate for this purpose.
#
# The equations are given as
#
# Sdot = mu - beta(t)*S*(I+i) - nu*S        (Susceptibles)
# Edot = beta(t)*S*(I+i) - (sigma+nu)*E     (Exposed)
# Idot = sigma*E - (gamma+nu)*I              (Infectious)
# 
# Here beta(t) - the infection rate - is parameterized by a sinusoidal function
# plus a constant.
#
# Traditionally, there is an additional state
#
# Rdot = gamma*I - nu*R
#
# However, since we only observe I, R contributes nothing to the data fit and
# we have removed it from the system. 
#
# Other parameters are
# i - a visiting process
# nu - death rate
# sigma - the rate of movement from Exposed to Infectious.
# gamma - the rate of recovery from infection.
#
# It is generally more stable to solve the equations for the log states rather
# than the states themselves. CollocInfer contains a number of useful tools
# that let you make this transition without needing to re-code all of your
# differential equations.

library('CollocInfer')
            
            
#### Get some data and parameters

SEIRvarnames = SEIRvarnames
SEIRparnames = SEIRparnames

SEIRtimes = SEIRtimes
SEIRdata = SEIRdata

SEIRpars = SEIRpars


#### Now format the data so that S and E measurements are listed as NA

data = cbind(matrix(NA,length(SEIRdata),2),SEIRdata)

# We'll also look at the log observations

logdata = log(data)



#### define the right side evaluation function

SEIRfn = make.SEIR()


#### A couple of functions to define the infection rate

beta.fun = function(t,p,more){
    return( p['b0'] + p['b1']*sin(2*pi*t) + p['b2']*cos(2*pi*t) )
}

beta.dfdp = function(t,p,more){
    dfdp =  cbind(rep(1,length(t)), sin(2*pi*t), cos(2*pi*t)) 
    colnames(dfdp) = c('b0','b1','b2')
    return(dfdp)
}

betamore = list(beta.fun=beta.fun,
                beta.dfdp=beta.dfdp,
                beta.ind=c('b0','b1','b2'))


# Create a collocation basis to represent the state vector

rr = range(SEIRtimes)
knots = seq(rr[1],rr[2],2/52)
norder = 3
nbasis = length(knots)+norder-2

bbasis = create.bspline.basis(range=rr,norder=norder,nbasis=nbasis,breaks=knots)

# To get an initial estimate of the states we smooth the observed I component
# and set the other coefficients to zero.  

DEfd = smooth.basis(SEIRtimes,logdata[,3],fdPar(bbasis,1,0.1))

plotfit.fd(log(SEIRdata),SEIRtimes,DEfd$fd)

coefs = cbind(matrix(0,bbasis$nbasis,2),DEfd$fd$coefs)
DEfd = fd(coefs,bbasis)


# We will want to represent the state variables on the log scale so that they remain
# always positive. To do this, we set the 'posproc' component of LS.setup to 1. We
# will also compare the logstate to the log of the data directly, and therefore set
# 'poslik' to 0.

# We call LS.setup first and use the outputted lik and proc objects to pull the 
# coefficients for the unobserved state variables into line with the 
# differential equation. 

objs = LS.setup(SEIRpars,fn=SEIRfn,fd.obj=DEfd,more=betamore,data=data,times=SEIRtimes,
  posproc=1,poslik=0,names=SEIRvarnames,lambda=c(100,1,1))

proc = objs$proc
lik = objs$lik

res1 = FitMatchOpt(coefs=coefs,which=1:2,proc=proc,pars=SEIRpars,meth='nlminb')


# Let's have a look at the result

DEfd1 = fd(res1$coefs,bbasis)
plot(DEfd1,ylim=c(5,13))
points(SEIRtimes,logdata[,3])

# We can now run an initial smooth using the estimated coefficients as starting points. 

res2 = inneropt(data=logdata,times=SEIRtimes,pars=SEIRpars,proc=proc,lik=lik,coefs=res1$coefs,in.meth='nlminb')

# Has this changed much?

DEfd2 = fd(res2$coefs,bbasis)

plot(DEfd2,lwd=2,ylim=c(5,13))
lines(DEfd1)


# And we can call the optimizing functions. 

res3 = outeropt(data=logdata,times=SEIRtimes,pars=SEIRpars,proc=proc,lik=lik,coefs=res2$coefs,
  active=c('i','b0','b1','b2'))


# Some plots

## First, the estimated trajectories

DEfd3 = fd(res3$coefs,bbasis)
plot(DEfd3,lwd=2,ylim=c(5,14))


# Let's compare this to the data

plotfit.fd(logdata[,3],SEIRtimes,DEfd3[3],ylab='Fit to Data')


# We can also look at the discrepancy between the estimated trajectory and the
# differential equation

traj = eval.fd(SEIRtimes,DEfd3)     # estimated trajectory
colnames(traj) = SEIRvarnames

dtraj = eval.fd(SEIRtimes,DEfd3,1)  # derivative of the estimated trajectory

ftraj = proc$more$fn(SEIRtimes,traj,res3$pars,proc$more$more)   # Trajectory predicted by ODE


X11()
matplot(SEIRtimes,dtraj,type='l',lty=2,ylim =c(-10,10),ylab='SEIR derivatives' )
matplot(SEIRtimes,ftraj,type='l',lty=1,add=TRUE)

X11()
matplot(SEIRtimes,dtraj-ftraj,type='l',ylim=c(-4,4),ylab='Fit to Model')





## The alternative is to exponentiate the state before we compare to original data.
# This can take a very long time and is only recommended if you really need to do
# it, or have a couple of hours to wait. 


objs2 = LS.setup(SEIRpars,fn=SEIRfn,fd.obj=DEfd,more=betamore,data=data,times=SEIRtimes,
  posproc=1,poslik=1,names=SEIRvarnames,SEIRparnames,lambda=c(100,1,1))

lik2 = objs2$lik
proc2 = objs2$proc

res2 = inneropt(data=data,times=SEIRtimes,pars=res3$pars,proc=proc2,lik=lik2,coefs=res3$coefs)

res3 = outeropt(data=data,times=SEIRtimes,pars=res3$pars,proc=proc2,lik=lik2,coefs=res3$coefs,
  active=c('i','b0','b1','b2'))



###############################################################################
# Some more basic setup operations
#
# Here we go through the steps necessary to set up the SEIR equations manually. 
# This allows us several options for dealing with the positivity of the 
# state vector. 
#
# 1. We can ignore it and hope for the best. 
#
# 2. We can take a log transform of the ODE and then exponentiate the solutions
# to compare to the data
#
# 3. We can take a log transform of the ODE and compare this to the log data. 
###############################################################################


# First of all, we need the values of the basis at the observation times and the 
# quadrature points. 


qpts = 0.5*(knots[1:(length(knots)-1)] + knots[2:length(knots)])

bvals.obs = Matrix(eval.basis(times,bbasis),sparse=TRUE)

bvals = list(bvals = Matrix(eval.basis(qpts,bbasis),sparse=TRUE),
            dbvals = Matrix(eval.basis(qpts,bbasis,1),sparse=TRUE))


### This proc object is just the standard squared error deviation
# from the right hand side of a differential equation

sproc = make.SSEproc()
sproc$bvals = bvals
sproc$more = make.SEIR()
sproc$more$more = betamore
sproc$more$qpts = qpts
sproc$more$weights = matrix(1,length(qpts),3)%*%diag(c(1e2,1e0,1e0))
sproc$more$names = SEIRvarnames
sproc$more$parnames = SEIRparnames


### However, ODEs are often much more numerically stable if represented on a
# log scale. The make.logtrans() function will take the right hand side functions
# and derivatives defined for any differential equation and provide the equivalent
# system for the log state. Note that this does affect the way  you represent
# your state when considering the observations. 

lsproc = make.SSEproc()
lsproc$bvals = bvals
lsproc$more = make.logtrans()
lsproc$more$more = make.SEIR()
lsproc$more$more$more = betamore
lsproc$more$qpts = qpts
lsproc$more$weights = matrix(1,length(qpts),3)%*%diag(c(1e2,1e0,1e0))
lsproc$more$names = SEIRvarnames
lsproc$more$parnames = SEIRparnames

### Lik objects, this is the standard squared error. 

slik = make.SSElik()
slik$bvals = eval.basis(times,bbasis)
slik$more = make.id()
slik$more$weights = array(1,dim(data))
slik$more$names = SEIRvarnames
slik$more$parnames = SEIRparnames

# Log transform transformation for the lik object. For this, we note that we 
# have represented the trajectory on the log scale and will need to transform 
# back, this makes the numerics much harder and it can take a very long time to 
# converge.  

lslik = make.logstate.lik()
lslik$bvals = slik$bvals
lslik$more$weights = slik$more$weights
lslik$more = slik
lslik$more$parnames = SEIRparnames


# Numerically things work much better on the log scale

dfd = data2fd(logdata[,3],times,bbasis)

coefs = matrix(0,nbasis,3)
coefs[,3] = dfd$coefs
 
res = FitMatchOpt(coefs=coefs,which=1:2,proc=lsproc,pars=SEIRpars,meth='nlminb')

res2 = inneropt(data=logdata,times=SEIRtimes,pars=SEIRpars,proc=lsproc,lik=lslik,coefs=res$coefs)

res3 = outeropt(data=data,times=SEIRtimes,pars=SEIRpars,proc=lsproc,lik=lslik,coefs=res2$coefs,
  active=c('i','b0','b1','b2'))


## Or just log the observations, this is faster. 

res2 = inneropt(data=logdata,times=SEIRtimes,pars=SEIRpars,proc=lsproc,lik=slik,coefs=res$coefs)

res3 = outeropt(data=logdata,times=SEIRtimes,pars=SEIRpars,proc=lsproc,lik=slik,coefs=res2$coefs,
  active=c('i','b0','b1','b2'))

