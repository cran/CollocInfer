\name{Profiling Routines}
\alias{Profile.LS}
\alias{Profile.multinorm}
\title{Profile Estimation Functions}
\description{These functions are wrappers that create lik and proc functions
and run generalized profiling.}
\usage{
Profile.LS(fn,data,times,pars,coefs=NULL,basisvals=NULL,lambda,
                        fd.obj=NULL,more=NULL,weights=NULL,quadrature=NULL,
                        likfn = make.id(), likmore = NULL,
                        in.meth='nlminb',out.meth='nls',
                        control.in,control.out,eps=1e-6,active=NULL,posproc=FALSE,
                        poslik=FALSE,discrete=FALSE,names=NULL,sparse=FALSE)

Profile.multinorm(fn,data,times,pars,coefs=NULL,basisvals=NULL,var=c(1,0.01),
                        fd.obj=NULL,more=NULL,quadrature=NULL,
                        in.meth='nlminb',out.meth='optim',
                        control.in,control.out,eps=1e-6,active=NULL,
                        posproc=FALSE,poslik=FALSE,discrete=FALSE,names=NULL,sparse=FALSE)
}
\arguments{
\item{fn}{ A function giving the right hand side of a differential/difference equation.  The function should have arguments
\itemize{
  \item{times}{The times at which the RHS is being evaluated.}
  \item{x}{The state values at those times.}
  \item{p}{Parameters to be entered in the system.}
  \item{more}{An object containing additional inputs to \code{fn}
}}
It should return a matrix of the same dimension of \code{x} giving the right hand side values.

If \code{fn} is given as a single function, its derivatives are estimated by finite-differencing with
stepsize \code{eps}. Alternatively, a list can be supplied with elements:
\itemize{
  \item{fn}{Function to calculate the right hand side should accept a matrix of state values at .}
  \item{dfdx}{Function to calculate the derivative with respect to \code{x}}
  \item{dfdp}{Function to calculate the derivative with respect to \code{p}}
  \item{d2fdx2}{Function to calculate the second derivative with respect to \code{x}}
  \item{d2fdxdp}{Function to calculate the second derivative with respect to \code{x} and \code{p}}
}
These functions take the same arguments as \code{fn} and should output multidimensional arrays with
the dimensions ordered according to time, state, deriv1, deriv2; here derivatives with respect to \code{x}
always precede derivatives with respect to \code{p}.}
\item{data}{  Matrix of observed data values. }
\item{times}{ Vector observation times for the data.}
\item{pars}{ Initial values of parameters to be estimated processes. }
\item{coefs}{ Vector giving the current estimate of the coefficients in the spline. }
\item{basisvals}{Values of the collocation basis to be used. This can either be a basis object from the \code{fda} package,
or a list elements:
\itemize{
  \item{bvals.obs}{A matrix giving the values of the basis at the observation times}
  \item{bvals}{A matrix giving the values of the basis at the quadrature times}
  \item{dbvals}{A matrix giving the derivative of the basis at the quadrature times}
}}
\item{lambda}{(\code{Profile.LS} only) Penalty value trading off fidelity to data with fidelity to differential equations.}
\item{var}{(\code{profile.Cproc} or \code{profile.Dproc}) A vector of length 2, giving  }
\item{fd.obj}{(Optional) A functional data object; if this is non-null, \code{coefs} and \code{basisvals} is extracted from here. }
\item{more}{An object specifying additional arguments to \code{fn}. }
\item{weights}{(\code{Profile.LS} only)  }
\item{quadrature}{ Quadrature points, should contain two elements (if not NULL)
\itemize{
  \item{qpts}{Quadrature points; defaults to midpoints between knots}
  \item{qwts}{Quadrature weights; defaults to normalizing by the length of \code{qpts}.   }
}}
\item{in.meth}{ Inner optimization function to be used, currently one of 'nlminb', 'MaxNR', 'optim' or 'house'.
The last calls \code{SplineEst.NewtRaph}. This is fast but has poor convergence.  }
\item{out.meth}{ Outer optimization function to be used, depending on the type of method
\itemize{
  \item{\code{Profile.LS}}{One of 'nls' or 'ProfileGN'; the latter calls \code{Profile.GausNewt} which is fast but
may have poor convergence.}
  \item{\code{Profile.multinorm}}{ One of 'optim' (defaults to  BFGS routine in \code{optim} unless \code{control.out$meth}
specifies otherwise), 'nlminb', or 'maxNR'.}
}}
\item{control.in}{ Control object for inner optimization function. }
\item{control.out}{ Control object for outer optimization function. }
\item{eps}{ Finite differencing step size, if needed. }
\item{active}{ Incides indicating which parameters of \code{pars} should be estimated; defaults to all of them.  }
\item{posproc}{ Should the state vector be constrained to be positive? If this is the case, the state is represented by
an exponentiated basis expansion in the \code{proc} object. }
\item{poslik}{ Should the state be exponentiated before being compared to the data? When the state is represented
on the log scale (\code{posproc=TRUE}), this is an alternative to taking the log of the data. }
\item{discrete}{ Is this a discrete-time or a continuous-time system? If discrete, the derivative is instead
taken to be the value at the next time point. }
\item{names}{ The names of the state variables if not given by the column names of \code{coefs}.}
\item{sparse}{ Should sparse matrices be used for basis values? This option can save memory when 
\code{ProfileGausNewt} and \code{SplineEstNewtRaph} are called. Otherwise sparse matrices will be converted to 
full matrices and this can slow the code down.}
\item{likfn}{ Defines a map from the trajectory to the observations. This should be in the same form as
\code{fn}. If a function is given, derivatives are estimated by finite differencing, otherwise a list
is expected to provide the same derivatives as \code{fn}. If \code{poslik=TRUE}, the states are
exponentiated before the \code{likfn} is evaluated and the derivatives are updated to account for this.
Defaults to the identity transform. }
\item{likmore}{ A list containing additional inputs to \code{likfn} if needed, otherwise set to \code{NULL} }
}
\value{A list with elements
\item{pars}{Optimized parameters}
\item{coefs}{Optimized coefficients at \code{pars}}
\item{lik}{The \code{lik} object generated}
\item{proc}{The \code{proc} item generated}
\item{data}{The data used in doing the fitting.}
\item{times}{The vector of times at which the observations were made}
}
\details{These functional all carry out the profiled optimization method of Ramsay et al 2007.
\code{Profile.LS} uses a sum of squared errors criteria for both fit to data and the fit of the derivatives
to a differential equation. \code{Profile.multinorm} uses multivariate normal approximations.
\code{discrete} changes the system to a discrete-time difference equation with the right hand side function
giving the transition function.

Note that these all call \code{outeropt}, which creates
temporary files 'curcoefs.tmp' and 'optcoefs.tmp' to update coefficients as \code{pars} evolves. These overwrite
existing files of those names and are deleted before the function terminates.
}
\seealso{outeropt, ProfileErr, ProfileSSE, LS.setup, multinorm.setup}
\examples{

###############################
####   Data             #######
###############################

data(FhNdata)

###############################
####  Basis Object      #######
###############################

knots = seq(0,20,0.2)
norder = 3
nbasis = length(knots) + norder - 2
range = c(0,20)

bbasis = create.bspline.basis(range=range(FhNtimes),nbasis=nbasis,
	norder=norder,breaks=knots)


#### Start from pre-estimated values to speed up optimization

data(FhNest)

spars = FhNestPars
coefs = FhNestCoefs

lambda = 10000

res1 = Profile.LS(make.fhn(),data=FhNdata,times=FhNtimes,pars=spars,coefs=coefs,
  basisvals=bbasis,lambda=lambda,in.meth='nlminb',out.meth='nls')

Covar = Profile.covariance(pars=res1$pars,times=FhNtimes,data=FhNdata,
  coefs=res1$coefs,lik=res1$lik,proc=res1$proc) 


\dontrun{
## Alternative, starting from perturbed coefficients -- takes too long for 
# automatic checks in CRAN

# Initial values for coefficients will be obtained by smoothing

DEfd = smooth.basis(FhNtimes,FhNdata,fdPar(bbasis,1,0.5))   # Smooth to estimate
                                                            # coefficients first
coefs = DEfd$fd$coefs
colnames(coefs) = FhNvarnames


###############################
####     Optimization       ###
###############################

spars = c(0.25,0.15,2.5)          # Perturbed parameters
names(spars)=FhNparnames
lambda = 10000

res1 = Profile.LS(make.fhn(),data=FhNdata,times=FhNtimes,pars=spars,coefs=coefs,
  basisvals=bbasis,lambda=lambda,in.meth='nlminb',out.meth='nls')

par(mfrow=c(2,1))
plotfit.fd(FhNdata,FhNtimes,fd(res1$coefs,bbasis))
}  
  
 
 
  
\dontrun{
####################################################
###  An Explicitly Multivariate Normal Formation ### 
####################################################

var = c(1,0.0001)

res2 = Profile.multinorm(make.fhn(),FhNdata,FhNtimes,pars=res1$pars,
          res1$coefs,bbasis,var=var,out.meth='nlminb', in.meth='nlminb')
}}