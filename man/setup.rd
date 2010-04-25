\name{setup}
\alias{setup}
\alias{LS.setup}
\alias{multinorm.setup}
\title{Setup Functions for proc and lik objects}
\description{These functions set up lik and proc objects of squared error
and multinormal processes.}
\usage{
LS.setup(pars,coefs=NULL,fn,basisvals=NULL,lambda,fd.obj=NULL,
        more=NULL,data=NULL,weights=NULL,times=NULL,quadrature=NULL,eps=1e-6,
        posproc=0,poslik=0,discrete=0,names=NULL,sparse=FALSE)

multinorm.setup(pars,coefs=NULL,fn,basisvals=NULL,var=c(1,0.01),fd.obj=NULL,
        more=NULL,data=NULL,times=NULL,quadrature=NULL,eps=1e-6,posproc=0,
        poslik=0,discrete=0,names=NULL,sparse=FALSE)
}
\arguments{
\item{pars}{ Initial values of parameters to be estimated processes. }
\item{coefs}{ Vector giving the current estimate of the coefficients in the spline. }
\item{fn}{ A function giving the right hand side of a differential/difference equation.  The function should have arguments
\itemize{
  \item{times}{The times at which the RHS is being evaluated.}
  \item{x}{The state values at those times.}
  \item{p}{Parameters to be entered in the system.}
  \item{more}{An object containing additional inputs to \code{fn} }}
It should return a matrix of the same dimension of \code{x} giving the right hand side values.

If \code{fn} is given as a single function, its derivatives are estimated by finite-differencing with
stepsize \code{eps}. Alternatively, a list can be supplied with elements:
\itemize{
  \item{fn}{Function to calculate the right hand side should accept a matrix of state values at .}
  \item{dfdx}{Function to calculate the derivative with respect to \code{x}}
  \item{dfdp}{Function to calculate the derivative with respect to \code{p}}
  \item{d2fdx2}{Function to calculate the second derivative with respect to \code{x}}
  \item{d2fdxdp}{Function to calculate the second derivative with respect to \code{x} and \code{p}
}}
These functions take the same arguments as \code{fn} and should output multidimensional arrays with
the dimensions ordered according to time, state, deriv1, deriv2; here derivatives with respect to \code{x}
always precede derivatives with respect to \code{p}. 

\code{fn} can also be given as a \code{pomp} object (see the \code{pomp} package), in which case it is 
interfaced to \code{CollocInfer} through \code{pomp.skeleton} using a finite differencing. 
}
\item{basisvals}{Values of the collocation basis to be used. This can either be a basis object from the \code{fda} package,
or a list elements:
\itemize{
  \item{bvals.obs}{A matrix giving the values of the basis at the observation times}
  \item{bvals}{A matrix giving the values of the basis at the quadrature times}
  \item{dbvals}{A matrix giving the derivative of the basis at the quadrature times}
}
For discrete systems, it may also be specified as a matrix, in which case \code{bvals$bvals} is obtained by deleting the last row
and \code{bvals$dbvals} is obtained by deleting the first/  

If left as NULL, it is taken from \code{fd.obj} for \code{discrete = 0} and defaults to an identity matrix
of the same dimension as the number of observations for \code{discrete=1} systems. 
}
\item{lambda}{(\code{LS.setup} only) Penalty value trading off fidelity to data with fidelity to differential equations.}
\item{var}{(\code{profile.Cproc} or \code{profile.Dproc}) A vector of length 2, giving  }
\item{fd.obj}{(Optional) A functional data object; if this is non-null, \code{coefs} and \code{basisvals} is extracted from here. }
\item{more}{An object specifying additional arguments to \code{fn}. }
\item{data}{The data to be used, this can be a matrix, or a three-dimensional array. If the latter, the middle
dimension is taken to be replicates. The data are returned, if replicated they are returned in a concatenated form.}
\item{weights}{(\code{LS.setup} only)  }
\item{times}{ Vector observation times for the data. If the data are replicated, times are returned in a concatenated form.}
\item{quadrature}{ Quadrature points, should contain two elements (if not NULL)
\itemize{
  \item{qpts}{Quadrature points; defaults to midpoints between knots}
  \item{qwts}{Quadrature weights; defaults to normalizing by the length of \code{qpts}.   } 
}}
\item{eps}{ Finite differencing step size, if needed. }
\item{posproc}{ Should the state vector be constrained to be positive? If this is the case, the state is represented by
an exponentiated basis expansion in the \code{proc} object. }
\item{poslik}{ Should the state be exponentiated before being compared to the data? When the state is represented
on the log scale \code{posproc=1}, this is an alternative to taking the log of the data. }
\item{discrete}{ Is this a discrete or continuous-time system?}
\item{names}{ The names of the state variables if not given by the column names of \code{coefs}.}
\item{sparse}{ Should sparse matrices be used for basis values? This option can save memory when 
\code{ProfileGausNewt} and \code{SplineEstNewtRaph} are called. Otherwise sparse matrices will be converted to 
full matrices and this can slow the code down.}
}
\value{A list with elements
\item{coefs}{Starting values for \code{coefs}}
\item{lik}{The \code{lik} object generated}
\item{proc}{The \code{proc} item generated}
\item{data}{The data matrix, concatenated if from a 3d array.}
\item{times}{The vector of observation times, concatenated if data is a 3d array.}
}
\details{These functions provide basic setup utilities for the collocation inference methods. They define
\code{lik} and \code{proc} objects for sum of squared errors and multivariate normal log likelihoods with
nonlinear transfer functions describing the evolution of the state vector.
\itemize{
  \item{LS.setup}{Creates sum of squares functions}
  \item{multinorm.setup}{Creates multinormal log likelihoods for a continuous-time system.}
}}
\seealso{inneropt, outeropt, Profile.LS, Profile.multinorm, Smooth.LS, Smooth.multinorm}
\examples{

# FitzHugh-Nagumo

t = seq(0,20,0.05)            # Observation times

pars = c(0.2,0.2,3)           # Parameter vector
names(pars) = c('a','b','c')

knots = seq(0,20,0.2)         # Create a basis
norder = 3
nbasis = length(knots) + norder - 2
range = c(0,20)

bbasis = create.bspline.basis(range=range,nbasis=nbasis,norder=norder,breaks=knots)

lambda = 10000               # Penalty value

coefs = matrix(0,nbasis,2)   # Coefficient matrix

profile.obj = LS.setup(pars=pars,coefs=coefs,fn=make.fhn(),basisvals=bbasis,lambda=lambda,times=t)


# Using multinorm

var = c(1,0.01)

profile.obj = multinorm.setup(pars=pars,coefs=coefs,fn=make.fhn(),basisvals=bbasis,var=var,times=t)


# Henon - discrete

hpars = c(1.4,0.3)
t = 1:200

coefs = matrix(0,200,2)
lambda = 10000

profile.obj = LS.setup(pars=hpars,coefs=coefs,fn=make.Henon(),basisvals=NULL,lambda=lambda,times=t,discrete=1)
}