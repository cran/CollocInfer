\name{Smooth.LS}
\alias{Smooth.LS}
\alias{Smooth.multinorm}
\title{Model-Based Smoothing Functions}
\description{Perform the inner optimization to estimate coefficients given parameters.}
\usage{
Smooth.LS(fn,data,times,pars,coefs=NULL,basisvals=NULL,lambda,fd.obj=NULL,
        more=NULL,weights=NULL,quadrature=NULL,in.meth='nlminb',
        control.in,eps=1e-6,posproc=FALSE,poslik=FALSE,discrete=FALSE,names=NULL,sparse=FALSE)
        
Smooth.multinorm(fn,data,times,pars,coefs=NULL,basisvals=NULL,var=c(1,0.01),
        fd.obj=NULL,more=NULL,quadrature=NULL,in.meth='nlminb',
        control.in,eps=1e-6,posproc=FALSE,poslik=FALSE,discrete=FALSE,names=NULL,sparse=FALSE)
}
\arguments{
\item{fn}{ A function giving the right hand side of a differential/difference equation.  The function should have arguments
\itemize{
  \item{times}{The times at which the RHS is being evaluated.}
  \item{x}{The state values at those times.}
  \item{p}{Parameters to be entered in the system.}
  \item{more}{An object containing additional inputs to \code{fn}}
}
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
always precede derivatives with respect to \code{p}. }
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
\item{lambda}{(\code{Smooth.LS} only) Penalty value trading off fidelity to data with fidelity to differential equations.}
\item{var}{(\code{Smooth.multinorm}) A vector of length 2, giving  }
\item{fd.obj}{(Optional) A functional data object; if this is non-null, \code{coefs} and \code{basisvals} is extracted from here. }
\item{more}{An object specifying additional arguments to \code{fn}. }
\item{weights}{(\code{Smooth.LS} only)  }
\item{quadrature}{ Quadrature points, should contain two elements (if not NULL)
\itemize{
  \item{qpts}{Quadrature points; defaults to midpoints between knots}
  \item{qwts}{Quadrature weights; defaults to normalizing by the length of \code{qpts}.   } 
}}
\item{in.meth}{ Inner optimization function to be used, currently one of 'nlminb', 'MaxNR', 'optim' or 'SplineEst'.
The last calls \code{SplineEst.NewtRaph}. This is fast but has poor convergence.  }
\item{control.in}{ Control object for inner optimization function. }
\item{eps}{ Finite differencing step size, if needed. }
\item{posproc}{ Should the state vector be constrained to be positive? If this is the case, the state is represented by
an exponentiated basis expansion in the \code{proc} object. }
\item{poslik}{ Should the state be exponentiated before being compared to the data? When the state is represented
on the log scale (\code{posproc=TRUE}), this is an alternative to taking the log of the data. }
\item{discrete}{ Is this a discrete or continuous-time system?}
\item{names}{ The names of the state variables if not given by the column names of \code{coefs}.}
\item{sparse}{ Should sparse matrices be used for basis values? This option can save memory when 
\code{ProfileGausNewt} and \code{SplineEstNewtRaph} are called. Otherwise sparse matrices will be converted to 
full matrices and this can slow the code down.}
}
\value{A list with elements
\item{coefs}{Optimized coefficients at \code{pars}}
\item{lik}{The \code{lik} object generated}
\item{proc}{The \code{proc} item generated}
\item{res}{The result of the optimization method}
\item{data}{The data used in doing the fitting.}
\item{times}{The vector of times at which the observations were made}
}
\details{These routines create \code{lik} and \code{proc} objects and call \code{inneropt}.}
\seealso{inneropt, LS.setup, multinorm.setup, SplineCoefsErr}
\examples{
# Henon system

hpars = c(1.4,0.3)              # Parameters
t = 1:200

x = c(-1,1)                     # Create some dataa
X = matrix(0,200+20,2)
X[1,] = x

for(i in 2:(200+20)){ X[i,] = make.Henon()$ode(i,X[i-1,],hpars,NULL) }

X = X[20+1:200,]

Y = X + 0.05*matrix(rnorm(200*2),200,2)

basisvals = diag(rep(1,200))    # Basis is just identiy
coefs = matrix(0,200,2)


# For sum of squared errors

lambda = 10000

res1	= Smooth.LS(make.Henon(),data=Y,times=t,pars=hpars,coefs,basisvals=basisvals,
  lambda=lambda,in.meth='nlminb',discrete=TRUE)

# For multinormal transitions

var = c(1,0.01)

res2 = Smooth.multinorm(make.Henon(),data=Y,t,pars=hpars,coefs,basisvals=NULL,
  var=var,in.meth='nlminb',discrete=TRUE)
  
}