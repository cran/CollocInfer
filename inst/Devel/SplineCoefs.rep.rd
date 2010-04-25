\name{SplineEst.rep}
\alias{SplineEst.NewtRaph.rep}
\alias{SplineCoefsList.rep}
\alias{SplineCoefsErr.rep}
\alias{SplineCoefsDC.rep}
\alias{SplineCoefsDP.rep}
\alias{SplineCoefsDC2.rep}
\alias{SplineCoefsDCDP.rep}
\title{Spline Estimation Functions for Repeated Processes}
\description{Model-based smoothing; estimation,  objective criterion and derivatives. }
\useage{
SplineEst.NewtRaph.rep(coefs,times,data,lik,proc,pars,whichpars,control=list(reltol=1e-12,maxit=1000,maxtry=10,trace=0))  
SplineCoefsList.rep(coefs,times,data,lik,proc,pars,whichpars,sgn=1)
SplineCoefsErr.rep(coefs,times,data,lik,proc,pars,whichpars,sgn=1)
SplineCoefsDC.rep(coefs,times,data,lik,proc,pars,whichpars,sgn=1)
SplineCoefsDP.rep(coefs,times,data,lik,proc,pars,whichpars,sgn=1)
SplineCoefsDC2.rep(coefs,times,data,lik,proc,pars,whichpars,sgn=1)
SplineCoefsDCDP.rep(coefs,times,data,lik,proc,pars,whichpars,sgn=1)
}
\arguments{
\item{coefs}{ Vector giving the current estimate of the coefficients in the spline. }
\item{times}{ Vector observation times for the data.} 
\item{data}{  Matrix of observed data values. }
\item{lik} { \code{lik} object defining the observation process. }
\item{proc}{ \code{proc} object defining the state process. }
\item{pars}{ Parameters to be used for the processes. }
\item{whichpars}{ A list of length the number of replicates, each elements specifies the indices of the vector \code{pars} that correspond to that replicate; defaults to the whole vector applying to each replicate. }
\item{sgn}{ Is the minimizing (1) or maximizing (0)? }
\item{control}{ A list giving control parameters for Newton-Raphson optimization. It should contain
  \item{reltol}{ Relative tollerance criterion for the gradient and improvement before termination. }
  \item{maxit}{ Maximum number of iterations. }
  \item{maxtry}{ Maximum number of halving-steps to try before declaring no improvement is possible. }
  \itme{trace}{ How much iteration history to output; 0 surpresses all output, a positive value outputs parameters and improvement at each 
iteration. }
}
}
\value{
\item{SplineEst.NewtRaph.rep}{Returns a list that is the result of the optimization with elements
  \item{value}{The final objective criterion.}
  \item{coefs}{The optimizing value of the coefficients.}
  \item{g}{The gradient at the optimizing value.}
  \item{H}{The Hessian at the optimizing value.}
}
\item{SplineCoefsList.rep}{Collates the gradient calculations and returns a list with elements
  \item{value}{Output of \code{SplineCoefsErr.rep}}
  \item{gradient}{Output of \code{SplineCoefsDC.rep}}
  \item{Hessian}{Output of \code{SplineCoefsDC2.rep}}
}
\item{SplineCoefsErr.rep}{The complete data log likelihood for the smooth; the inner optimization objective. }
\item{SplineCoefsDC.rep}{The derivative fo \code{SplineCoefsErr.rep} with respect to \code{coefs}.}
\item{SplineCoefsDP.rep}{The derivative fo \code{SplineCoefsErr.rep} with respect to \code{pars}.}
\item{SplineCoefsDC2.rep}{The second derivative fo \code{SplineCoefsErr.rep} with respect to \code{coefs}.}
\item{SplineCoefsDCDP.rep}{The second derivative fo \code{SplineCoefsErr.rep} with respect to \code{coefs} and \code{pars}.}
The output of gradients is in terms of an array with dimensions corresponding to derivatives. Derivatives with with respect to coefficients
are given in dimensions before those that give derivatives with respect to parameters.
}
\details{\code{SplineEst.NewtRaph.rep} performs a simple Newton-Raphson estimate for the optimal value of the coefficients.
This estimate lacks the convergence checks of other estimation packages, but may yeild a fast solution when needed.}
\seealso{}
\examples{}