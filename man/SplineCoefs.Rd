\name{SplineEst}
\alias{SplineEst.NewtRaph}
\alias{SplineCoefsList}
\alias{SplineCoefsErr}
\alias{SplineCoefsDC}
\alias{SplineCoefsDP}
\alias{SplineCoefsDC2}
\alias{SplineCoefsDCDP}
\title{Spline Estimation Functions}
\description{Model-based smoothing; estimation,  objective criterion and derivatives. }
\usage{
SplineEst.NewtRaph(coefs,times,data,lik,proc,pars,
      control=list(reltol=1e-12,maxit=1000,maxtry=10,trace=0))
      
SplineCoefsList(coefs,times,data,lik,proc,pars,sgn=1)

SplineCoefsErr(coefs,times,data,lik,proc,pars,sgn=1)

SplineCoefsDC(coefs,times,data,lik,proc,pars,sgn=1)

SplineCoefsDP(coefs,times,data,lik,proc,pars,sgn=1)

SplineCoefsDC2(coefs,times,data,lik,proc,pars,sgn=1)

SplineCoefsDCDP(coefs,times,data,lik,proc,pars,sgn=1)
}
\arguments{
\item{coefs}{ Vector giving the current estimate of the coefficients in the spline. }
\item{times}{ Vector observation times for the data.} 
\item{data}{  Matrix of observed data values. }
\item{lik}{ \code{lik} object defining the observation process. }
\item{proc}{ \code{proc} object defining the state process. }
\item{pars}{ Parameters to be used for the processes. }
\item{sgn}{ Is the minimizing (1) or maximizing (0)? }
\item{control}{ A list giving control parameters for Newton-Raphson optimization. It should contain
\describe{
  \item{reltol}{ Relative tollerance criterion for the gradient and improvement before termination. }
  \item{maxit}{ Maximum number of iterations. }
  \item{maxtry}{ Maximum number of halving-steps to try before declaring no improvement is possible. }
  \item{trace}{ How much iteration history to output; 0 surpresses all output, a positive value outputs parameters and improvement at each
iteration. }
}}
}
\value{
\item{SplineEst.NewtRaph}{Returns a list that is the result of the optimization with elements
\describe{
  \item{value}{ The final objective criterion.}
  \item{coefs}{ The optimizing value of the coefficients.}
  \item{g}{ The gradient at the optimizing value.}
  \item{H}{ The Hessian at the optimizing value.}
}}
\item{SplineCoefsList}{Collates the gradient calculations and returns a list with elements
\describe{
  \item{value}{ Output of \code{SplineCoefsErr}}
  \item{gradient}{ Output of \code{SplineCoefsDC}}
  \item{Hessian}{ Output of \code{SplineCoefsDC2}}
}}
\item{SplineCoefsErr}{The complete data log likelihood for the smooth; the inner optimization objective. }
\item{SplineCoefsDC}{The derivative of \code{SplineCoefsErr} with respect to \code{coefs}.}
\item{SplineCoefsDP}{The derivative of \code{SplineCoefsErr} with respect to \code{pars}.}
\item{SplineCoefsDC2}{The second derivative of \code{SplineCoefsErr} with respect to \code{coefs}.}
\item{SplineCoefsDCDP}{The second derivative of \code{SplineCoefsErr} with respect to \code{coefs} and \code{pars}.}
The output of gradients is in terms of an array with dimensions corresponding to derivatives. Derivatives with with respect to coefficients
are given in dimensions before those that give derivatives with respect to parameters.
}
\details{\code{SplineEst.NewtRaph} performs a simple Newton-Raphson estimate for the optimal value of the coefficients.
This estimate lacks the convergence checks of other estimation packages, but may yeild a fast solution when needed.}
\seealso{\code{\link{inneropt}}, \code{\link{Smooth.LS}}}
