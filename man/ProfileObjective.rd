\name{ProfileObjective}
\alias{Profile.GausNewt}
\alias{ProfileSSE}
\alias{ProfileList}
\alias{ProfileErr}
\alias{ProfileDP}
\title{Profile Estimation with Collocation Inference}
\description{ Profile estimation and objective functions for collocation estimation of parameters in continuous-time
stochastic processes.  }
\usage{
Profile.GausNewt(pars,times,data,coefs,lik,proc,in.meth="nlminb",
      control.in=NULL,active=1:length(pars),
      control=list(reltol=1e-6,maxit=50,maxtry=15,trace=1))

ProfileSSE(pars,allpars,times,data,coefs,lik,proc,in.meth='nlminb',
      control.in=NULL,active=1:length(pars),dcdp=NULL,oldpars=NULL,
      use.nls=TRUE,sgn=1)
      
ProfileErr(pars,allpars,times,data,coefs,lik,proc,in.meth = "house",
      control.in=NULL,sgn=1,active=1:length(allpars))

ProfileDP(pars,allpars,times,data,coefs,lik,proc,in.meth = "house",
      control.in=NULL,sgn=1,sumlik=1,active=1:length(allpars))

ProfileList(pars,allpars,times,data,coefs,lik,proc,in.meth = "house",
      control.in=NULL,sgn=1,active=1:length(allpars))
}
\arguments{
\item{pars}{ Initial values of parameters to be estimated processes. }
\item{allpars}{ Full parameter vector including \code{pars}. Assignment \code{allpars[active] = pars} is always made.}
\item{times}{ Vector observation times for the data.} 
\item{data}{  Matrix of observed data values. }
\item{coefs}{ Vector giving the current estimate of the coefficients in the spline. }
\item{lik}{ \code{lik} object defining the observation process. }
\item{proc}{ \code{proc} object defining the state process. }
\item{in.meth}{ Inner optimization function currently one of 'nlminb', 'MaxNR', 'optim' or 'house'. The last calls \code{SplineEst.NewtRaph}.
        This is fast but has poor convergence.  }
\item{control.in}{ Control object for inner optimization function. }
\item{sgn}{ Is the minimizing (1) or maximizing (0)? }
\item{active}{ Incides indicating which parameters of \code{pars} should be estimated; defaults to all of them.  } 
\item{oldpars}{ Starting parameter values for the Q-function in the EM algorithm.}
\item{dcdp}{ Estimate for the gradient of the optimized coefficients with respect to parameters; used internally.  }
\item{use.nls}{ In ProfileSSE, is 'nls' being used in the outer-optimization?}
\item{sumlik}{ In ProfileDP and ProfileDP.AllPar; should the gradient be given for each observation, or summed over them? }
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
\item{Profile.GausNewt}{Output of a simple Gaus-Newton iteration to optimizing the objective function when the
observation likelihood takes the form of a sum of squared errors. Returns a list with the following elements:
\describe{
  \item{pars}{ The optimized value of the parameters.}
  \item{in.res}{ The result of the inner optimization.}
  \item{value}{ The value of the optimized sum of squared errors.}
}}
\item{ProfileSSE}{Output for the outer optimization when the observation likelihood is given by squared error. This is a list
with the following elements
\describe{
  \item{value}{ The value of the outer optimization criterion.}
  \item{gradient}{ The derivative of \code{f} with respect to \code{pars}.}
  \item{coefs}{ The optimized value of \code{coefs} for the current value of \code{pars}.}
  \item{dcdp}{ The derivative of the optimized value of \code{coefs} at the current value of \code{pars}.}
}}
\item{ProfileErr}{The outer optimization criterion in the general case: the value of the observation likelihood at the profiled
estimates of \code{coefs}.}
\item{ProfileDP}{The derivative of \code{ProfileErr} with respect to \code{allpars[active]}.}
\item{ProfileList}{Returns the results of ProfileErr and ProfileDP as a list with elements \code{value} and \code{gradient}}
}
\details{ \code{Profile.GausNewt} provides a simple implementation of Gauss-Newton optimization and may
not result in optimized values that are as good as other optimizers in \code{R}.

When  \code{Profile.GausNewt} is not being used, \code{ProfileSEE} and \code{ProfileERR} create the
following temporary files:
\describe{
\item{counter.tmp}{The number of halving-steps taken on the current optimization step.}
\item{curcoefs.tmp}{The current estimates of the coefficients.}
\item{optcoefs.tmp}{The optimal estimates of the coefficients in the iteration history.}  }
These need to be removed manually when the optimization is finished. \code{optcoefs.tmp} will contain
the optimal value of \code{coefs} for plotting the estimated trajectories.
}
\seealso{\code{\link{outeropt}}, \code{\link{Profile.LS}}, \code{\link{Profile.multinorm}}, \code{\link{LS.setup}}, \code{\link{multinorm.setup}}}
