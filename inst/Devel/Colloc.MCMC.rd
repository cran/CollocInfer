\name{Colloc.MCMC}
\alias{Colloc.MCMC}
\title{Bayesian MCMC for Collocation Inference}
\description{ An MCMC algorithm for Collocation Inference  }
\usage{
Colloc.MCMC(times,data,pars,coefs,lik,proc,prior,walk.var,
  nstep,in.meth='house',control.in=NULL)
}
\arguments{
  \item{times}{ Vector observation times for the data.}
  \item{data}{  Matrix of observed data values. }
    \item{pars}{ Initial values of parameters to be estimated processes. }
  \item{coefs}{ Vector giving the current estimate of the coefficients in the spline. }
  \item{lik} { \code{lik} object defining the observation process. }
  \item{proc}{ \code{proc} object defining the state process. }
  \item{niter}{ Number of MCMC steps to take. }
  \item{walk.var}{ Random walk variance for parameters. }
  \item{in.meth}{ Inner optimization function to be used, currently one of 'nlminb', 'MaxNR', 'optim' or 'house'.
    The last calls \code{SplineEst.NewtRaph}. This is fast but has poor convergence.  }
  \item{control.in}{ Control object for inner optimization functions. }
}
\value{
    \item{parhist}{The Markov chain history of the parameters}
    \item{coefhist}{The Markov chain history of the coefficients. }
  }
}
\details{
The MCMC algorithm proceeds by taking a mean-zero Gaussian random walk step in the 
parameter space, minimizes  the inner optimization in coefficients and then takes a random
walk in coefficient space with variance given by the second derivative matrix of the 
inner optimization. 
}
\seealso{sse.setup, multinorm.setup, Profile.sse, Profile.multinorm, inneropt}
\examples{}