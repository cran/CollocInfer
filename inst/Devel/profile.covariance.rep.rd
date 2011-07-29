\name{profile.covariance.rep}
\alias{profile.covariance.rep}
\title{profile.covariance.rep}
\description{ Newey-West estimate of covariance of parameter estimates from profiling for repeated processes. }
\useage{
profile.covariance.rep(pars,active=NULL,times,data,coefs,lik,proc,whichpars,in.meth,control.in,eps=1e-6)
}
\arguments{
\item{pars}{ Initial values of parameters to be estimated processes. }
\item{active}{ Incides indicating which parameters of \code{pars} should be estimated; defaults to all of them.  } 
\item{times}{ Vector observation times for the data.} 
\item{data}{  Matrix of observed data values. }
\item{coefs}{ Vector giving the current estimate of the coefficients in the spline. }
\item{lik} { \code{lik} object defining the observation process. }
\item{proc}{ \code{proc} object defining the state process. }
\item{whichpars}{ A list of length the number of replicates, each elements specifies the indices of the vector \code{allpar} that correspond to that replicate; defaults to the whole vector applying to each replicate. }
\item{in.meth}{ Inner optimization function currently one of 'nlminb', 'MaxNR', 'optim' or 'house'. The last calls \code{SplineEst.NewtRaph}. This is fast but has poor convergence.  }
\item{control.in}{ Control object for inner optimization functions. }
\item{eps}{ Step-size for finite difference estimate of second derivatives.  }
}
\value{Returns a Newey-West estimate of the covariance matrix of the parameter estimates. }
\details{Currently assumes a lag-5 auto-correlation among observation vectors. }
\seealso{}
\examples{}