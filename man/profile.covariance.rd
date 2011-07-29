\name{Profile.covariance}
\alias{Profile.covariance}
\title{Profile.covariance}
\description{ Newey-West estimate of covariance of parameter estimates from profiling. }
\usage{
Profile.covariance(pars,active=NULL,times,data,coefs,lik,proc,
          in.meth='nlminb',control.in=NULL,eps=1e-6,GN=FALSE)
}
\arguments{
\item{pars}{ Initial values of parameters to be estimated processes. }
\item{active}{ Incides indicating which parameters of \code{pars} should be estimated; defaults to all of them.  } 
\item{times}{ Vector observation times for the data.} 
\item{data}{  Matrix of observed data values. }
\item{coefs}{ Vector giving the current estimate of the coefficients in the spline. }
\item{lik}{ \code{lik} object defining the observation process. }
\item{proc}{ \code{proc} object defining the state process. }
\item{in.meth}{ Inner optimization function currently one of 'nlminb', 'MaxNR', 'optim' or 'house'. The last calls \code{SplineEst.NewtRaph}. This is fast but has poor convergence.  }
\item{control.in}{ Control object for inner optimization functions. }
\item{eps}{ Step-size for finite difference estimate of second derivatives.  }
\item{GN}{ Indicator of whether a Gauss-Newton approximation for the Hessian should be employed. Only valid for
least-squares methods. }
}
\value{Returns a Newey-West estimate of the covariance matrix of the parameter estimates. }
\details{Currently assumes a lag-5 auto-correlation among observation vectors. }
\seealso{ProfileErr, ProfileSSE, Profile.LS, Profile.multinorm}
\examples{

# See example in Profile.LS

}