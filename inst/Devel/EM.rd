\name{EM}
\alias{EM}
\alias{Qfn}
\title{EM-Estimation for Collocation Inference}
\description{ EM algorithm for collocation estimation of parameters in continuous-time
stochastic processes.  }
\usage{
EM(pars,times,data,coefs,lik,proc,niter=10,
    in.meth,control.in,tfn=NULL,active=1:length(pars))

Qfn(pars,times,data,coefs,lik,proc,oldpars,
    tfn=NULL,active=1:length(pars))
}
\arguments{
  \item{pars}{ Initial values of parameters to be estimated processes. }
  \item{times}{ Vector observation times for the data.}
  \item{data}{  Matrix of observed data values. }
  \item{coefs}{ Vector giving the current estimate of the coefficients in the spline. }
  \item{lik} { \code{lik} object defining the observation process. }
  \item{proc}{ \code{proc} object defining the state process. }
  \item{niter}{ Number of EM steps to take. }
  \item{in.meth}{ Inner optimization function to be used, currently one of 'nlminb', 'MaxNR', 'optim' or 'house'.
    The last calls \code{SplineEst.NewtRaph}. This is fast but has poor convergence.  }
  \item{control.in}{ Control object for inner optimization functions. }
  \item{tfn}{ Name of a transformation function for parameters. }
  \item{active}{ Incides indicating which parameters of \code{pars} should be estimated; defaults to all of them.  }
  \item{oldpars}{ Starting parameter values for the Q-function in the EM algorithm.}
}
\value{
  \item{Qfn}{Returns the value of the Q-function in the Laplace-EM algorithm; entered into EM}
  \item{EM}{Returns the result of the EM algorithm after \code{niter} steps. It results in a list
  \itemize{
    \item{pars}{The (complete) optimized parameter vector}
    \item{coefs}{The coefficients optimized for \code{pars}}
    \item{Pres}{Optimization with respect to \code{pars[active]}}
    \item{Cres}{Optimization with respect to \code{coefs}}
  }}
}
\details{The Laplace-EM algorith proceeds by laplace approximations to the distribution of
\code{coefs} given \code{pars}. That is, we maximize the complete-data log likelihood for \code{coefs} and then
approximate this log likelihood with a second-order Taylor expansion about \code{coefs} at the maximizing value. The iteration has
the following steps

1. \code{C} = argmin \code{CDLL(pars,coefs)} = argmin \code{lik(pars,coefs)} + \code{proc(pars,coefs)}

2. \code{Qfn(pars,newpars)} =  \code{CDLL(pars,C)} + [d2\code{CDLL(pars,C)}/d\code{C}2][d2\code{CDLL(newpars,C)}/d\code{C}2]

3. \code{newpars} = argmin \code{Qfn(pars,newpars)}
}
\seealso{sse.setup, multinorm.setup, Profile.sse, Profile.multinorm, inneropt}
\examples{}