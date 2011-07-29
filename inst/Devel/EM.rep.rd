\name{EM.rep}
\alias{EM.rep}
\alias{Qfn.rep}
\title{EM-Estimation from Collocation Inference}
\description{ EM algorithm for collocation estimation of parameters in repeated continuous-time
stochastic processes.  }
\useage{
EM.rep(pars,times,data,coefs,lik,proc,niter=10,control.in,tfn=NULL,active=1:length(pars),whichpars)
Qfn.rep(pars,times,data,coefs,lik,proc,oldpars,tfn=NULL,active=1:length(pars),whichpars)
}
\arguments{
\item{pars}{ Initial values of parameters to be estimated processes. }
\item{times}{ Vector observation times for the data.} 
\item{data}{  Matrix of observed data values. }
\item{coefs}{ Vector giving the current estimate of the coefficients in the spline. }
\item{lik} { \code{lik} object defining the observation process. }
\item{proc}{ \code{proc} object defining the state process. }
\item{niter}{ Number of EM steps to take. }
\item{control.in}{ Control object for inner optimization functions. }
\item{tfn}{ Name of a transformation function for parameters. }
\item{active}{ Incides indicating which parameters of \code{pars} should be estimated; defaults to all of them.  } 
\item{oldpars}{ Starting parameter values for the Q-function in the EM algorithm.}
\item{whichpars}{ A list of length the number of replicates, each elements specifies the indices of the vector \code{allpar} that correspond to that replicate; defaults to the whole vector applying to each replicate. }
}
\value{
  \item{Qfn.rep}{Returns the value of the Q-function in the Laplace-EM algorithm; entered into EM}
  \item{EM.rep}{Returns the result of the EM algorithm after \code{niter} steps. It results in a list
    \item{Pres}{Optimization with respect to \code{pars[active]}}
    \item{Cres}{Optimization with respect to \code{coefs}}
  }
}
\details{The Laplace-EM algorith proceeds by laplace approximations to the distribution of
\code{coefs} given \code{pars}. That is, we maximize the complete-data log likelihood for \code{coefs} and then
approximate this log likelihood with a second-order Taylor expansion about \code{coefs} at the maximizing value. The iteration has
the following steps

1. \code{C} = argmin \code{CDLL(pars,coefs)} = argmin \code{lik(pars,coefs)} + \code{proc(pars,coefs)}
2. \code{Qfn.rep(pars,newpars)} =  \code{CDLL(pars,C)} + [d2\code{CDLL(pars,C)}/d\code{C}2]\[d2\code{CDLL(newpars,C)}/d\code{C}2]
3. \code{newpars} = argmin \code{Qfn.rep(pars,newpars)}
}
\seealso{}
\examples{}