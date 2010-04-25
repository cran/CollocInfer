\name{ParsMatch.rep}
\alias{ParsMatchErr.rep}
\alias{ParsMatchDP.rep}
\title{Estimate of Parameters from Smooth}
\description{ Objective function and derivatives to estimate parameters with a fixed smooth of repeated processes.  }
\useage{
ParsMatchErr.rep(pars,coefs,proc,active=1:length(pars),allpars,whichpars,sgn=1)
ParsMatchDP.rep(pars,coefs,proc,active=1:length(pars),allpars,whichpars,sgn=1)
}
\arguments{
\item{pars}{ Initial values of parameters to be estimated processes. }
\item{coefs}{ Vector giving the current estimate of the coefficients in the spline. }
\item{proc}{ \code{proc} object defining the state process. }
\item{active}{ Incides indicating which parameters of \code{allpar} should be estimated; defaults to all of them.  } 
\item{allpar}{  Vector of all parameters, the assignment \code{allpar[active]=pars} is made initially.   }
\item{whichpars}{ A list of length the number of replicates, each elements specifies the indices of the vector \code{allpar} that correspond to that replicate; defaults to the whole vector applying to each replicate. }
\item{sgn}{ Is the minimizing (1) or maximizing (0)? }
}
\value{
\item{ParsMatchErr.rep}{The value of the process likelihood at the current estimated states.}
\item{ParsMatchDP.rep}{The derivative fo \code{ParsMatchErr.rep} with respect to \code{pars[active]}.}
}
\details{}
\seealso{}
\examples{}