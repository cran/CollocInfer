\name{FitMatch.rep}
\alias{FitMatchErr.rep}
\alias{FitMatchDC.rep}
\alias{FitMatchDC2.rep}
\title{Estimating Hidden States}
\description{Estimating hidden states of repeated processes to maximize agreement with the process.  }
\useage{
FitMatchErr.rep(coefs,allcoefs,which,pars,whichpars,proc,sgn=1)
FitMatchDC.rep(coefs,allcoefs,which,pars,whichpars,proc,sgn=1)
FitMatchDC2.rep(coefs,allcoefs,which,pars,whichpars,proc,sgn=1)
}
\arguments{
\item{coefs}{ Vector giving the current estimate of the coefficients for the hidden states. }
\item{allcoefs}{ Matrix giving the coefficients of all the states including initial values for \code{coefs}.} 
\item{which}{  Vector of indices of states to be estimated. }
\item{pars}{ Parameters to be used for the processes. }
\item{whichpars}{ A list of length the number of replicates, each elements specifies the indices of the vector \code{pars} that correspond to that replicate; defaults to the whole vector applying to each replicate. }
\item{proc}{ \code{proc} object defining the state process. }
\item{sgn}{ Is the minimizing (1) or maximizing (0)? }
}
\value{
\item{FitMatchErr.rep}{The value of the process likelihood at the current estimated states.}
\item{FitMatchDC.rep}{The derivative fo \code{FitMatchErr.rep} with respect to the elements \code{coefs} for the states being estimated.}
\item{FitMatchDC2.rep}{The second derivative fo \code{FitMatchErr.rep} with respect to the elements \code{coefs} for the states being estimated.}
}
\details{}
\seealso{}
\examples{}