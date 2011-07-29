\name{make.findif}
\alias{make.findif}
\alias{make.findif.ode}
\alias{make.findif.loglik}
\alias{make.findif.var}
\title{Finite Difference Functions}
\description{Returns a list of functions that calculate finite difference derivatives.  }
\usage{
make.findif.ode()

make.findif.loglik()

make.findif.var()
}
\value{A list of functions that calculate the derivatives via finite differencing schemes. 

\item{make.findif.ode}{ calculates finite differences of a transform.}
\item{make.findif.loglik}{ returns the finite differences to a calculated log likelihood; used
within \code{lik} objects, or as \code{more} arguments to \code{Cproc} or \code{Dproc}.}
\item{make.findif.var}{ finite difference approximations to variances; mostly used in
the \code{Multinorm} functions.}
}
\details{All these functions require the sepcification of \code{more$eps} to give the size of the
finite differencing step. They also require \code{more} to specify the original object (ODE right hand side functions,
definitions of \code{lik} and \code{proc} objects).}
\seealso{LS.setup, multinorm.setup}
\examples{

# Sum of squared errors with finite differencing to get right-hand-side derivatives

proc = make.SSEproc()
proc$more = make.findif.ode()


# Finite differencing for the log likelihood

lik = make.findif.loglik()
lik$more = make.SSElik()


# Multivariate normal transitions with finite differencing for mean and variance functions

lik = make.multinorm()
lik$more = c(make.findif.ode,make.findif.var)

# Finite differencing for transition density of a discrete time system

proc = make.Dproc()
proc$more = make.findif.loglik()

}