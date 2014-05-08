\name{make.logtrans}
\alias{make.logtrans}
\alias{make.logstate.lik}
\alias{make.exp.Cproc}
\alias{make.exp.Dproc}
\title{Log Transforms}
\description{Functions to modify liklihood, transform, \code{lik} and \code{proc}
objects so that the operate with the state defined on a log scale.}
\usage{
make.logtrans()

make.exptrans()

make.logstate.lik()

make.exp.Cproc()

make.exp.Dproc()
}
\value{A list of functions that calculate log transforms and derivatives in various contexts.

\item{make.logtrans}{ modifies the right hand side of a differential equation and its
derivatives for a loged state vector. }
\item{make.exptrans}{ modfies a map from states to observations to a map from logged states
to observations along with its derivatives. }
\item{make.logstate.lik}{ modifies a \code{lik} object for state vectors given on the log scale. }
\item{make.exp.Cproc}{ \code{Cproc} with the state given on the log scale.}
\item{make.exp.Dproc}{ \code{Dproc} with the state given on the log scale.}
}
\details{ All functions require \code{more} to specify the original object (ODE right hand side functions,
definitions of \code{lik} and \code{proc} objects).}
\seealso{\code{\link{LS.setup}}, \code{\link{make.Cproc}}, \code{\link{make.Dproc}}}
\examples{

# Model the log of an SEIR process

proc = make.SSEproc()
proc$more = make.logtrans()
proc$more$more = make.SEIR()

# Observe a linear combination  of

lik = make.logstate.lik()
lik$more = make.SSElik()
lik$more$more = make.genlin()

# SEIR Model with multivariate transition densities

proc = make.exp.Cproc()
proc$more = make.multinorm()
proc$more$more = c(make.SEIR(),make.cvar())

}