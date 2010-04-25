\name{make.proc}
\alias{make.Dproc}
\alias{make.Cproc}
\alias{make.SSEproc}
\title{Process Distributions}
\description{Functions to define process distributions in the collocation inference package.}
\usage{
make.Dproc()

make.Cproc()

make.SSEproc()
}
\arguments{}
\value{A list of functions that the process distribution
\item{make.Cproc}{ creates functions to evaluate the distribution of the derivative of
the state vector given the current state for continuous-time systems.  }       
\item{make.Dproc}{ creates functions to evaluate the distribution of the next time point of
the state vector given the current state for discrete-state systems.  }
\item{make.SSEproc}{ treats the distribution of the derivative as an independent gaussian and
cacluates weighted sums of squared errors between derivatives and the prediction from the current state.}
}
\details{All functions require \code{more} to specify this distribution. This should be a list containing
\itemize{
\item{\code{fn}}{ The distribution specified.}
\item{\code{dfdx}}{ The derivative of \code{fn} with respect to states.}
\item{\code{dfdp}}{ The derivative of \code{fn} with respect to parameters. }
\item{\code{d2fdx2}}{ The second derivative of \code{fn} with respect to states.}
\item{\code{d2fdxdp}}{ The cross derivative of \code{fn} with respect to states and parameters.}
}
For \code{Cproc} and \code{Dproc} this should specify the distribution; for \code{SSEproc} it
should specify the right hand side of a differential equation. 
}
\seealso{LS.setup, multinorm.setup}
\examples{

# FitzHugh-Nagumo Equations

proc = make.SSEproc()
proc$more = make.fhn()

# Henon Map

proc = make.Dproc()
proc$more = make.Henon


# SEIR with multivariate normal transitions

proc = make.Cproc()
proc$more = make.multinorm()
proc$more$more = c(make.SEIR(),make.var.SEIR())

}