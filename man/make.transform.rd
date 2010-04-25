\name{make.transfer}
\alias{make.transfer}
\alias{make.id}
\alias{make.exp}
\alias{make.genlin}
\alias{make.fhn}
\alias{make.Henon}
\alias{make.SEIR}
\alias{make.NS}
\alias{make.diagnostics}
\alias{chemo.fun}
\alias{pomp.skeleton}
\title{Transfer Functions}
\description{Returns a list of functions that calculate the transform and its derivatives. }
\usage{
make.id()

make.exp()

make.genlin()

make.fhn()

make.Henon()

make.SEIR()

make.NS()

chemo.fun(times,y,p,more=NULL)

pomp.skeleton(times,y,p,more)
}
\arguments{
All the functions 
created by \code{make...} functions, require the arguments needed by  \code{chemo.fun}
\item{times}{ the evaluation times}
\item{y}{ values of the state at the evaluation times }
\item{p}{ parameters to be used }
\item{more}{ a list of additional arguments, in this case \code{NULL}, for 
\code{pomp.sekelton}, \code{more} should be a list containing a \code{pomp} object
in the element \code{pomp.obj}. }
}
\value{A list of functions that calculate the transform and its derivatives,
in a form compatible with the collocation inference functions.
\item{make.id}{ returns the identity transform.}
\item{make.exp}{ returns the exponential transform.}
\item{make.genlin}{ returns a linear combination transform -- see details section below.}
\item{make.fhn}{ returns the FitzHugh-Nagumo equations.}
\item{make.Henon}{ reutrns the Henon map.}
\item{make.SEIR}{ returns SEIR equations for estimating the shape of a seasonal forcing component. }
\item{make.diagnostics}{ functions to perform forcing function diagnostics. }
}
\details{\code{make.genlin} requires the specification of further elements in the list. In particular
the element \code{more} should be a list containing
\itemize{
  \item{\code{mat}}{ a matrix defining the linear transform before any parameters are added.
This may be all zero, but it may also specify fixed elements, if desired. }
  \item{\code{sub}}{ a k-by-3 matrix indicating which parameters should be entered into
which elements of \code{mat}. Each row is a triple giving the row and colum of \code{mat} to be 
specified and the element of the parameter vector that should be substituted. \code{sub} over-rides
any values in \code{mat}.}
  \item{\code{force}}{ if input functions are given, these are given as a list.}
  \item{\code{force.mat}}{ specifying the influence of the elements of \code{force} on the state
variables. Defined as in \code{mat}.}
  \item{\code{force.sub}}{ defined as in \code{sub}, over-rides the elements of \code{force.mat} with
parameter values.}
}

\code{make.diagnostics} estimates forcing-function diagnostics as in Hooker, 2009 for 
goodness-of-fit assessment. It requires
\itemize{
  \item{\code{psi}}{ Values of a basis expansion for forcing functions at the quadrature points. }
  \item{\code{which}}{ Which states are to be forced? }
  \item{\code{fn}, \code{dfdx}, \code{d2fdx2}}{ Functions and derivatives as would be used to estimate 
  parameters for the original equations. }
  \item{pars}{ Parameters to go into \code{more$fn}. }
}

\code{make.SEIR} estimates parameters and a seasonal variation in the infection rate in an 
SEIR model.  It requires the specification of the seasonal change rate in \code{more} by
a list with objects
\itemize{
  \item{\code{beta.fun}}{ A function to calculate beta, it should have arguments \code{t}, 
    \code{p} and \code{betadef} and return a matrix giving the value of beta at times \code{t} 
    with parameters \code{p}. }
  \item{\code{beta.dfdp}}{ Should calculate the derivative of \code{beta.fun} with respect to \code{p},
    at times \code{t} returning a matrix. The matrix should be of size \code{length(t)} by
    \code{length(p)} where \code{p} is the entire parameter vector. }
  \item{\code{betadef}}{ Additional inputs (eg bases) to \code{beta.fun} and \code{beta.dfdp}.}
}

\code{make.NS} provides functions for the North Shore example. This is a possibly time-varying
forced linear system of one dimension. It requires \code{more} to specify \code{betabasis} to
describe the autoregressive coefficient, and \code{alphabasis} to provide a contant in front of
the functional data object \code{rainfd}.

\code{chemo.fun} Is a five-state predator-prey-resources model used as an example. It stands
alone as a function and should be used with the \code{findif.ode} functions. 
}
\seealso{LS.setup, multinorm.setup}
\examples{

# Observe the FitzHugh-Nagumo equations

proc = make.SSEproc()
proc$more = make.fhn()

lik = make.SSElik()
lik$more = make.id()

# Observe an unknown scalar transform of each component of a Henon map, given
# in the first two elements of the parameter vector:

proc = make.Dproc()
proc$more = make.multinorm()
proc$more$more = c(make.Henon,make.cvar)

lik = make.multinorm()
lik$more = c(make.genlin,make.cvar)
lik$more$more = list(mat = matrix(0,2,2),sub=matrix(c(1,1,1,2,2,2),2,3,byrow=TRUE))

# Model SEIR equations on the log scale and then exponentiate

lik = make.SSElik()
lik$more = make.exp()

proc = make.SSEproc()
proc$more = make.logtrans()
proc$more$more = make.SEIR()

}