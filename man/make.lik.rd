\name{make.lik}
\alias{make.lik}
\alias{make.SSElik}
\alias{make.multinorm}
\alias{pomp.dmeasure}
\title{Observation Process Distribution Function}
\description{Returns a list of functions that calculate the observation process distribution and its derivatives;
designed to be used with the collocation inference functions. }
\usage{
make.SSElik()

make.multinorm()
}
\value{A list of functions that calculate the log observation distribution and its derivatives.
\item{make.SSElik}{ calculates weighted squared error between predictions
(given by \code{fn} in \code{more})  and observations}
\item{make.Multinorm}{ calculates a multivariate normal distribution.}
}
\details{These functions require \code{more} to be a list with elements:
\itemize{
\item{\code{fn}}{ The transform function of the states to observations, or to their derivatives.}
\item{\code{dfdx}}{ The derivative of \code{fn} with respect to states.}
\item{\code{dfdp}}{ The derivative of \code{fn} with respect to parameters. }
\item{\code{d2fdx2}}{ The second derivative of \code{fn} with respect to states.}
\item{\code{d2fdxdp}}{ The cross derivative of \code{fn} with respect to states and parameters.}
}

\code{make.Multinorm} further requires:
\itemize{
\item{\code{var.fn}}{ The variance given in terms of states and parameters. }
\item{\code{var.dfdx}}{ The derivative of \code{var.fn} with respect to states.}
\item{\code{var.dfdp}}{ The derivative of \code{var.fn} with respect to parameters. }
\item{\code{var.d2fdx2}}{ The second derivative of \code{var.fn} with respect to states.}
\item{\code{var.d2fdxdp}}{ The cross derivative of \code{var.fn} with respect to states and parameters.}
}

\code{make.SSElik} further requres \code{weights} giving weights to each observation.
}
\seealso{\code{\link{LS.setup}}, \code{\link{multinorm.setup}}}
\examples{

# Straightforward sum of squares:

lik = make.SSElik()
lik$more = make.id()

# Multivariate normal about an exponentiated state with constant variance

lik = make.multinorm()
lik$more = c(make.exp(),make.cvar())

}