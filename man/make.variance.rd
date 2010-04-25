\name{make.variance}
\alias{make.variance}
\alias{make.cvar}
\alias{make.var.SEIR}
\title{Variance Functions}
\description{Returns a list of functions that calculate a (possibly state and parameter dependent) variance. }
\usage{
make.cvar()

make.var.SEIR()
}
\arguments{}
\value{A list of functions that calculate a variance function and its derivatives,
in a form compatible with the collocation inference functions.

\item{make.cvar}{ returns a variance that is constant but may depend on parameters }

\item{make.var.SEIR}{ returns a state-dependent transition covariance matrix calculated for the SEIR equations.}
}
\details{\code{make.cvar} requires the specification of further elements in the list. In particular
the element \code{more} should be a list containing
}
\seealso{make.multinorm}
\examples{

# Multivariate normal observation of the state vector.

lik = make.multinorm()
lik$more = c(make.id(),make.cvar())

}