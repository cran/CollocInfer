\name{forward.prediction.error}
\alias{forward.prediction.error}
\title{forward.prediction.error}
\description{ Forward prediction error objective for choice of lambda in square error criteria.  }
\usage{
forward.prediction.error(times,data,coefs,lik,proc,pars,whichtimes=NULL)
}
\arguments{
  \item{times}{ Vector observation times for the data.}
  \item{data}{  Matrix of observed data values. }
  \item{coefs}{ Vector giving the current estimate of the coefficients in the spline. }
  \item{lik}{ \code{lik} object defining the observation process. }
  \item{proc}{ \code{proc} object defining the state process. }
  \item{pars}{ Initial values of parameters to be estimated processes. }
  \item{whichtimes}{ Specifies the start and end times for forward prediction, given by indeces of \code{times}. This can be one of
    \describe{
    \item{list}{ each element of the list is itself a list of length 2; the first element gives the
      starting time to use and the second is a vector giving the prediction times.}
    \item{matrix}{ the first column giving the starting times and the second giving the ending times.}
    }
  If left NULL, \code{whichtimes} defaults to predicting one observation ahead from each observation.
  }
}
\value{The forwards prediction error from the estimates.  }
\details{Forward prediction error can be used to choose values of \code{lambda} in the profiled
estimation routines. The ordinary differential equation is solved starting from the starting
times specified in \code{whichtimes} and measured at the corresponding measurement times. The error is then recorded.
This should then be minimized by a grid search. }
\seealso{ \code{\link{ProfileSSE}}, \code{\link{outeropt}}}
