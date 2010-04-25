\name{FhNdata}
\alias{FhNdata}
\alias{FhNtimes}
\alias{FhNpars}
\alias{FhNvarnames}
\alias{FhNparnames}
\title{FitzHugh-Nagumo data}
\description{Data generated for FitzHugh-Nagumo Examples}
\usage{
FhNdata
}
\format{
\itemize{
\item{FhNdata}{ A 41 by 2 matrix of data generated from the FitzHugh Nagumo equations.}
\item{FhNtimes}{ A vector of 41 observation times corresponding to FhNdata.}
\item{FhNpars}{ Named parameter vector used to generate \code{FhNdata}.}
\item{FhNvarnames}{ \code{c('V','R')}: the state variable names for the FitzHugh Nagumo system.}
\item{FhNparnames}{ \code{c('a','b','c')} parameter names for the FitzHugh Nagumo system.}
}}
\source{
James Ramsay, Giles Hooker David Campbell and Jiguo Cao, 2007.
"Parameter Estimation for Differential Equations: A Generalized Smoothing Approach".
Journal of the Royal Statistical Society Vol 69 No 5.
}
\seealso{}