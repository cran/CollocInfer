\name{ChemoRMData}
\alias{ChemoRMData}
\alias{ChemoRMTime}
\alias{ChemoRMPars}
\alias{RMvarnames}
\alias{RMparnames}
\title{Rosenzweig-MacArthur Model Applied to Chemostat Data}
\description{Two-Species Rosenzweig-MacArthur Model}
\usage{
ChemoRMData
}
\format{
\describe{
\item{ChemoRMData}{ A 108 by 2 matrix of data observed in a chemostat.}
\item{ChemoRMPars}{ Named parameter vector as a starting point for estimation in \code{ChemoRMData}.}
\item{ChemoRMTime}{ A vector of 108 observation times corresponding to ChemoData.}
\item{RMparnames}{ parameter names for the Rosenzweig-MacArthur system.}
\item{RMvarnames}{  the state variable names for the  Rosenzweig-MacArthur system.}
}}
\source{
Becks, L., S. P. Ellner, L. E. Jones,  and N. G. Hairston, 2010,
"Reduction of adaptive genetic diversity radically alters eco-evolutionary community dynamics", Ecology Letters, 13,
pp. 989-997.
}
