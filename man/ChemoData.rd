\name{ChemoData}
\alias{ChemoData}
\alias{ChemoTime}
\alias{ChemoPars}
\alias{ChemoVarnames}
\alias{ChemoParnames}
\title{Chemostat Example Data}
\description{Five-species Chemostat Model}
\usage{
ChemoData
}
\format{
\itemize{
\item{ChemoData}{ A 61 by 2 matrix of data observed in a chemostat.}
\item{ChemoTime}{ A vector of 61 observation times corresponding to ChemoData.}
\item{ChemoPars}{ Named parameter vector as a starting point for estimation \code{ChemoData}.}
\item{ChemoVarnames}{ \code{c('N','C1','C2','B','S')}: the state variable names for the chemostat system.}
\item{ChemoParnames}{ parameter names for the chemostat system.}
}}
\source{
Yoshida, T.,  L. E. Jones, S. P. Ellner, G. F. Fussmann and N. G. Hairston, 2003,
"Rapid evolution drives ecological dynamics in a predator-prey system", Nature, 424,
pp. 303-306.
}
\seealso{}