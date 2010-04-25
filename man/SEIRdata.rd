\name{SEIRdata}
\alias{SEIRdata}
\alias{SEIRtimes}
\alias{SEIRpars}
\alias{SEIRvarnames}
\alias{SEIRparnames}
\title{SEIR data}
\description{Data generated for SEIR Examples}
\usage{
SEIRdata
}
\format{
\itemize{
\item{SEIRdata}{ A 261 by 1 matrix of data generated from the SEIR Gillespie process run over 5 years.}
\item{SEIRtimes}{ A vector of 261 observation times corresponding to SEIRdata.}
\item{SEIRpars}{ Named parameter vector used to generate \code{SEIRdata}.}
\item{SEIRvarnames}{ \code{c('V','R')}: the state variable names for the SEIR system.}
\item{SEIRparnames}{  parameter names for the SEIR system.}
}}
\source{
Giles Hooker, Stephen P. Ellner, David Earn and Laura Roditi, 2010.
"Parameterizing State-space Models for Infectious Disease Dynamics by Generalized Profiling: Measles in Ontario",
Technical Report, Cornell University. 
}
\seealso{}