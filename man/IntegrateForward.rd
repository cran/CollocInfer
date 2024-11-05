\name{IntegrateForward}
\alias{IntegrateForward}
\title{IntegrateForward}
\description{ Solves a differential equation going forward based on a \code{proc} object. }
\usage{
IntegrateForward(y0,ts,pars,proc,more)
}
\arguments{
\item{y0}{ Initial conditions to start from. }
\item{ts}{ Vector of time points at which to report values of the differential equation solution. }
\item{pars}{ Initial values of parameters to be estimated processes. }
\item{proc}{ Object defining the state process. This can either be a function evaluating
the right hand side of the differential equation or a \code{proc} object. If a \code{proc} object is given, \code{proc$more$fn} is assumed to give the right
hand side of the differential equation. }
\item{more}{ If \code{proc} is a function, this contains a list of additional inputs. }
}
\value{Returns the output from solving the differential equation using the \code{lsoda} routines. 
Specifically, it returns a list with elements
\describe{
\item{times}{ The output times. }
\item{states}{ The output states. }
}
}

\seealso{\code{\link{Profile.LS}}, \code{\link{Profile.multinorm}}}
\examples{
proc = make.SSEproc()
proc$more = make.fhn()
proc$more$names = c('V','R')

y0 = c(-1,1)
names(y0) = c('V','R')

pars = c(0.2,0.2,3)
names(pars) = c('a','b','c')

ts = seq(0,20,0.5)

value = IntegrateForward(y0,ts,pars,proc)

matplot(value$times,value$states)
}

