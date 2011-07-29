\name{inneropt}
\alias{inneropt}
\title{Inner Optimization Functions}
\description{Estmates coefficients given parameters.}
\usage{
inneropt(data,times,pars,coefs,lik,proc,in.meth='nlminb',control.in=list())
}
\arguments{
\item{data}{  Matrix of observed data values. }
\item{times}{ Vector observation times for the data.}
\item{pars}{ Initial values of parameters to be estimated processes. }
\item{coefs}{ Vector giving the current estimate of the coefficients in the spline. }
\item{lik}{ \code{lik} object defining the observation process. }
\item{proc}{ \code{proc} object defining the state process. }
\item{in.meth}{ Inner optimization function currently one of 'nlminb', 'maxNR', 'optim', 'trust' or 'SplineEst'.
The last calls \code{SplineEst.NewtRaph}. This is fast but has poor convergence.  }
\item{control.in}{ Control object for inner optimization function. }
}
\value{A list with elements
\item{coefs}{A matrix giving he optimized coefficients.}
\item{res}{The results of the inner optimization function.} }
\details{This minimizes the objective function defined by the addition of the \code{lik}
and \code{proc} objectives with respect to the coefficients. A number of generic
optimization routines can be used and some experimentation is recommended. }
\seealso{outeropt, smooth.LS, smooth.Cproc, smooth.Dproc, LS.setup, multinorm.setup, SplineCoefsErr}
\examples{
\dontrun{
# FitzHugh-Nagumo Equations

data(FhNdata)   # Some data

knots = seq(0,20,0.2)         # Create a basis
norder = 3
nbasis = length(knots) + norder - 2
range = c(0,20)

bbasis = create.bspline.basis(range=range(FhNtimes),nbasis=nbasis,norder=norder,breaks=knots)

lambda = 10000               # Penalty value

DEfd = smooth.basis(FhNtimes,FhNdata,fdPar(bbasis,1,0.5))   # Smooth to estimate
                                                            # coefficients first
coefs = DEfd$fd$coefs
colnames(coefs) = FhNvarnames

profile.obj = LS.setup(pars=FhNpars,coefs=coefs,fn=make.fhn(),basisvals=bbasis,lambda=lambda,times=FhNtimes)

lik = profile.obj$lik
proc= profile.obj$proc

res = inneropt(FhNdata,times=FhNtimes,FhNpars,coefs,lik,proc,in.meth='nlminb')

plot(fd(res$coefs,bbasis))
}}