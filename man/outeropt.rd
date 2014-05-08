\name{outeropt}
\alias{outeropt}
\title{Outer Optimization Functions}
\description{Outer optimization; performs profiled estimation.}
\usage{
outeropt(data,times,pars,coefs,lik,proc,
      in.meth='nlminb',out.meth='nlminb',
      control.in=list(),control.out=list(),active=1:length(pars))
}
\arguments{
\item{data}{  Matrix of observed data values. }
\item{times}{ Vector observation times for the data.}
\item{pars}{ Initial values of parameters to be estimated processes. }
\item{coefs}{ Vector giving the current estimate of the coefficients in the spline. }
\item{lik}{ \code{lik} object defining the observation process. }
\item{proc}{ \code{proc} object defining the state process. }
\item{in.meth}{ Inner optimization function currently one of 'nlminb', 'maxNR', 'optim' or 'SplineEst'. The last calls \code{SplineEst.NewtRaph}.
This is fast but has poor convergence.  }
\item{out.meth}{ Outer optimization function to be used, one of 'optim' (defaults to  BFGS routine in \code{optim} unless \code{control.out$meth}
specifies otherwise), 'nlminb', 'maxNR' #, 'trust'
or 'subplex'. When squared error is being used, 'ProfileGN' and 'nls' can also be given. The former of these calls \code{Profile.GausNewt},
a fast but naive Gauss-Newton solver. }
\item{control.in}{ Control object for inner optimization function. }
\item{control.out}{ Control object for outer optimization function. }
\item{active}{ Indices indicating which parameters of \code{pars} should be estimated; defaults to all of them.  }
}
\value{A list containing
\item{pars}{Optimized parameters}
\item{coefs}{Optimized coefficients at \code{pars}}
\item{res}{The result of the outer optimization.}
\item{counter}{A set of parameters and objective values for each successful iteration.}
}
\details{The outer optimization for parameters looks only at the objective defined by the \code{lik}
object. For every parameter value, \code{coefs} are optimized by \code{inneropt} and then the value of
\code{lik} for these coefficients is computed.

A number of optimization routines can be used here, some experimentation is recommended.  Libraries
for these optimization routines are not pre-loaded.   Where these functions take options as explicit arguments
instead of a list, they should be listed in \code{control.out} and will be called by their names.

The routine creates
temporary files 'curcoefs.tmp' and 'optcoefs.tmp' to update coefficients as \code{pars} evolves. These overwrite
existing files of those names and are deleted before the function terminates.
}
\seealso{\code{\link{inneropt}}, \code{\link{Profile.LS}}, \code{\link{ProfileSSE}}, \code{\link{ProfileErr}}, \code{\link{LS.setup}}, \code{\link{multinorm.setup}}}
\examples{

\dontrun{
data(FhNdata)

knots = seq(0,20,0.2)         # Create a basis
norder = 3
nbasis = length(knots) + norder - 2
range = c(0,20)

bbasis = create.bspline.basis(range=range,nbasis=nbasis,norder=norder,breaks=knots)

lambda = 10000               # Penalty value

DEfd = smooth.basis(FhNtimes,FhNdata,fdPar(bbasis,1,0.5))   # Smooth to estimate
                                                            # coefficients first
coefs = DEfd$fd$coefs
colnames(coefs) = FhNvarnames

profile.obj = LS.setup(pars=FhNpars,coefs=coefs,fn=make.fhn(),basisvals=bbasis,
      lambda=lambda,times=FhNtimes)

lik = profile.obj$lik
proc= profile.obj$proc

res = outeropt(data=FhNdata,times=FhNtimes,pars=FhNpars,coefs=coefs,lik=lik,proc=proc,
    in.meth="nlminb",out.meth="nlminb",control.in=NULL,control.out=NULL)


plot(res$coefs,main='outeropt')
print(blah)
}}