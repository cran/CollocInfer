\name{ParsMatch}
\alias{ParsMatchOpt}
\alias{ParsMatchErr}
\alias{ParsMatchDP}
\alias{ParsMatchList}
\title{Estimate of Parameters from Smooth}
\description{ Objective function and derivatives to estimate parameters with a fixed smooth.  }
\usage{
ParsMatchOpt(pars,coefs,proc,active=1:length(pars),meth='nlminb',control=list())

ParsMatchErr(pars,coefs,proc,active=1:length(pars),allpars,sgn=1)

ParsMatchDP(pars,coefs,proc,active=1:length(pars),allpars,sgn=1)

ParsMatchList(pars,coefs,proc,active=1:length(pars),allpars,sgn=1)
}
\arguments{
\item{pars}{ Initial values of parameters to be estimated processes. }
\item{coefs}{ Vector giving the current estimate of the coefficients in the spline. }
\item{proc}{ \code{proc} object defining the state process. }
\item{active}{ Incides indicating which parameters of \code{allpar} should be estimated; defaults to all of them.  } 
\item{allpars}{  Vector of all parameters, the assignment \code{allpar[active]=pars} is made initially.   }
\item{sgn}{ Is the minimizing (1) or maximizing (0)? }
\item{meth}{ Optimization function currently one of 'nlminb', 'MaxNR', 'optim' or 'trust'. }
\item{control}{ Control object for optimization function. }
}
\value{
\item{ParsMatchOpt}{A list containing:
\itemize{
  \item{pars}{The entire parameter vector after optimization.}
  \item{res}{The output of the optimization routine.}
}}
\item{ParsMatchErr}{The value of the process likelihood at the current estimated states.}
\item{ParsMatchDP}{The derivative fo \code{ParsMatchErr} with respect to \code{pars[active]}.}
\item{ParsMatchList}{A list with entries \code{value} and \code{gradient} given by the output
of \code{ParsMatchErr} and \code{ParsMatchDP} respectively.}
}
\details{These routines fix the estimated states at the value given by \code{coefs}
and estimate \code{pars} to maximize agreement between the fixed state and the objective
given by the \code{proc} object.

A number of optimization routines have been implemented in \code{FitMatchOpt}, some experimentation is advised.
}
\seealso{FitMatchErr, SplineCoefsErr, inneropt}
\examples{

data(FhNdata)

###############################
####  Basis Object      #######
###############################

knots = seq(0,20,0.2)
norder = 3
nbasis = length(knots) + norder - 2
range = c(0,20)

bbasis = create.bspline.basis(range=range(FhNtimes),nbasis=nbasis,
	norder=norder,breaks=knots)


# Initial values for coefficients will be obtained by smoothing

DEfd = smooth.basis(FhNtimes,FhNdata,fdPar(bbasis,1,0.5))   # Smooth to estimate
                                                            # coefficients first
coefs = DEfd$fd$coefs
colnames(coefs) = FhNvarnames



#################################
### Initial Parameter Guesses ###
#################################

profile.obj = LS.setup(pars=FhNpars,coefs=coefs,fn=make.fhn(),basisvals=bbasis,
    lambda=1000,times=FhNtimes)
lik = profile.obj$lik
proc= profile.obj$proc

pres = ParsMatchOpt(FhNpars,coefs,proc)

npars = pres$pars

}