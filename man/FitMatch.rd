\name{FitMatch}
\alias{FitMatchOpt}
\alias{FitMatchErr}
\alias{FitMatchDC}
\alias{FitMatchDC2}
\alias{FitMatchList}
\title{Estimating Hidden States}
\description{Estimating hidden states to maximize agreement with the process.  }
\usage{
FitMatchOpt(coefs,which,pars,proc,meth='nlminb',control=list())

FitMatchErr(coefs,allcoefs,which,pars,proc,sgn=1)

FitMatchDC(coefs,allcoefs,which,pars,proc,sgn=1)

FitMatchDC2(coefs,allcoefs,which,pars,proc,sgn=1)

FitMatchList(coefs,allcoefs,which,pars,proc,sgn=1)
}
\arguments{
\item{coefs}{ Vector giving the current estimate of the coefficients for the hidden states. }
\item{allcoefs}{ Matrix giving the coefficients of all the states including initial values for \code{coefs}.} 
\item{which}{  Vector of indices of states to be estimated. }
\item{pars}{ Parameters to be used for the processes. }
\item{proc}{ \code{proc} object defining the state process. }
\item{sgn}{ Is the minimizing (1) or maximizing (0)? }
\item{meth}{ Optimization function currently one of 'nlminb', 'MaxNR', 'optim' or 'trust'. }
\item{control}{ Control object for optimization function. }
}
\value{
\item{FitMatchOpt}{A list containing
\describe{
  \item{coefs}{ The optimized coefficients for all states.}
  \item{res}{ The output of the optimization routine.}
}}
\item{FitMatchErr}{The value of the process likelihood at the current estimated states.}
\item{FitMatchDC}{The derivative of \code{FitMatchErr} with respect to the elements \code{coefs} for the states being estimated.}
\item{FitMatchDC2}{The second derivative of \code{FitMatchErr} with respect to the elements \code{coefs} for the states being estimated.}
\item{FitMatchList}{Returns a list with elements \code{value}, \code{gradient} and \code{hessian} given by the output of
\code{FitMatchErr}, \code{FitMatchDC} and \code{FitMatchDC2}.}
}
\details{These routines allow the values of coefficients for some states to be optimized relative to the others. That is, the
objective defined by \code{proc} is minimized over those states specified in \code{which} leaving the others constant. This
would be typically done, for example, a smooth is taken to estimate some states non-parametrically, but data is not available on all
of them.

A number of optimization routines have been implemented in \code{FitMatchOpt}, some experimentation is advised.}
\seealso{\code{\link{ParsMatchErr}}, \code{\link{SplineCoefsErr}}, \code{\link{inneropt}}}
\examples{
###############################
#### Some Data            #####
###############################

data(FhNdata)

# And parameter estimates

data(FhNest)


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

fd.data = FhNdata[,1]

DEfd = smooth.basis(FhNtimes,fd.data,fdPar(bbasis,1,0.5))   

coefs = cbind(DEfd$fd$coefs,rep(0,nbasis))
colnames(coefs) = FhNvarnames

#############################################################
### If We Only Observe One State, We Can Re-Smooth Others ### 
#############################################################

profile.obj = LS.setup(pars=FhNpars,coefs=coefs,fn=make.fhn(),
                      basisvals=bbasis,lambda=1000,times=FhNtimes)
lik = profile.obj$lik
proc= profile.obj$proc

#  DD = Matrix(diag(1,200),sparse=TRUE)
#  tDD = t(DD)

fres = FitMatchOpt(coefs=coefs,which=2,pars=FhNpars,proc)

plot(fd(fres$coefs,bbasis))
}
