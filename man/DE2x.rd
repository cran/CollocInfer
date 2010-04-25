\name{DE2x}
\alias{DE2x}
\title{Solve ODEs}
\description{ Solve the differential equation given in \code{proc$more$fn} }
\usage{
DE2x(y0,times,pars,proc)
}
\arguments{
\item{y0}{ Vector of initial conditions for the ODE.  }
\item{times}{ Vector of time points to evaluate the ODE.  }
\item{pars}{ Vector of parameters } 
\item{proc}{  \code{proc} object for the system }
}
\value{
A matrix giving the value of the states at each time point. 
}
\details{This is a wrapper function for \code{lsoda}.}
\seealso{forward.prediction.error}
\examples{
##################################################################
#### Set up lik and proc objects for FitzHugh-Nagumo System  #####
##################################################################

data(FhNdata)

knots = seq(0,20,0.2)
norder = 3
nbasis = length(knots) + norder - 2
range = c(0,20)

bbasis = create.bspline.basis(range=range(FhNtimes),nbasis=nbasis,
	norder=norder,breaks=knots)

coefs = matrix(0,bbasis$nbasis,2)
colnames(coefs) = FhNvarnames

profile.obj = LS.setup(pars=FhNpars,coefs=coefs,fn=make.fhn(),basisvals=bbasis,lambda=1000,times=FhNtimes)
lik = profile.obj$lik
proc= profile.obj$proc

##############################################################
### Now we can call the ODE solver and compare to the data ###
##############################################################

y0 = c(-1,1)
x = DE2x(y0,FhNtimes,FhNpars,proc)

matplot(FhNtimes,FhNdata)
matplot(FhNtimes,x,type='l',add=TRUE)

}