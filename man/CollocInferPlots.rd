\name{CollocInferPlots}
\alias{CollocInferPlots}
\title{Diagnostic PLots for CollocInfer}
\description{Diagnostic Plots on the Results of CollocInfer}
\usage{
CollocInferPlots(coefs,pars,lik,proc,times=NULL,data=NULL,
      cols=NULL,datacols=NULL,datanames=NULL,ObsPlot=TRUE,DerivPlot=TRUE,
      DerivResid=TRUE,newplot=FALSE,cex.axis=1.5,cex.lab=1.5,cex=1.5,lwd=2)
}
\arguments{
\item{coefs}{ Vector giving the current estimate of the coefficients. }
\item{pars}{ Vector of estimated parameters. }
\item{lik}{ \code{lik} object defining the observation process. }
\item{proc}{ \code{proc} object defining the state process. }
\item{times}{ Vector observation times for the data.}
\item{data}{  Matrix of observed data values. }
\item{cols}{ Optional vector specifying a color for each state variable. }
\item{datacols}{ Optional vector specifying a color for each observation dimension. }
\item{datanames}{ Optional character vector specifying a glyph to plot the data. Taken from the column-names
of \code{data} if not given. }
\item{ObsPlot}{ Should a plot of predictions and observations be given? }
\item{DerivPlot}{ Should derivative diagnostics be produced? }
\item{DerivResid}{ Should a plot of the difference between Dx and f(x) be produced?}
\item{newplot}{ Should plots be opened in new devices? Note if this is FALSE, DerivePlot will overwrite ObsPlot. }
\item{cex.axis}{ Axis font size. }
\item{cex.lab}{ Label font size. }
\item{cex}{ Plotting point font size }
\item{lwd}{ Plotting line width }
}
\value{
A list containing elements used in plotting:
\item{timevec}{ Times at which the trajectories etc were evaluated. }
\item{traj}{ Estimated value of the trajectory. }
\item{dtraj}{ Derivative of the estimated trajectory. }
\item{ftraj}{ Value of the derivative of the trajectory predicted by \code{proc} }
\item{otraj}{ Predicted values of the observations from \code{lik}. }
}
\details{ Timevec is taken to be the quadrature values. Three plots can be produced:

If \code{ObsPlot=TRUE} a plot is given of the predicted values of the observations along with
the observations themselves (if given). 

If \code{DerivPlot=TRUE} two plots are produced. The first gives the value of the derivative of 
the estimated trajectory (dashed) and the value of the right-hand-side of the ordinary differential equation
in \code{proc} (hence the predicted derivative) (solid).  The second plot gives their difference in the first panel as well as the estimated trajectory in the second panel. 
}

