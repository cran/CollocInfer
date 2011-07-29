CollocInferPlots = function(coefs,pars,lik,proc,times=NULL,data=NULL,
      cols=NULL,datacols=NULL,datanames=NULL,ObsPlot=TRUE,DerivPlot=TRUE,
      cex.axis=1.5,cex.lab=1.5,cex=1.5,lwd=2)
{
  if(is.null(cols)){ cols = 1:ncol(coefs) }



 timevec = proc$more$qpts

 traj = proc$bvals$bvals%*%coefs
 colnames(traj) = proc$more$names
 
 dtraj = proc$bvals$dbvals%*%coefs

 ftraj = proc$more$fn(timevec,traj,pars,proc$more$more)
 
 otraj = lik$more$fn(timevec,traj,pars,lik$more$more)
 
 
 
 if(ObsPlot){
  if(is.null(times)){ times = 1:nrow(data) }
  if(is.null(datacols)){  
    if(ncol(otraj) == ncol(traj)){ datacols = cols}
    else{ datacols = 1:ncol(otraj) }
  }
  
  dev.new()
  if(cex.lab > 1.5){ par(mar=c(5,5,1,1)) }
  matplot(timevec,otraj,type='l',lwd=lwd,xlab='time',ylab='Observations',
        cex.lab=cex.lab,cex.axis=cex.axis,col=datacols)
  if(!is.null(data)){
    if(is.null(datanames)){ datanames = substring(colnames(data),1,1) }
    if(length(datanames)==0){ datanames = as.character(1:ncol(data)) }
    if(is.null(times)){ times = 1:nrow(data) }
    matplot(times,data,col=datacols,add=TRUE,pch=datanames,cex=cex)
  }
 }
  
 if(DerivPlot){
  dev.new()
  if(cex.lab>1.5){par(mar=c(5,5,1,1))}
  matplot(timevec,ftraj,type='l',lwd=lwd,xlab='time',ylab='f(x): -, dx: --',
        cex.lab=cex.lab,cex.axis=cex.axis,col=cols,lty=1) 
  matplot(timevec,dtraj,type='l',lwd=lwd,lty=2,add=TRUE,col=cols)
  abline(h = 0)
  legend(x='topright',legend=proc$more$names,lwd=2,lty=1,col=cols)
 
  dev.new()
  if(cex.lab>1.5){par(mar=c(5,5,1,1))}
  par(mfrow=c(2,1))
  matplot(timevec,dtraj-ftraj,type='l',lwd=lwd,xlab='time',ylab='dx-f(x)',
        cex.lab=cex.lab,cex.axis=cex.axis,col=cols,lty=1)
  abline(h=0)     
  matplot(timevec,traj,type='l',lwd=lwd,xlab='time',ylab='x',
        cex.lab=cex.lab,cex.axis=cex.axis,col=cols,lty=1)    
  legend(x='topright',legend=proc$more$names,lwd=lwd,lty=1,col=cols,cex=cex)
 }

 return(list(timevec=timevec,traj=traj,dtraj=dtraj,ftraj=ftraj,otraj=otraj))

}