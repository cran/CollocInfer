IntegrateForward = function(y0,ts,pars,proc)
{
  parms = list(pars=pars,proc=proc)
  out = lsoda(y=y0,times=ts,func=oderhs,parms=parms)
  
  return(list(times=out[,1],states=out[,2:ncol(out)]))
}



oderhs = function(t,r,parms)
{
  proc = parms$proc
  pars = parms$pars
  r = matrix(r,1,length(r))
  colnames(r) = proc$more$names

  y = as.vector(proc$more$fn(t,r,pars,proc$more$more))

  return(list(y))
}