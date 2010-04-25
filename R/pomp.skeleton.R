pomp.skeleton = function(times,y,p,more)
{
  # Turns a skeleton function from a 'pomp' object into the right hand side
  # of and ODE for use in CollocInfer
  x <- array(t(y),dim=c(ncol(y),1,nrow(y)),dimnames=list(colnames(y),NULL,NULL))
  params <- array(data=p,dim=c(length(p),1),dimnames=list(names(p),NULL))

  y = skeleton(more$pomp.obj,x=x,params=params,t=times)
  return(t(y[,1,]))
}
