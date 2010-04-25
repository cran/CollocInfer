##########
# This file contains definitions for the functions to obtain Newey-West based
# estimates of parameter estimates resulting for the generalized profile
# proceedure. 
#
# The main function is 'Profile.covariance'. Some utility functions are also defined. 

Profile.covariance <- function(pars,active=NULL,times,data,coefs,lik,proc,in.meth='nlminb',control.in=NULL,eps=1e-6)
{
    check.lik.proc.data.coefs(lik,proc,data,times,coefs)

    if(is.null(active)){ active = 1:length(pars) }

    apars = pars[active]
    
    H = matrix(0,length(apars),length(apars))

    g = ProfileDP(pars=apars,allpars=pars,times=times,data=data,coefs=coefs,lik=lik,proc=proc,active=active,sumlik=FALSE)
    gg = apply(g,2,sum)

    for(i in 1:length(apars)){        
      if(file.exists('curcoefs.tmp')){file.remove('curcoefs.tmp')}
      if(file.exists('optcoefs.tmp')){file.remove('optcoefs.tmp')}
      if(file.exists('counter.tmp')){file.remove('counter.tmp')}

        tpars = apars
        tpars[i] = tpars[i] + eps

        tf = ProfileErr(tpars,pars,times,data,coefs,lik,proc,in.meth=in.meth,control.in=control.in,active=active)  
        tg = ProfileDP(tpars,pars,times,data,coefs,lik,proc,active=active,sumlik=TRUE)                            

        H[,i] = (tg-gg)/eps

      if(file.exists('curcoefs.tmp')){file.remove('curcoefs.tmp')}
      if(file.exists('optcoefs.tmp')){file.remove('optcoefs.tmp')}
      if(file.exists('counter.tmp')){file.remove('counter.tmp')}
    }
               
    Covar = NeweyWest.Var( 0.5*(t(H)+H) ,g,5)
        
    return( Covar )
}



############################################################################################
#
# Some Utilities
#
############################################################################################


blocks2mat = function(H)  # List of matrices -> large matrix
{

    out = c()

    for(i in 1:length(H)){
      tout = H[[i]][[1]]
      if(length(H[[i]])>1){
        for(j in 2:length(H[[i]])){
            if(inherits(H[[i]][[j]],'dgCMatrix')|inherits(H[[i]][[j]],'dgeMatrix')){ tout=cBind(tout,H[[i]][[j]]) }
            else{ tout = cbind(tout,H[[i]][[j]]) }
        }
      }
      if(i > 1){ 
        if(inherits(tout,'dgCMatrix')|inherits(tout,'dgeMatrix')){ out=rBind(out,tout) }
        else{ out = rbind(out,tout) }
      } 
      else{ out = tout }
    }

    return(out)
}

## Newey West Calculations, with thanks to Steve Ellner


## GAUSS trimr function: trims n1 rows from the start and n2 rows from the end
## of a matrix or vector 
trimr <- function (a,n1,n2) {
        da<-dim(a); 
        if(is.null(da)) {a[(n1+1):(length(a)-n2)]}
        else {a[(n1+1):(da[1]-n2),]}
}

Newey.West <-function(x,y,maxlag) {
        w=1-(1:maxlag)/(maxlag+1); w=w/length(x); 
        out=mean(x*y); 
        for(i in 1:maxlag) {
            out=out+w[i]*sum(trimr(x,i,0)*trimr(y,0,i))+w[i]*sum(trimr(y,i,0)*trimr(x,0,i))
        }
        return(out)     
} 

NeweyWest.Var = function(H,g,maxlag)       
{
    V = solve(H)

    I = 0*H
    
    if(is.null(maxlag)){ 
        n = nrow(g)
        maxlag = max(5,n^(0.25))
    }

    if(maxlag > 0){
        for(i in 1:(ncol(g)-1)){
            for(j in (i+1):ncol(g)){
                I[i,j] = Newey.West(g[,i],g[,j],maxlag)
                I[j,i] = I[i,j]
            }    
        }
     
    }

    return( V%*%(I+ t(g)%*%g)%*%V )

}

