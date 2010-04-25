Colloc.MCMC = function(times,data,pars,coefs,lik,proc,prior,walk.var,cscale,nstep,in.meth='SplineEst',control.in=NULL)
{
    parhist = matrix(0,nstep,length(pars))                     # Store parameter history
    coefhist = matrix(0,nstep,length(coefs))
    bestcoefhist = matrix(0,nstep,length(coefs))
    acc = rep(0,nstep)
    
 

    
    f = SplineCoefsErr(coefs,times,data,lik,proc,pars,sgn=-1)   # Complete data log posterior
#print(proc$fn(coefs,bvals,pars,more))
#    H = SplineCoefsDC2(coefs,times,data,lik,proc,pars,sgn=-1)
#    G = SplineCoefsDCDP(coefs,times,data,lik,proc,pars,sgn=-1)
    
#    dcdp = ginv(H)%*%G
 
  tcoefs = coefs
   
        
#     ei = eigen(tH)
 #    eta = rnorm(length(coefs))
  #   logq = sum(eta^2)/2   

#     coefs = coefs + ei$vectors%*%diag(1/sqrt(ei$values*(ei$values>0)))%*%eta            
 

               # Initialize
    bestcoefhist[1,]=tcoefs
    coefhist[1,] = tcoefs




    eta = 0        # cscale = 0 defaults to Collocation MCMC
    logq = 0
    for(i in 2:nstep){
     #   tcoefs=coefs

                                                   # Now a Gaussian error about the minimum
       tH = -1*SplineCoefsDC2(tcoefs,times,data,lik,proc,pars,sgn=-1)
         
          
       if(cscale>0){
         ei = eigen(tH)
         eta = rnorm(length(coefs))
  #       logq = -sum(eta^2)/2 - sum(log(ei$values[ei$values>0])/2)
          logq = sum(dnorm(eta,log=TRUE))
         
         eta = cscale*ei$vectors%*%diag(1/sqrt(ei$values*(ei$values>0)))%*%eta
       }
       
       tf =  SplineCoefsErr(as.vector(tcoefs)+eta,times,data,lik,proc,pars,sgn=-1) 
         
                      
#                 print(proc$fn(matrix(as.vector(tcoefs),nrow=ncol(lik$bvals)),bvals,pars,more)) 
print(c(tf,f ))
       if( runif(1) < exp( tf-f)){             # Acceptance probability
          coefs = as.vector(tcoefs)+eta
         
                    
          f = tf
#          H = SplineCoefsDC2(tcoefs,times,data,lik,proc,pars,sgn=-1)
#          G = SplineCoefsDCDP(tcoefs,times,data,lik,proc,pars,sgn=-1)
#          dcdp = ginv(H)%*%G

          acc[i] = 1
       }
       
       
       coefhist[i,] = as.vector(tcoefs)+eta
       bestcoefhist[i,] = coefs
    }

   return(list(parhist=parhist, coefhist=coefhist, acc = acc, bch=bestcoefhist))
}