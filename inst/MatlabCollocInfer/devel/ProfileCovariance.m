%%%%%%%%%%
% This file conta=s def=itions for the functions to obta= Newey-West based
% estimates of parameter estimates result=g for the generalized profile
% proceedure_ 
%
% The ma= function is 'Profile_covariance'_ Some utility functions are also def=ed_ 

function Covar = ProfileCovariance(pars,active,times,data,coefs,lik,proc,in_meth,control_in,eps,GN)
%%% CHECK DEFAULTS

check_lik_proc_data_coefs(lik,proc,data,times,coefs)

    if isempty(active)
        active = 1:length(pars); 
    end

    apars = pars(active);
   
    g = ProfileDP(apars,pars,times,data,coefs,lik,proc,active,0);  %%% CHECK ARGUMENTS HERE
    if isnumeric(g) == 0,
        g = matrix(g,length(g),length(active));
    end

    if GN == 0

      H = zeros(length(apars),length(apars));
      gg = ProfileDP(apars,pars,times,data,coefs,lik,proc,active,1);
      for i = 1:length(apars)
          
          %%%% NEED TO DEAL WITH OPTIMIZATION HERE
          
          tpars = apars;
          tpars(i) = tpars(i) + eps;

          tf = ProfileErr(tpars,pars,times,data,coefs,lik,proc,in_meth,control_in,active);
          tg = ProfileDP(tpars,pars,times,data,coefs,lik,proc,active,1);

          H(:,i) = (tg-gg)/eps;
      end
    else
        H = g(:,active)'*g(:,active);
    end
           
    Covar = NeweyWest_Var( 0.5*(t(H)+H) ,g,5);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Some Utilities
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Newey West Calculations, with thanks to Steve Ellner


%% GAUSS trimr function: trims n1 rows from the start and n2 rows from the end
%% of a matrix or vector 
function a = trimr(a,n1,n2) 
        da=size(a); 
        if isempty(da) 
            a = a((n1+1):(length(a)-n2));
        else
          a = a((n1+1):(da(1)-n2),:);
        end
end

function out = Newey_West(x,y,maxlag) 
        w=1-(1:maxlag)/(maxlag+1); 
        w=w/length(x); 
        out=mean(x*y); 
        for i = 1:maxlag 
            out=out+w(i)*sum(trimr(x,i,0)*trimr(y,0,i))+w(i)*sum(trimr(y,i,0)*trimr(x,0,i));
        end  
end 

function out = NeweyWest_Var(H,g,maxlag)       

    if nargin < 3, maxlag = []; end
      
    I = 0*H;
    
    if isempty(maxlag); 
        n = size(g,1);
        maxlag = max(5,n^(0.25));
    end
              
    if maxlag > 0
        for i = 1:size(g,2)
            for j = i:ncol(g)
                I(i,j) = Newey_West(g(:,i),g(:,j),maxlag);
                I(j,i) = I(i,j);
            end    
        end
    end
   out =  H\(I + g'*g)/H;

end


