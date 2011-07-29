function y = oderhs(t,r,parms)
% ODERHS makes a proc object to a function that can be employed by ode45 to
% solve differential equations. 
%
% Arguments:
% 
% (t,r)     ... time and state variable as used by functions in ode45. 
% parms     ... a struct containing elements
%                - proc: a proc object for SSE calculations
%                - pars: a set of parameter values
% 
% This format is used for compatibility with 'odesolve' in R. 

% Last modified 29 July, 2011

proc = parms.proc;
pars = parms.pars;

r = reshape(r,1,length(r));
    
y = squeeze(proc.more.fn(t,r,pars,proc.more.more));
y = reshape(y,length(y),1);

end