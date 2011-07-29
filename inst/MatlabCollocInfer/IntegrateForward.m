function [times,states] = IntegrateForward(y0,ts,pars,proc,more)
% INTEGRATEFORWARD solves a differential equation using the proc object
% provided to it. This function is written because CollocInfer requires 
% a slightly idfferent format that ode45 in order to solve differential
% equations
%
% Arugments:
%
% y0       ... initial conditions for the ODE
% ts       ... time points at which ODE solution is required
% pars     ... parameter values to go into the right hand side of the ODE
% proc     ... definition of the right hand side of the ODE. This can
%              either be a proc object as used by CollocInfer, or a right
%              hand side function that would be input into CollocInfer. 
% more     ... if proc is a function rather than a proc object, more is a 
%              list allowing additional inputs to proc. If a proc object is
%              given, more is ignored. 
%
% Ouput:
% times    ... list of times at which a solution to the ODE is given. 
% states   ... matrix of solution values at times. 

% Last modified 29 July, 2011. 


if nargin < 5, more = []; end

if strcmp(class(proc), 'function_handle')
    proc2.fn = proc;
    proc.more = proc2;
    proc.more.more = more;
end

    parms.pars = pars;
    parms.proc = proc;
    parms.more = more;
    
    [times,states] = ode45(@oderhs,ts,y0,[],parms);
    
end



