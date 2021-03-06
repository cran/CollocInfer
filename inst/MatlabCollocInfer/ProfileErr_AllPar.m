function [value, gradient, ncoefs] = ...
    ProfileErr_AllPar(pars, times, data, coefs, lik, proc, ...
                      in_method, options_in, sumlik)
%  Inner optimization of coefficient values 
%
%  Arguments:
%    The first six arguments are required.  The remaining arguments each
%    have default values, can be omitted.  All arguments after the first
%    omitted argument must be also omitted.
%  PARS       ...  A vector of parameter values.  These are initial values
%                  that are used to start the outer iterations.  Some
%                  of these may be fixed, and the remainder are optimized.
%                  The indices of the optimized parameters are contained
%                  in argument ACTIVE.
%  TIMES      ...  The times of observation of those variables that are
%                  measured.  These times are assumed to be common to all
%                  measured variables.
%  DATA       ...  Data matrix.  The number of columns is equal to the
%                  number of variables (including unobserved variables).
%                  The number of columns is equal to the number of 
%                  observation points per variable times the number of
%                  replications.
%  COEFS0     ...  A matrix containing coefficients for the expansion.  
%                  It has NVAR columns, and its number of rows is
%                  NBASIS*NVAR.
%  LIK        ...  A struct object set up by function LS_SETUP that
%                  defines the loss function for fitting the data
%  PROC       ...  A struct object set up by function LS_SETUP that
%                  defines the loss function for fitting the 
%                  differential equation.
%  IN_METHOD  ...  A string identifing the optimization function to be
%                  used in the outer optimization loop.  It defaults to [].
%  OPTIONS_IN ...  A struct object containing values controlling the 
%                  inner optimization process.
%  SUMLIK     ...  Logical constant indicating whether or not the 
%                  likelihood is summed across times.
%                  Default is 1.

%  Last modified 13 January 2011

if nargin < 9, sumlik     = 1;   end
if nargin < 8, options_in = [];  end
if nargin < 7, in_method  = [];  end

global INNEROPT_COEFS0

INNEROPT_COEFS0 = coefs; 

nbasis = size(lik.bvals,2);

ncoefs = inneropt(times, data, pars, lik, proc, in_method, options_in);

devals = lik.bvals * ncoefs;

value = sum(lik.fn(data, times, devals, pars, lik.more));

if nargout > 1
    
[~, ~, d2fdc2, d2fdcdp] = ...
    SplineCoefs(coefs, times, data, pars, lik, proc);
dcdp = -d2fdc2\d2fdcdp;

if sumlik
    dlikdc = lik.bvals'*lik.dfdx(data, times, devals, pars, lik.more);
    dlikdc = dlikdc(:);
    df = dcdp'*dlikdc + ...
         sum(lik.dfdp(data, times, devals, pars, lik.more),1)';
    gradient = df(:);
else
    dlikdx = lik.dfdx(data, times, devals, pars, lik.more);    
    dlikdp = lik.dfdp(data, times, devals, pars, lik.more);    
    [m,n] = size(dlikdx);
    dlikdc = zeros(m,n*nbasis);
    m2 = 0;
    for i = 1:n
        m1 = m2 + 1;
        m2 = m2 + nbasis;
        dlikdc(:,m1:m2) = diag(dlikdx(:,i))*lik.bvals;
    end
    gradient = dlikdc * dcdp + dlikdp;
end

end

end
