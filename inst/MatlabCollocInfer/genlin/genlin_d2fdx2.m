function d2fdx2val = genlin_d2fdx2(t,y,p,more)
% GENLIN_D2FDX2 computes second partial derivatives of values of fits to 
%  observations at time t with respect to state values in X.
%  Observations are related to variables by the equation
%  Fit = Ay + Bf where y is the column vector of state values at time t,
%  and f is a matrix of forcing function values at time t.
%  Nonzero elements of A and B are taken from parameter vector P.

if nargin < 4,  more = []; end
%  check argument MORE
[m,n] = size(y);
npar  = length(p);
more  = checkmore(more,n,npar);
%  set up mapping matrix A
Arow = size(more.mat,1);
%  compute contribution to fit values from Y
d2fdx2val = zeros(m,Arow,n,n);
%  contributions to fit derivatives from each forcing vector are zero

end
