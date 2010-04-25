function d2fdxdpval = genlin_d2fdxdp(t,y,p,more)
% GENLIN_D2FDXDP computes second partial derivatives of values of fits to 
%  observations at time t with respect to state values in Y and
%  parameter values in P
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
%  compute contribution to fit derivatives from Y
d2fdxdpval = zeros(m,Arow,n,npar);
psub = more.sub;
for i = 1:Arow;
    d2fdxdpval(:,psub(i,1),psub(i,2),psub(i,3)) = 1;
end
%  contributions to fit derivatives from each forcing vector are zero

end

