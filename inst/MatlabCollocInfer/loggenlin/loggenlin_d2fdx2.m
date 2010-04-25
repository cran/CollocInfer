function d2fdx2val = loggenlin_d2fdx2(t,y,p,more)
%  LOGGENLIN_D2FDX2 computes second partial derivatives of values of fits  
%  to logged observations at time t with respect to state values in X.
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
[Arow,Acol] = size(more.mat);
Amat  = zeros(Arow,Acol);
nAsub = size(more.sub,1);
for i=1:nAsub
    if more.sub(i,1) > Arow || ...
       more.sub(i,2) > Acol || ...
       more.sub(i,3) > npar
        error(['Subscripts in row ', num2str(i), ...
               ' of more.sub are inconsistent ', ...
               'with the size of more.mat.']);
    end
    Amat(more.sub(i,1),more.sub(1,2)) = p(more.sub(i,3));
end
%  compute contribution to fit values from Y
fitval = y * Amat';
dfitdx = zeros(m,Arow,n);
for i = 1:Arow
    for j = 1:n
        dfitdx(:,i,j) = Amat(i,j);
    end
end
%  compute contribution to fit values from Y
d2fdx2val = zeros(m,Arow,n,n) - kron(dfitdx/fitval, dfitdx'/fitval);
%  contributions to fit derivatives from each forcing vector are zero

end
