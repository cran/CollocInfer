function d2fdxdpval = loggenlin_d2fdxdp(t,y,p,more)
% LOGGENLIN_D2FDXDP computes second partial derivatives of values of fits  
%  to logged observations at time t with respect to state values in Y and
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
[Arow,Acol] = size(more.mat);
Amat  = zeros(Arow,Acol);
nAsub = size(more.sub,1);
for i=1:nAsub
    Amat(more.sub(i,1),more.sub(i,2)) = p(more.sub(i,3));
end
%  compute contribution to fit values from Y
fitval = y * Amat';
dfitdx = zeros(m,Arow,n);
for i = 1:Arow
    for j = 1:n
        dfitdx(:,i,j) = Amat(i,j);
    end
end
dfitdp = zeros(m,size(more.mat,1),length(p));
for i = 1:size(more.sub,1)
    dfitdp(:,more.sub(i,1),more.sub(i,3)) = y(:,more.sub(i,2));
end
%  compute contributions to fit derivatives with respect to P 
%     from each forcing vector
if ~isempty(more.force)
    %  set up matrix of forcing function values at time t
    fs = zeros(length(t),length(more.force));
    for i = 1:length(more.force)
        if is_fd(more.force{i})
            fs(:,i) = eval_fd(t,more.force{i});
        else
            fs(:,i) = more.force{i}(t,more.force.input);
        end
    end
    %  Add contributions to fit from forcing functions
    for i = 1:size(more.force.sub,1)
        dfitdp(:,more.force.sub(i,1),more.force.sub(i,3)) =  ...
        dfitdp(:,more.force.sub(i,1),more.force.sub(i,3)) +  ...
             fs(:,more.force.sub(i,2));
    end
end
%  compute contribution to fit derivatives from Y
d2fitval = zeros(m,Arow,n,npar);
psub = more.sub;
for i = 1:Arow;
    d2fitval(:,psub(i,1),psub(i,2),psub(i,3)) = 1;
end
d2fdxdpval = d2fitval - kron(dfitdp/fitval, dfitx/fitval);
%  contributions to fit derivatives from each forcing vector are zero

end

