function fnval = loggenlin_fn(t,y,p,more)  
% LOGGENLIN_FN computes values of fits to logged observations 
%   at time t.
%  Observations are related to variables by the equation
%  Fit = Ay + Bf where y is the column vector of state values at time t,
%  and f is a matrix of forcing function values at time t.
%  Nonzero elements of A and B are taken from parameter vector P.

y = exp(y(:))';
if nargin < 4,  more = [];  end
%  check argument MORE
n    = size(y,2);
npar = length(p);
more = checkmore(more,n,npar);
%  Set up mapping matrix A
[Arow,Acol] = size(more.mat);
Amat  = zeros(Arow,Acol);
nAsub = size(more.sub,1);
for i=1:nAsub
    Amat(more.sub(i,1),more.sub(i,2)) = p(more.sub(i,3));
end
%  compute contribution to fit values from Y
fitval = y * Amat';
%  compute contributions to fit values from each forcing function vector
if  ~isempty(more.force)
    %  set up matrix of forcing function values at time t
    fs = zeros(length(t),length(more.force));
    for i = 1:length(more.force)
        if isa_fd(more.force{i})
            fs(:,i) = eval_fd(t,more.force{i});
        else
            fs(:,i) = more.force{i}(t,more.force_input);
        end
    end
    %  set up mapping matrix B
    [Brow,Bcol] = size(more.force.mat);
    Bmat  = zeros(Brow,Bcol);
    nBsub = size(more.force.sub,1);
    for i=1:nBsub
        Bmat(more.force.sub(i,1),more.force.sub(i,2)) = ...
            p(more.force.sub(:,3));
    end
    %  Add contributions to fit from forcing functions
    fitval = fitval + fs * Bmat';
end
disp(fitval)
fnval = log(fitval);

end
