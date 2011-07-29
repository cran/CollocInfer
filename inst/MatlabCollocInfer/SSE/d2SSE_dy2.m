%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dy2val = d2SSE_dy2(data,times,devals,pars,more)

[m,n] = size(data);
r = zeros(m,n,n);
ind = [repmat((1:m)',1,n), kron([1:n,1:n],ones(m,1))];
weights = checkweights(more.weights,more.whichobs,difs);
r(ind)  = weights;
dy2val = 2*r;

end

