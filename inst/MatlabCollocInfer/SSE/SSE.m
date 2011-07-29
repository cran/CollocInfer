%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fnval = SSE(data,times,devals,pars,more)

if ~isfield(more, 'more')
    more.more = [];
end
fdevals = more.fn(times, devals, pars, more.more);
%  Define the values of the process at observation times

difs = data - fdevals;
%  if process is unobserved, data(i) is NaN, and so is difs(i)
difs(isnan(difs)) = 0;
if isfield(more,'which')
    whichobs = more.which;
else
    whichobs = [];
end
weights = checkweights(more.weights,whichobs,difs);
fnval = sum(weights.*difs.^2,2);

end

