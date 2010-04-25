%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dyval = dSSE_dy(data,times,devals,pars,more)

if ~isfield(more,'more')
    more.more = [];
end
fdevals = more.fn(times,devals,pars,more.more);
difs    = data-fdevals;
difs(isnan(difs)) = 0;
weights = checkweights(more.weights,more.whichobs,difs);
difs    = weights.*difs;
dyval   = 2*difs;

end

