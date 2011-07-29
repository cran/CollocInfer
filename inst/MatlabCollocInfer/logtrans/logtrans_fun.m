function xoverexpy = logtrans_fun(times,y,p,more)

if ~isfield(more, 'more')
    more.more = [];
end
expy = exp(y);
x = more.fn(times,expy,p,more.more);
xoverexpy = x./expy;

end
