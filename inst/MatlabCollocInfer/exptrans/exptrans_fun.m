%%%%%%%%%%% exponentiates-then-function composition  %%%%%%%%%%%%%%%%%%

function x = exptrans_fun(times,y,p,more)

y = exp(y);
x = more.fn(times,y,p,more.more);

end