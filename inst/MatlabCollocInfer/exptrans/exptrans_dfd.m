function x = exptrans_dfdp(times,y,p,more)

y = exp(y);
x = more.dfdp(times,y,p,more.more);

end
