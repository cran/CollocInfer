function x = exptrans_dfdx(times,y,p,more)

y  = exp(y);
x  = more.dfdx(times,y,p,more.more);

for i = 1:size(x,2)
    for j = 1:size(x,3)
        x(:,i,j) = x(:,i,j).*y(:,j);
    end
end

end
