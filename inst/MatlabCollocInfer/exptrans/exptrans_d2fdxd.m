function x = exptrans_d2fdxd(times,y,p,more)

y = exp(y);
x = more.d2fdxdp(times,y,p,more.more);

for i = 1:size(x,2)
    for j = 1:size(x,3)
        for k = 1:size(x,4)
            x(:,i,j,k) = x(:,i,j,k).*y(:,j);
        end
    end
end

end
