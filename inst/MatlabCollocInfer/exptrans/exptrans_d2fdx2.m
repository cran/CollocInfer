function x = exptrans_d2fdx2(times,y,p,more)

x1 = exptrans_dfdx(times,y,p,more);
y = exp(y);
x = more.d2fdx2(times,y,p,more.more);

for i = 1:size(x,2)
    for j = 1:size(x,3)
        for k = 1:size(x,4)
            x(:,i,j,k) = x(:,i,j,k).*y(:,j).*y(:,k);
        end
        x(:,i,j,j) = x(:,i,j,j) + x1(:,i,j);
    end
end

end