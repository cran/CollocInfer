function fnval = fhn_fun(times,y,p,more)

r = y;
r(:,1) = p(3).*(y(:,1) - (y(:,1).^3)./3 + y(:,2));
r(:,2) = -(y(:,1) - p(1) + p(2).*y(:,2))./p(3);
fnval = r;

end
