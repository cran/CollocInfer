function dfdpval = SEIR_dfdp(t,y,p,more)

t       = t(:);
beta    = more.beta_fun(t,p,more);
beta    = beta(:);
dbetadp = more.beta_dfdp(t,p,more);
tmpmat = diag(y(:,1).*(p(2)+y(:,3)))*dbetadp;
betaYS = beta.*y(:,1);
r = zeros(length(t),3,length(p));
r(:,1,more.beta_ind) = -tmpmat;
r(:,2,more.beta_ind) =  tmpmat;
r(:,1,1) = 1;
r(:,1,2) = -betaYS;
r(:,2,2) =  betaYS;
r(:,1,3) = -y(:,1);
r(:,2,3) = -y(:,2);
r(:,3,3) = -y(:,3);
r(:,2,4) = -y(:,2);
r(:,3,4) =  y(:,2);
r(:,3,5) = -y(:,3);
dfdpval  = r;

end

