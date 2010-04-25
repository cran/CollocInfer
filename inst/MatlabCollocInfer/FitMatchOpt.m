function coefs = FitMatchOpt(coefs0, which, pars, proc, options_in)

if nargin < 5, options_in = [];  end

allcoefs   = coefs0;
coefs_strt = coefs0(:,which);
coefs_strt = coefs_strt(:);

if isempty(options_in)
    options_in = optimset('LargeScale', 'on', 'GradObj', 'on', ...
                          'Hessian', 'on', ...
                          'Display', 'iter', 'MaxIter', 100);
end

coefs_opt  = fminunc(@FitMatchCoefs, coefs_strt, options_in, ...
                     allcoefs, which, pars, proc); 
                 
npars     = length(coefs_opt);
n         = size(allcoefs,1);  
coefs_opt = reshape(coefs_opt,n,npars/n);
coefs     = coefs0;
coefs(:,which) = coefs_opt;

end


