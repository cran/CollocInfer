addpath('F:/work/External/matlab_packages/fdaM')
addpath('../SSE')
addpath('../exptrans')
addpath('../logtrans')
addpath('../findif')
addpath('../ChemoStat')

% The Chemostat equations represent a four-species Chemostat plus the resource 
% of Nitrogen. There are two species of Algae with varying defenses against
% Rotifers. The Rotifers themselves are divided into two class -- breeding
% and senescent, although these two are very tightly coupled. 
%
% A full description of these equations can be found in the user manual. 
% The five state variables for the equations are
% 
% N - nitrogen content in the Chemostat
% C1 - Algal type 1
% C2 - Algal type 2 
% B - Breeding Rotifers
% S - Senescent Rotifers
%
% The system has 16 parameters, also described in the user manual. Notable
% features include that only the sums C1+C2 and B+S can be observed. Further, 
% an unknown fraction of each is counted at each time. This requires us to 
% set up a model for the observation process along with the ODE. 


% First we generate some data 

fid  = fopen('ChemoData.txt','rt');
ChemoData = reshape(fscanf(fid,'%f'),2,101)';
fclose(fid);

% The first two of these parameters give the fractions of Algae and 
% Rotifers that are counted. The remaining parameters are all positive and  
% using their logged values is helpful. 

ChemoTime = 0:100;
ChemoPars = [1.000e+02 1.000e+00 5.000e-02 1.000e+00 1.600e+02 ...
             3.000e-01 5.500e-02 4.000e-01 2.025e-03 4.400e+00 ...
             2.200e+00 2.700e+02 1.100e-02 1.700e+02 1.500e-01 ...
             1.500e-01];
ChemoVarNames = ['N '; 'C1'; 'C2'; 'B '; 'S '];
ChemoParNames = ['a1    '; ...    
                 'a2    '; ...    
                 'p1    '; ...    
                 'p2    '; ...    
                 'NI    '; ... 
                 'delta '; ...         
                 'm     '; ...
                 'lambda'; ...
                 'XC    '; ...   
                 'KC1   '; ...    
                 'KC2   '; ...    
                 'rho   '; ...      
                 'G     '; ...   
                 'XB    '; ...    
                 'KB    '; ... 
                 'Qstar '];


% Parameters 'p1' and 'p2' represent relative palatability of the two algal
% clones, as such only one can be estimated and we fix p2 = 0. 

active = [1:2,5,7:16];     

% We'll choose a fairly large value of lambda. 

lambda = 10000*ones(5,1); 

% We need some basis functions

rr = [min(ChemoTime) max(ChemoTime)];
knots = rr(1):0.5:rr(2);

bbasis = create_bspline_basis(rr,[],4,knots);

% We now need to define two functions. First, chemo_fun defines a map from the
% trajectory to its derivative

fn = @chemo_fun_tmp;

% Second, we need a map from the trajectory to the observations. We'll use 
% log observations in this case since the observations are positive and this
% helps to stabilize scales. 

likfn = @chemo_obs;

%%% From here we will set-up the objects needed to run profiling. Below we give
% LS_setup the parameters, maps from trajectories to derivatives and to
% observations, along with data, times, basis functions, variable names and
% trade-off parameter lambda. 
%
% Additionally, posproc and poslik indicate that we are modelling the log of the
% trajectory, so that we should exponentiate before mapping to derivatives and
% before mapping to observations. 

% We'll start by defining some dummy coefficients:

coefs0 = zeros(getnbasis(bbasis),5);

[lik, proc, coefs] = LS_setup(fn,ChemoTime,coefs0,bbasis,lambda,[],[],[],[],1e-6,1,1,0,likfn,[]);

% Because we don't have direct observations of any state, we'll use a starting
% smooth obtained from generating some ODE solutions. 

y0 = log([2,0.1,0.4,0.2,0.1]);
logpars = [ChemoPars(1:2) log(ChemoPars(3:16))];

[odetrajt,odetrajy] = IntegrateForward(y0,ChemoTime,logpars,proc);

DEfd   = smooth_basis(odetrajt,odetrajy,fdPar(bbasis,int2Lfd(2),1e-6));
coefs0 = getcoef(DEfd);

obsvals = lik.more.fn(odetrajt,odetrajy,ChemoPars,lik.more.more);

plot(odetrajt,obsvals)
hold on
plot(ChemoTime,log(ChemoData))
hold off

ChemoData = exp( obsvals + 0.5*randn(size(obsvals)) );

%% We can now call the inner optimization -- this just tries to fined the 
% coefficients that both fit the data and the differential equation well. 

global INNEROPT_COEFS0;

INNEROPT_COEFS0 = coefs0;

control_in = optimset('LargeScale', 'off', 'GradObj', 'on', ...
                      'Hessian', 'off', 'Diagnostics', 'off', ...
                      'Display', 'iter');

coefs1 = inneropt(ChemoTime, log(ChemoData), logpars, lik, proc, [], control_in);

% We'll form the trajectory and also the appropriate sum of exponentiated
% states to compare to the data. 

traj    = lik.bvals * coefs1;
obstraj = lik.more.fn(ChemoTime,traj,logpars,lik.more.more);

% Plot these against the data

figure(1)
subplot(2,1,1)
plot(ChemoTime, obstraj(:,1), '-', ChemoTime, log(ChemoData(:,1)), 'bo' )
ylabel('Algae')
subplot(2,1,2)
plot(ChemoTime, obstraj(:,2), '-', ChemoTime, log(ChemoData(:,2)), 'bo')
ylabel('Brachionus')
xlabel('days')

% The function outeropt now allows parameter estimation to take place through 
% profiling procedure.

control_in = optimset('LargeScale', 'on', 'GradObj', 'on', ...
                      'Hessian', 'on', 'Diagnostics', 'off', ...
                      'Display', 'iter', ...
                      'MaxIter', 50);


control_out = optimset('LargeScale', 'off', 'GradObj', 'off', ...
                      'Hessian', 'off', 'Diagnostics', 'off', ...
                      'MaxIter', 0, ...
                      'Display', 'iter');

tic;
[pars_opt, coesfs2, value, gradient] = ...
         outeropt(ChemoTime, log(ChemoData), coefs1, logpars, ...
                  lik, proc, active, ...
                  [], control_in, [], control_out);
toc

% And obtain an estimated trajectory and the exponentiated sum to comprare
% to the data. 

traj = lik.bvals*coefs2;
ptraj = lik.more.fn(ChemoTime,traj,pars_opt,lik.more.more);

% Lets have a look at how much we changed our parameters on the original 
% scale. 

new_pars = pars_opt;
new_pars(3:16) = exp(new_pars(3:16));

disp(ChemoPars)
disp(new_pars)
disp(new_pars/ChemoPars)

% Now we can produce a set of diagnostic plots. 

% Firstly, a representation of the trajectory compared to the data. 

subplot(2,1,1)
plot(ChemoTime,ptraj(:,1))
ylabel('Chlamy')
xlabel('')
points(ChemoTime,ChemoData(:,1))
subplot(2,1,2)
plot(ChemoTime,ptraj(:,2))
ylabel('Brachionus')
xlabel('days')
points(ChemoTime,ChemoData(:,2))

% Now we'll plot both the derivative of the trajectory and the value of the
% differential equation right hand side at each point. This represents the 
% fit to the model. 

traj2 = proc.bvals.bvals*coefs2;
dtraj2 = proc.bvals.dbvals*coefs2;

colnames(traj2) = ChemoVarnames;
ftraj2 = proc.more.fn(proc2.more.qpts,traj2,pars_opt,proc.more.more);


for i = 1:5
  subplot(5,1,i)  
  plot(mids,dtraj2(:,i))
  xlabel('')
  ylabel(ChemoVarnames(i))
  lines(mids,ftraj2(:,i))
end
legend('Smooth','Model','location','NorthWest')

% Solving the differential equation from the estiamted initial condition
% of the trajectory allows us to compare the qualitative behavior of 
% our estimate to that of the differential equation. 

y0 = traj(1,:);
names(y0) = ChemoVarnames;
[odetrajt,odetrajy] = IntegrateForward(y0,ChemoTime,pars_opt,proc);


subplot(2,1,1)
plot(ChemoTime,traj)
ylabel('')
title('Reconstructed Trajectories')
legend(ChemoVarnames,'location','NorthEast')
subplot(2,1,2)
plot(ChemoTime,odetrajy)
ylabel('')
title('ODE Solution')

% We can also compare the pattern of observations predicted by the differential
% equation and that estimated by our methods. 

otraj = lik.more.fn(ChemoTime,odetraj.states,pars_opt,lik.more.more);

subplot(2,1,1)
plot(ChemoTime,ptraj)
xlabel('days')
ylabel('')
title('Predicted Observations -- Smooth')
plot(ChemoTime,ChemoData,'.')
legend(['Algae','Rotifers'],'location','NorthEast')
subplot(2,1,2)
plot(ChemoTime,otraj)
xlabel('days')
ylabel('')
title('Predicted Observations -- ODE')
legend(['Algae','Rotifers'],'location','NorthEast')
