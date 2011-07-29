function F = forward_prediction_error(times,data,coefs,lik,proc,pars,whichtimes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward Prediction Error is a means of estimating smoothing parameters
% for generalized profiling. It defines an objective function that
% integrates an exact solution of an ODE over a specified time window and
% comparing these forecasts to the observed data.  The solution has initial 
% conditions given by the value of the profile-estimated trajector at the 
% start of the window and uses the profile-estiamted parameters. Smoothing
% parameters can be chosen to minimize this criterion. 
%
% Arguments:
% times          ... vector of observation times
% data           ... matrix of observed data values with rows corresponding
%                    to entries of times. 
% coefs          ... estimated coefficients from profiling
% lik            ... lik object for the system
% proc           ... proc object for the system
% pars           ... parameters estimated from profiling
% whichtimes     ... indicates which time points to be used. This can have
%                    a number of formats.
%                    - cell array: each cell gives the indices of a set of
%                     observations to be used. y0 is chosen for the first
%                     index and ODE is integrated and compared to the
%                     observations listed in the subsequent indices. 
%                    - k-by-2 matrix of indices. y0 is chosen to be the 
%                      index in column 2 and the ODE is integrated to the
%                      data corresponding to the index in column 2. A
%                      comparison is made to all data points in between. 
%                    - empty or not given, we start at each data point and 
%                      go to the next. 
%
% Output:
%
% F - the forward prediction error: log likelihood (or squared error)
% summed over all predictions and observations. 

% Last modified August 4, 2011. 

  if nargin<7, whichtimes=[]; end
  
  if isempty(whichtimes) 
      whichtimes = [ (1:(length(times)-1))',(2:length(times))' ]; 
  end
  
  if isnumeric(whichtimes)
    twhich = whichtimes;
    whichtimes = cell(size(twhich,1),1);
    for i = 1:size(twhich,1)
      whichtimes{i} = twhich(i,1):twhich(i,2);
    end
  end

  traj = lik.bvals*coefs;
  ptraj = [];
  pdata = [];
  pts = [];
  whichobs = [];
  for i = 1:length(whichtimes) 

      temptimes = whichtimes{i};
      y0 = traj(temptimes(1),:);
      ts = times(temptimes);
 
      [odet,odey] = IntegrateForward(y0,ts,pars,proc);  
      if length(ts) == 2,
         odey = odey([1 size(odey,1)],:); 
      end
      
      temptimes = temptimes(2:length(temptimes));
      
      ptraj = [ptraj; odey(2:size(odey,1),:)];
      pdata = [pdata; data(temptimes,:)];
      pts  =  [pts; times(temptimes)];
      whichobs = [whichobs, temptimes];
  end

  lik.more.which = whichobs;
  
  F = sum(lik.fn(pdata,pts,ptraj,pars,lik.more));
end

