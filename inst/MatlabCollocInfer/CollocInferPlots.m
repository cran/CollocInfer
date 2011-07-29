function [timevec,traj,dtraj,ftraj,otraj] = CollocInferPlots(coefs,pars,lik,proc,times,data,...
                    ObsPlot,DerivPlot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CollocInferPlots -- diagnostic plots for CollocInfer
%
% These plots examine the goodness of fit of the results of profiling with
% CollocInfer. In particular, they plot the departure of the estimated
% trajectory from the data as well as the departure of the estimated
% trajectory from the differential equation. 
%
% Arguments:
% coefs            ... estimated coefficients
% pars             ... estimated parameters
% lik              ... lik object defining the likelihood of observations
% proc             ... proc object defining the ODE
% times            ... observation times
% data             ... observed data
% ObsPlot          ... indicator of whether the observation diagnostic plot
%                      should be given.
% DerivPlot        ... indicator of whether the ODE diagnostic plots should
%                      be given.
%
% Observation diagnostic is simply a plot of the observations along with
% the predicted observations at the trajectory defined by coefs. 
%
% ODE diagnostic plots compare the derivative of the estimated trajectory
% to the right hand side of the differential equation. These should be the
% same if the ODE is correct. This produces to plots
%     - An overlay of the predicted derivatives from the ODE (solid) along
%     with those taken from the trajectory (dashed). 
%     - A two-panel plot looking at the difference of the two derivatives
%     above. The lower panel plots the estimated trajectory as a means of
%     visually comparing if there appears to be a systematic relationship
%     between the differences in the derivative and the value of the state
%     variables. 

% Last updated: August 4, 2011. 

if nargin < 8, DerivPlot = 1; end
if nargin < 7, ObsPlot = 1; end

%   if isempty(cols) 
%       cols = 1:size(coefs,2); 
%   end

 timevec = proc.more.qpts;

 traj = proc.bvals.bvals*coefs;
 
 dtraj = proc.bvals.dbvals*coefs;

 ftraj = proc.more.fn(timevec,traj,pars,proc.more.more);
 
 otraj = lik.more.fn(timevec,traj,pars,lik.more.more);
 
 
 
 if(ObsPlot)
  if isempty(times),  times = 1:nrow(data); end
%   if isempty(datacols)  
%       if size(otraj,2) == size(traj,2),
%           datacols = cols;
%       else
%           datacols = 1:size(otraj,2);
%       end
%  end
  
  figure(1)
  plot(timevec,otraj);
  xlabel('time')
  ylabel('Observations')
  if isempty(data) == 0
%    if isempty(datanames), datanames = substring(colnames(data),1,1); end
    if isempty(times), times = 1:size(data,2); end
    hold on
    plot(times,data,'.')
    hold off
  end
 end
  
 if(DerivPlot)
  figure(2)
  plot(timevec,ftraj,'-')
  xlabel('time')
  ylabel('f(x): -, dx: --')
  
  hold on
  plot(timevec,dtraj,'--')
  plot([min(timevec) max(timevec)],[0 0],'k')
  hold off
%  legend(x='topright',legend=proc.more.names,lwd=2,lty=1,col=cols)
 
  figure(3)
  subplot(2,1,1)
  plot(timevec,dtraj-ftraj)
  xlabel('time')
  ylabel('dx-f(x)')
  hold on
  plot([min(timevec) max(timevec)],[0 0],'k')
  hold off
  
  subplot(2,1,2)
  plot(timevec,traj)
  xlabel('time')
  ylabel('x')
  
%  legend(x='topright',legend=proc.more.names,lwd=2,lty=1,col=cols)
 end
end