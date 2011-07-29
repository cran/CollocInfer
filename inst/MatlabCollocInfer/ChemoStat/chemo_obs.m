function r = chemo_obs(t,y,p,more)

  Algae = p(1)*(y(:,2)+y(:,3));
  Rotifers = p(2)*(y(:,4)+y(:,5));
  
  r = log([Algae, Rotifers]);

end

