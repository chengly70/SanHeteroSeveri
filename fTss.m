function ss=fTss(v)    % steady state for inactivation gating variable fT in T type Ca current
ss=1./(1+exp((v+58.7)./3.8));
end