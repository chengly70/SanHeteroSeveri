function ss=fLss(v)     % steady state for inactivation gating variable fL in L type Ca current
ss=1./(1+exp((v+37.4)./5.3));
end