function ss=dLss(v)       % steady state for activation gating variable dL in L type Ca current
ss=1./(1+exp(-(v+20.3)./4.2));
end