function ss=dTss(v)       % steady state for activation gating variable dT in T type Ca current
ss=1./(1+exp(-(v+38.3)./5.5));
end