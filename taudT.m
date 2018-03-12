function tau=taudT(v)       % time constant for activation gating variable dT in T type Ca current
tau=0.001./(1.068*exp((v+38.3)./30)+1.068*exp(-(v+38.3)./30));
end