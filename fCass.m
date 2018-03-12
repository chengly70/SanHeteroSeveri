function ss=fCass(casub)       % steady state for inactivation gating variable fCa in L type Ca current
KmfCa = 0.00035; % mM Dissociation constant of Ca-dependent IcaL inactivation
ss=KmfCa./(KmfCa+casub);
end