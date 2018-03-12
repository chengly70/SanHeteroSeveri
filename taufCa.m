function tau=taufCa(casub)   % time constant for inactivation gating variable fCa in L type Ca current
afCa = 0.01; % m/s Ca dissociation rate constant for ICaL
tau=(0.001*fCass(casub))./afCa;
end