function tau=taupaS(v)      % time constant for activation gating variable paS in Kr current
tau=0.84655./(4.2*exp(v./17)+0.15*exp(-v./21.6));
end