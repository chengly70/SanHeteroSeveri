function tau=taupaF(v)      % time constant for activation gating variable paF in Kr current
tau=1./(30*exp(v./10)+exp(-v./12));
end