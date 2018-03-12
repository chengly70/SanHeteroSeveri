function tau=tauyf(v) % time constant for gating variable y in funny current
tau=0.7./(0.0708*exp(-(v+5)./20.28)+10.6*exp(v./18));
end