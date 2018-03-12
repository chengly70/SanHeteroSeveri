function tau=taufL(v)     % time constant for inactivation gating variable fL in L type Ca current
tau=0.001*(44.3+230*exp(-((v+36)./10).^2));
end