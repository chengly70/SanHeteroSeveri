function tau=taudL(v)       % time constant for activation gating variable dL in L type Ca current
tau=0.001./(adL(v)+bdL(v));
end